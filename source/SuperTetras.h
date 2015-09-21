#ifndef SUPERTETRAS_H
#define SUPERTETRAS_H

#include "Mesh.h"
#include <cmath>
#include "assert.h"
#include <algorithm>
#include <set>
#include <list>
#include <cfloat>
#include "Eigen/Dense"
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/QR"
#include "Eigen/Sparse"

#if __cplusplus > 201103L
#include <unordered_set>
#include <unordered_map>
#else
#include <tr1/unordered_set>
#include <tr1/unordered_map>
#endif

using namespace Eigen;
using namespace std;

typedef unsigned long long int gridkey;

namespace LibTetra {


template<class Coordinate, class FieldValue> class SuperTetras: public Mesh<Coordinate, FieldValue>{

private:
    //fields
    FieldValue *clusterDistance;   // Distance of the tetra to the closest centroid
    vector<Vertex<Coordinate, FieldValue> > centroids;    // Centroids of the Tetrahedra
    Coordinate expectedRegSize;
    int nClusters;
    double weightOfField;
    vector<Vertex<Coordinate, FieldValue> > regCenters;
    vector<Triangle> boundaries;
    vector<Triangle> strongBs;
    Vertex<Coordinate, FieldValue> *mBarycenter;
    std::tr1::unordered_map<int, TetraIndex> oldStartPoints; // Centroid from where we started in the last iteration
    int stdNorm; //ADD command line option to read this value
    FieldValue spanField;
    bool addIter;

    //functions

    /**
     * @brief retrieveCentroids Calculates and store the barycenters of each tetrahedron
     */
    inline void retrieveCentroids(){

        for(TetraIndex TI=0; TI < this->tetraNumber(); TI=(TetraIndex)(TI+1)){
            Tetra* T = this->getTetra(TI);
            Coordinate accX=0.0, accY=0.0, accZ=0.0;
            FieldValue accF=0.0;

            for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                VertexIndex VI=T->TV(pi);
                Vertex<Coordinate, FieldValue> *vert = this->getVertex(VI);
                accX += vert->x();
                accY += vert->y();
                accZ += vert->z();
                accF += vert->fieldValue();
            }
            centroids.push_back(Vertex<Coordinate, FieldValue>(accX/4.0, accY/4.0, accZ/4.0, accF/4.0));
        }
    }

    /**
     * @brief containsCoincident checks if a tetrahedron is degenerate and contains two or more coincident vertices
     * @param v1 first vertex of the tetrahedron
     * @param v2 second vertex of the tetrahedron
     * @param v3 third vertex of the tetrahedron
     * @param v4 fourth vertex of the tetrahedron
     * @return true if at least two vertices coincide, false otherwise
     */
    inline bool containsCoincident(VertexIndex v1, VertexIndex v2, VertexIndex v3, VertexIndex v4){

        if(this->getVertex(v1) == this->getVertex(v2) || this->getVertex(v1) == this->getVertex(v3) || this->getVertex(v1) == this->getVertex(v4) || this->getVertex(v2) == this->getVertex(v3) || this->getVertex(v2) == this->getVertex(v4) || this->getVertex(v3) == this->getVertex(v4))
            return true;
        else
            return false;
    }

    /**
     * @brief triangleNormal calculates the normal vector to a triangle (normal stored in res)
     * @param v1 Vector which represent a vertex
     * @param v2 Vector which represent a vertex
     * @param v3 Vector which represent a vertex
     * @param res The normal to the triangle (normalized)
     */
    inline void triangleNormal(Matrix<Coordinate,3,1> v1, Matrix<Coordinate,3,1> v2, Matrix<Coordinate,3,1> v3, Matrix<Coordinate,3,1> &res){

        Matrix<Coordinate, 3,1> a,b;
        a=v2-v1;
        b=v3-v1;

        res=a.cross(b);
        res.normalize();
    }

    /**
     * @brief getPlaneEquation retrieves the equation of a plane
     * @param v1 vertex
     * @param v2 vertex
     * @param v3 vertex
     * @param plane equation of the plane
     */
    inline void getPlaneEquation(VertexIndex v1, VertexIndex v2, VertexIndex v3, Matrix<Coordinate, 4,1> &plane){

        Matrix<Coordinate, 3, 1> params(0.0,0.0,0.0);
        Matrix<Coordinate, 3, 1> vec1(this->getVertex(v1)->x(), this->getVertex(v1)->y(), this->getVertex(v1)->z());
        Matrix<Coordinate, 3, 1> vec2(this->getVertex(v2)->x(), this->getVertex(v2)->y(), this->getVertex(v2)->z());
        Matrix<Coordinate, 3, 1> vec3(this->getVertex(v3)->x(), this->getVertex(v3)->y(), this->getVertex(v3)->z());

        triangleNormal(vec1, vec2, vec3, params);

        for(int k=0; k<3; k++)
            plane[k]=params[k];
        plane[3] = -params[0]*this->getVertex(v1)->x() - params[1]*this->getVertex(v1)->y() - params[2]*this->getVertex(v1)->z();
    }

    /**
     * @brief planarPoints check if 4 vertices are coplanar
     * @param v1 index of first vertex
     * @param v2 index of second vertex
     * @param v3 index of third vertex
     * @param v4 index of fourth vertex
     * @return true if the points lie on the same plane, false otherwise
     */
    inline bool planarPoints(VertexIndex v1, VertexIndex v2, VertexIndex v3, VertexIndex v4){

        Matrix<Coordinate, 4, 1> planeEq(0.0,0.0,0.0,0.0);
        Coordinate r=0.0, res=0.0;

        getPlaneEquation(v1, v2, v3, planeEq);
        r = (this->getVertex(v4)->x()*planeEq[0] + this->getVertex(v4)->y()*planeEq[1] + this->getVertex(v4)->z()*planeEq[2]);
        res = planeEq[3] + r;
        if(fabs(res) < 10E-5)
            return true;

        getPlaneEquation(v1, v2, v4, planeEq);
        r = this->getVertex(v3)->x()*planeEq[0] + this->getVertex(v3)->y()*planeEq[1] + this->getVertex(v3)->z()*planeEq[2];
        res = planeEq[3] + r;
        if(fabs(res) < 10E-5)
            return true;

        getPlaneEquation(v1, v3, v4, planeEq);
        r=this->getVertex(v2)->x()*planeEq[0] + this->getVertex(v2)->y()*planeEq[1] + this->getVertex(v2)->z()*planeEq[2];
        res = planeEq[3] + r;
        if(fabs(res) < 10E-5)
            return true;

        getPlaneEquation(v2, v3, v4, planeEq);
        r=this->getVertex(v1)->x()*planeEq[0] + this->getVertex(v1)->y()*planeEq[1] + this->getVertex(v1)->z()*planeEq[2];
        res = planeEq[3] + r;
        if(fabs(res) < 10E-5)
            return true;

        return false;
    }

    /**
     * @brief tetraUnsignedVolume calculates the (unsigned) volume of a tetrahedron
     * @param TI index of a tetrahedron
     * @return the volume
     */
    inline Coordinate tetraUnsignedVolume(TetraIndex TI){

        Tetra* T = this->getTetra(TI);

        if(containsCoincident(T->TV(a), T->TV(b), T->TV(c), T->TV(d)))
            return 0.0;
        if(planarPoints(T->TV(a), T->TV(b), T->TV(c), T->TV(d)))
            return 0.0;

        // not a degenerate case, computing volume
        Coordinate vol;
        Matrix<Coordinate, 3, 1> u,v,z,n;

        u = {this->getVertex(T->TV(b))->x()-this->getVertex(T->TV(a))->x(), this->getVertex(T->TV(b))->y()-this->getVertex(T->TV(a))->y(), this->getVertex(T->TV(b))->z()-this->getVertex(T->TV(a))->z()};
        v = {this->getVertex(T->TV(c))->x()-this->getVertex(T->TV(b))->x(), this->getVertex(T->TV(c))->y()-this->getVertex(T->TV(b))->y(), this->getVertex(T->TV(c))->z()-this->getVertex(T->TV(b))->z()};
        z = {this->getVertex(T->TV(d))->x()-this->getVertex(T->TV(c))->x(), this->getVertex(T->TV(d))->y()-this->getVertex(T->TV(c))->y(), this->getVertex(T->TV(d))->z()-this->getVertex(T->TV(c))->z()};

        n = u.cross(v);
        vol = n[0]*z[0] + n[1]*z[1] + n[2]*z[2];
        return fabs(vol/6.0);
    }

    /**
     * @brief meshVolume
     * @return the volume of the whole mesh
     */
    inline Coordinate meshVolume(){

        Coordinate accV=0.0;
        for(TetraIndex TI=0; TI < this->tetraNumber(); TI=(TetraIndex)(TI+1)){
            accV = accV + tetraUnsignedVolume(TI);
        }
        cout<<"Acc "<<accV<<endl;
        return accV;
    }

    inline FieldValue getValMinField(){

        FieldValue fld = this->getVertex(0)->fieldValue();

        for(int ii=0; ii<this->vertexNumber(); ii++){
            if(fld > this->getVertex(ii)->fieldValue())
                fld = this->getVertex(ii)->fieldValue();
        }
        return fld;
    }

    inline FieldValue getValMaxField(){

        FieldValue fld = this->getVertex(0)->fieldValue();

        for(int ii=0; ii<this->vertexNumber(); ii++){
            if(fld < this->getVertex(ii)->fieldValue())
                fld = this->getVertex(ii)->fieldValue();
        }
        return fld;
    }


    /**
     * @brief avgField
     * @return the average of the field values on the vertices
     */
    inline FieldValue avgField(){

        FieldValue acc=0.0;
        for(VertexIndex VI=0; VI<this->vertexNumber(); VI=(VertexIndex)(VI+1))
            acc=acc+this->getVertex(VI)->fieldValue();
        return acc/this->vertexNumber();
    }

    /**
     * @brief stdField
     * @return the standard deviation of the field value of the mesh vertices
     */
    inline FieldValue stdField(){

        FieldValue avgF = avgField();
        FieldValue acc=0.0;

        for(VertexIndex VI=0; VI<this->vertexNumber(); VI=(VertexIndex)(VI+1)){
            acc = acc + (this->getVertex(VI)->fieldValue()-avgF)*(this->getVertex(VI)->fieldValue()-avgF);
        }
        return sqrt(acc/(this->vertexNumber()-1));
    }

    /**
     * @brief rangeFieldNorm performs a range normalization of the scalar field (takes it in [0,1])
     */
    inline void rangeFieldNorm(){

        FieldValue mn = getValMinField();
        FieldValue mx = getValMaxField();

        for(VertexIndex VI=0; VI<this->vertexNumber(); VI=(VertexIndex)(VI+1)){
            this->getVertex(VI)->setFieldValue((this->getVertex(VI)->fieldValue()-mn)/(mx-mn));
        }
    }

    /**
     * @brief stdDevFieldNorm performs a standard deviation based normalization of the scalar field
     */
    inline void stdDevFieldNorm(){

        FieldValue stdf = stdField();
        cout<<"STDF "<<stdf<<endl;
        for(VertexIndex VI=0; VI<this->vertexNumber(); VI=(VertexIndex)(VI+1))
            this->getVertex(VI)->setFieldValue(this->getVertex(VI)->fieldValue()/stdf);
    }

    /**
     * @brief fieldNormalization normalizes the scalar field
     * @param option whether we want a range (0) or a std-dev based normalization (1)
     */
    inline void fieldNormalization(int option){

        if(option==0)
            rangeFieldNorm();
        else
            stdDevFieldNorm();
    }

    inline void setSpanField(){
        spanField = (getValMaxField()-getValMinField());
        cout<<"Span "<<spanField<<endl;
    }

    /**
     * @brief combinedDistance calculates the distance between tetrahedra (used in classification step)
     * @param SD Euclidean term
     * @param FD Field value based term
     * @return the distance between two tetrahedra
     */
    inline FieldValue combinedDistance(Coordinate SD, FieldValue FD){

        Coordinate denominator = expectedRegSize*expectedRegSize;
        return FD + weightOfField*weightOfField*(SD/denominator);
    }

    /**
     * @brief getGridKey get the key to access a grid cell
     * @param i
     * @param j
     * @param k
     * @return
     */
    inline gridkey getGridKey(short int i, short int j, short int k){

        return i << 16 | j << 8 | k;
    }

    inline gridkey getMapKey(gridkey i, gridkey j){

        if( i<= j)
            return i << 32 | j;
        return j << 32 | i;
    }

    /**
     * @brief meshBarycenter calculates the barycenter of the mesh
     */
    inline void meshBarycenter(){

        Coordinate xx=0.0, yy=0.0, zz=0.0;
        for(VertexIndex VI=0; VI < this->vertexNumber(); VI=(VertexIndex)(VI+1)){

            xx = xx + this->getVertex(VI)->x();
            yy = yy + this->getVertex(VI)->y();
            zz = zz + this->getVertex(VI)->z();
        }

        mBarycenter = Vertex<Coordinate, FieldValue>(xx/this->vertexNumber(), yy/this->vertexNumber(), zz/this->vertexNumber(), 0.0);
    }

    /**
     * @brief closestToMeshBarycenter
     * @return the index of the tetrahedron closest to the mesh barycenter
     */
    inline TetraIndex closestToMeshBarycenter(){

        TetraIndex winner=0;
        Coordinate distwin = mBarycenter->EuclideanDistance3D(centroids.at(0));

        for(TetraIndex TI=1; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){

            if(distwin > mBarycenter->EuclideanDistance3D(centroids.at(TI))){
                winner = TI;
                distwin = mBarycenter->EuclideanDistance3D(centroids.at(TI));
            }
        }
        return winner;
    }

    /**
     * @brief initSegmentation performs an initial subdivision of the mesh into a regular grid
     */
    inline void initSegmentation(){

        tr1::unordered_map<gridkey, int> hMap;
        int id = 0;

        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){

            short int i = floor((centroids.at((int)TI).x()-this->lowerLeft->x())/(expectedRegSize));
            short int j = floor((centroids.at((int)TI).y()-this->lowerLeft->y())/(expectedRegSize));
            short int k = floor((centroids.at((int)TI).z()-this->lowerLeft->z())/(expectedRegSize));

            //cout<<"Tetraindex "<<TI<<" Key "<<((centroids.at((int)TI).x()-this->lowerLeft->x())/(2*expectedRegSize))<<" "<<((centroids.at((int)TI).y()-this->lowerLeft->y())/(2*expectedRegSize))<<" "<<((centroids.at((int)TI).z()-this->lowerLeft->z())/(2*expectedRegSize))<<endl;

            if(hMap.count(getGridKey(i,j,k)) == 0){

                regCenters.push_back(centroids.at((int)TI));
                hMap[getGridKey(i,j,k)] = id++;
                oldStartPoints[hMap[getGridKey(i,j,k)]] = TI;
            }

            this->clusterIndex[TI] = hMap[getGridKey(i,j,k)];
        }
    }

    /**
     * @brief updateSegmentation after a classification step, takes the centroid of each region to the approximate center of the region
     * @return the average shift of the centroids
     */
    inline Coordinate updateSegmentation(){

        vector< Vertex<Coordinate, FieldValue> > auxRegCentroids;
        Coordinate accVol[regCenters.size()];
        Coordinate accCoords[regCenters.size()][3];
        FieldValue accField[regCenters.size()];
        Coordinate totError=0.0;

        for(unsigned int ii=0; ii<regCenters.size(); ii++){
            auxRegCentroids.push_back(regCenters.at(ii));
            accVol[ii] = 0.0;
            accField[ii] = 0.0;
            for(int j=0; j<3; j++)
                accCoords[ii][j]=0.0;
        }
        regCenters.clear();

        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){

            Coordinate tVol = this->tetraUnsignedVolume(TI);
            int reg = this->clusterIndex[(int)TI];
            accVol[reg] = accVol[reg] + tVol;
            accCoords[reg][0] = accCoords[reg][0] + tVol*centroids.at((int)TI).x();
            accCoords[reg][1] = accCoords[reg][1] + tVol*centroids.at((int)TI).y();
            accCoords[reg][2] = accCoords[reg][2] + tVol*centroids.at((int)TI).z();
            accField[reg] = accField[reg] + tVol*centroids.at((int)TI).fieldValue();
        }

        for(unsigned int ii=0; ii<auxRegCentroids.size(); ii++){

            if(accVol[ii]>0){
                Vertex<Coordinate, FieldValue> V = Vertex<Coordinate, FieldValue> (accCoords[ii][0]/accVol[ii], accCoords[ii][1]/accVol[ii], accCoords[ii][2]/accVol[ii], accField[ii]/accVol[ii]);
                regCenters.push_back(V);
                totError = totError + V.EuclideanDistance3D(auxRegCentroids.at(ii));
            }
            else{
                regCenters.push_back(auxRegCentroids.at(ii));
            }
        }
        nClusters = regCenters.size();
        cout<<"Total error "<<totError/nClusters<<", regions "<<nClusters<<endl;
        return totError/nClusters;
    }

    inline void buildMaptoAvgField(){

        std::tr1::unordered_map<int, FieldValue> indexToValue;

        for(int ii=0; ii<nClusters; ii++){
            indexToValue[ii] = avgSegmColor(ii);//regCenters.at(ii).fieldValue();
        }

        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){
            int indx = this->clusterIndex[(int)TI];
            this->clusterFieldValue[(int)(TI)] = indexToValue[indx];
        }
    }

    inline FieldValue avgSegmColor(int regIndex){

        FieldValue fv = 0.0;
        int counter=0;
        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){
            if(this->clusterIndex[(int)TI]==regIndex){
                counter++;
                fv=fv+centroids.at(TI).fieldValue();
            }

        }
        FieldValue toRet = fv/counter;
        return toRet;
    }

    /**
     * @brief neighClosest
     * @param TI index of a tetrahedron (old center)
     * @param regInd index of the region to which the tetrahedron belongs
     * @return The index of the tetrahedron closest to TI
     */
    inline TetraIndex neighClosest(TetraIndex TI, int regInd){

        Coordinate dd = regCenters.at(regInd).EuclideanDistance3D(centroids.at(TI));
        TetraIndex toret = TI;

        for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
            TetraIndex tn = this->getTetra(TI)->TT(pi);
            if(tn >= 0 && tn < this->tetraNumber() && regCenters.at(regInd).EuclideanDistance3D(centroids.at(tn)) < dd && this->clusterIndex[tn]==regInd){
                dd = regCenters.at(regInd).EuclideanDistance3D(centroids.at(tn));
                toret = tn;
            }
        }
        return toret;
    }

    inline TetraIndex closestToRC(int regionIndex, int &steps){

        TetraIndex Tindex = oldStartPoints[regionIndex];
        TetraIndex prev = Tindex;

        TetraIndex winner = neighClosest(Tindex, regionIndex);
        int countStep=0;
        while(winner != prev){
            prev=winner;
            winner = neighClosest(prev, regionIndex);
            countStep++;
        }
        steps=countStep;
        return winner;
    }

    inline TetraIndex bruteForceClosest(int regionIndex){

        TetraIndex Tindex = oldStartPoints[regionIndex];
        //TetraIndex winner=Tindex;
        FieldValue distance = DBL_MAX;//regCenters.at(Tindex).EuclideanDistance3D()
        for(TetraIndex TI=0;TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){

            if(this->clusterIndex[(int)TI]==regionIndex){
                if(distance > regCenters.at(regionIndex).EuclideanDistance3D(centroids.at((int)TI))){
                    Tindex=TI;
                    distance = regCenters.at(regionIndex).EuclideanDistance3D(centroids.at((int)TI));
                }
            }
        }
        return Tindex;
    }

    /**
     * @brief maxOf3
     * @param a
     * @param b
     * @param c
     * @return the maximum of three values
     */
    inline Coordinate maxOf3(Coordinate a, Coordinate b, Coordinate c){

        if(a > b && a > c)
            return a;
        //a is not the bigest
        if(b > c)
            return b;
        //if we're here, then it is c
        return c;
    }

    /**
     * @brief expandRegion Classification step for a region
     * @param NC Actual centroid of the region
     * @param region Index of the region
     */
    inline void expandRegion(TetraIndex NC, int region){

        //TetraIndex NC = closestToRC(region);
        if(region >= regCenters.size()){
            regCenters.push_back(centroids.at((int)NC));
            cout<<"Reg "<<region<<" size "<<regCenters.size()<<endl;
        }
        oldStartPoints[region] = NC;

        std::tr1::unordered_set<TetraIndex> visited;
        queue<TetraIndex> Q;
        Q.push(NC);

        while(!Q.empty()){
            TetraIndex first = Q.front();
            Q.pop();

            if(visited.count(first))
                continue;

            Coordinate ED = centroids.at(NC).EuclideanDistance3D(centroids.at(first));
            ED = ED*ED;
            FieldValue FD = (centroids.at(NC).fieldValue()-centroids.at(first).fieldValue())*(centroids.at(NC).fieldValue()-centroids.at(first).fieldValue());
            FieldValue dist = combinedDistance(ED, FD);

            if(dist < clusterDistance[first]){
                clusterDistance[first] = dist;
                this->clusterIndex[first] = region;
            }
            visited.insert(first);
            Coordinate latDist = maxOf3(fabs(centroids.at(first).x()-regCenters.at(region).x()), fabs(centroids.at(first).y()-regCenters.at(region).y()), fabs(centroids.at(first).z()-regCenters.at(region).z()));
            if(latDist < 2*expectedRegSize){
                for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                    if(this->getTetra(first)->TT(pi) >= 0 && this->getTetra(first)->TT(pi) < this->tetraNumber())
                        Q.push(this->getTetra(first)->TT(pi));
                }
            }
        }
    }

    /**
     * @brief classificationStep Assigns each tetrahedron to the closest centroid
     */
    inline void classificationStep(){

        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1))
            clusterDistance[(int)TI] = DBL_MAX;
        int accumStep = 0, numberS=0, countDiff=0;

        TetraIndex startPoints[regCenters.size()];
        for(unsigned int ii=0; ii< regCenters.size(); ii++){
            startPoints[ii] = bruteForceClosest(ii);//closestToRC(ii, numberS);
            //TetraIndex auxsp = bruteForceClosest(ii);

            //            if(startPoints[ii] != auxsp){
            //                countDiff++;
            //                //cout<<"Reg "<<ii<<" st "<<startPoints[ii]<<" ("<<this->clusterIndex[startPoints[ii]]<<"), stbis "<<auxsp<<" ("<<this->clusterIndex[auxsp]<<")"<<endl;
            //            }

            accumStep+=numberS;
        }
        for(unsigned int ii=0; ii < regCenters.size(); ii++){
            expandRegion(startPoints[ii],ii);
        }
        nClusters=regCenters.size();
    }

    /**
     * @brief postProcess enforce connectivity of the regions
     */
    inline void postProcess(){

        int* auxClindex = new int[this->tetraNumber()];
        for(int ii=0; ii<this->tetraNumber(); ii++){
            auxClindex[ii]=this->clusterIndex[ii];
            this->clusterIndex[ii]=-1;
        }
        std::tr1::unordered_set<TetraIndex> enqueued;
        for(int ii=0; ii<nClusters; ii++){

            TetraIndex cnt = oldStartPoints[ii];

            queue<TetraIndex> TQ;
            TQ.push(cnt);
            enqueued.insert(cnt);

            while(!TQ.empty()){

                TetraIndex Tindex = TQ.front();
                TQ.pop();
                this->clusterIndex[Tindex] = ii;

                for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){

                    if(this->getTetra(Tindex)->TT(pi) >= 0 && this->getTetra(Tindex)->TT(pi) < this->tetraNumber()){
                        if(this->clusterIndex[this->getTetra(Tindex)->TT(pi)]<0 && auxClindex[this->getTetra(Tindex)->TT(pi)] == ii && !enqueued.count(this->getTetra(Tindex)->TT(pi))){
                            TQ.push(this->getTetra(Tindex)->TT(pi));
                            enqueued.insert(this->getTetra(Tindex)->TT(pi));
                        }
                    }
                }
            }
        }

        enqueued.clear();
        queue<TetraIndex> orphaned;
        for(TetraIndex TI=0; TI<this->tetraNumber(); TI=(TetraIndex)(TI+1)){
            if(this->clusterIndex[(int)TI] < 0)
                orphaned.push(TI);
        }

        while(!orphaned.empty()){

            TetraIndex head = orphaned.front();
            orphaned.pop();

            int winnerReg = -1;
            FieldValue windist = DBL_MAX;
            for(PositionIndex PI=a; PI<=d; PI=(PositionIndex)(PI+1)){

                if(this->getTetra(head)->TT(PI) >= 0 && this->getTetra(head)->TT(PI) < this->tetraNumber() && this->clusterIndex[this->getTetra(head)->TT(PI)]>=0){
                    Coordinate ED = centroids.at(head).EuclideanDistance3D(centroids.at(this->getTetra(head)->TT(PI)));
                    ED=ED*ED;
                    FieldValue FD = (centroids.at(head).fieldValue()-centroids.at(this->getTetra(head)->TT(PI)).fieldValue())*(centroids.at(head).fieldValue()-centroids.at(this->getTetra(head)->TT(PI)).fieldValue());
                    FieldValue dist = combinedDistance(ED, FD);
                    if(dist < windist){
                        windist = dist;
                        winnerReg = this->clusterIndex[this->getTetra(head)->TT(PI)];
                    }
                }
            }
            if(winnerReg >= 0)
                this->clusterIndex[head] = winnerReg;
            else
                orphaned.push(head);
        }
    }

    inline void getSegmBoundaries(){

        for(TetraIndex TI=0; TI < this->tetraNumber(); TI=(TetraIndex)(TI+1)){

            int indexCurrent = this->clusterIndex[(int)TI];
            for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                TetraIndex NI = this->getTetra(TI)->TT(pi);
                if(NI>=0 && NI<this->tetraNumber()){
                    int indexN = this->clusterIndex[NI];
                    if(/*!visitedT.count(NI) && */indexCurrent != indexN)
                        boundaries.push_back(this->getTetra(TI)->TF(pi));
                }
            }
        }
    }

public:

    SuperTetras(){ //constructor
    }

    /**
         * @brief initStructs
         * Initializes the structures needed for the segmentation process
         */
    inline void initStructs(){

        this->clusterIndex = new int[this->tetraNumber()];
        this->clusterFieldValue = new FieldValue[this->tetraNumber()];
        clusterDistance = new FieldValue[this->tetraNumber()];
        for(int t=0; t<this->tetraNumber(); t++){
            this->clusterIndex[t] = -1;   // At the beginning, they are all uninitialized
            //clusterDistance[t] = -1;
        }

        this->setUpperRight();
        this->setLowerLeft();
        this->getDiagonal();
        cout<<"Volume "<<meshVolume()<<endl;
        expectedRegSize = pow(meshVolume()/(nClusters), 1.0/3.0);
        fieldNormalization(stdNorm);
        retrieveCentroids();
    }

    /**
     * @brief iteration performs the segmentation
     */
    inline void iteration(){
        int counter=0;
        Coordinate diffMove=updateSegmentation();
        Coordinate prevError = diffMove;
        Coordinate newerror;
        while(diffMove>0.05/*this->BBDiagonal */&& counter < 50){
            classificationStep();
            newerror=updateSegmentation();
            diffMove = fabs(prevError-newerror);
            cout<<"Difference "<<diffMove<<" threshold "<<0.05<<endl;
            prevError = newerror;
            counter++;
        }
    }

    /**
     * @brief triggerSegmentation initialize and updates the segmentation
     */
    inline void triggerSegmentation(){
        initSegmentation();
        setSpanField();
        iteration();

        cout<<"At the end of the process, regions are "<<nClusters<<" ("<<regCenters.size()<<")"<<endl;
        postProcess();
        buildMaptoAvgField();
        getSegmBoundaries();

    }

    /**
     * @brief setNormStrategy how do we want the scalar field to be normalized
     * @param choice (0 for range norm, 1 for std-dev based)
     */
    inline void setNormStrategy(int choice){
        stdNorm=choice;
    }

    inline Vertex<Coordinate, FieldValue> getNthCentroid(int index){
        return centroids.at(index);
    }

    inline void setNDesiredRegions(int NR){
        this->nClusters = NR;
    }

    inline int getNDesiredRegions(){
        return this->nClusters;
    }

    inline void setWeightOfField(double WoF){
        this->weightOfField = WoF;
    }


    inline double getWeightOfField(){
        return this->weightOfField;
    }

    inline int countReg(){
        return this->nClusters;
    }


    inline void writeSupertetrasBoundaries(char *filename){

        FILE* file = fopen(filename, "w");
        if(file==NULL){
            cout<<"NULL "<<endl;
            return;
        }

        fprintf(file, "# vtk DataFile Version 2.0\n\n");
        fprintf(file, "ASCII\n");
        fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
        fprintf(file, "POINTS %d float\n", this->vertexNumber());

        for(int i=0; i<this->vertexNumber(); i++)
            fprintf(file, "%f %f %f\n", this->getVertex(i)->x(), this->getVertex(i)->y(), this->getVertex(i)->z());
        fprintf(file, "\n\n");

        fprintf(file, "CELLS %d %d\n", boundaries.size(), 4*boundaries.size());
        for(int i=0; i<boundaries.size(); i++){
            //Triangle t=boundaries.at(i);
            fprintf(file, "3 %d %d %d\n", boundaries.at(i).FV(a), boundaries.at(i).FV(b), boundaries.at(i).FV(c));
        }
        fprintf(file, "\n");

        fprintf(file, "CELL_TYPES %d\n", boundaries.size());
        for(int i=0; i<boundaries.size(); i++)
            fprintf(file, "%d\n", 5);
        fprintf(file, "\n\n");

        fclose(file);
        //boundaries.clear();
    }

    /**
     * @brief output_segmentation writes the segmentation indices to a file
     * @param filename name of the file to be written
     * @param sz
     */
    inline void output_segmentation(const char* filename, int sz){
        FILE* file;
        file = fopen(filename, "w");

        for(int i=0; i<this->tetraNumber(); i++)
            fprintf(file, "%d\n", this->clusterIndex[i]);

        fclose(file);
    }

    inline void output_segAvgField(const char* filename){

        FILE* file=fopen(filename, "w");

        if(file==NULL){
            cout<<"Invalid file pointer"<<endl;
            return;
        }

        for(int ii=0; ii<regCenters.size(); ii++)
            fprintf(file, "%f %f %f %f\n", regCenters.at(ii).x(), regCenters.at(ii).y(), regCenters.at(ii).z(), regCenters.at(ii).fieldValue());

        fclose(file);
    }
};
}


#endif // SUPERTETRAS_H
