#ifndef MESH_H
#define MESH_H
#include <vector>
#include <map>
#include <queue>
#include "Tetra.h"

namespace LibTetra{

    template<class Coordinate, class FieldValue>
    // Class storing the whole tetrahedral mesh, along with operations on it
    class Mesh{
    protected:
        //bool printingForTesting;
        std::vector<Vertex<Coordinate,FieldValue>* >* vertexes; // std::vector storing vertices, with indices coherent with those stored in tetrahedra
        std::vector<Tetra*>* tetras; // std::vector storing tetrahedra

        Coordinate BBDiagonal;  //Diagonal of the mesh bounding box
        Vertex<Coordinate, FieldValue>* lowerLeft;
        Vertex<Coordinate, FieldValue>* upperRight;


        vector<int> critical_vertices;
        int *clusterIndex;  // Index of the region to which the tetra belongs
        FieldValue *clusterFieldValue;
        //Functions

        // Index sorting function
        void sortIndexes(PositionIndex *indexes, int numberOfIndexes, Tetra* tetra);

        // Affinity planarizing the given triple of vertices
        inline mc1::matrix<Coordinate> getPlanarizer(Triangle triangle); //TODO

        // Turn of a std::vector v with respect to oriented edge with tail edgeExtremeA and head egedeExtremeB calculated over the first 2 coordinates of vertices
        inline Turn getXYTurn(mc1::vector<Coordinate> v,mc1::vector<Coordinate> edgeExtremeA,mc1::vector<Coordinate> edgeExtremeB); // TODO

        // Turn of a vertex with respect to the plane by three given points (stored in face)
        inline Turn getTurn(Vertex<Coordinate,FieldValue>* vertex, Triangle face); //TODO

        // Turn of a vertex over a plane with respecto to a line on it, identified by two given distinct vertices
        inline Turn getTurn(Vertex<Coordinate,FieldValue>* vertex, Triangle face, PositionIndex edgePositionOnFace); //TODO

        // Relative position of a vertex with respect to a tetrahedron
        TVposition getRelativePosition(Vertex<Coordinate,FieldValue>* vertex, Tetra* tetra, void** storeForIncidentFaceOrEdge);

        inline void pushVertex(Vertex<Coordinate,FieldValue>* vertex){
            vertexes->push_back(vertex);
        }

        // Adding a tetrahedron to mesh
        inline void pushTetra(Tetra* tetra){
            tetras->push_back(tetra);
        }

        // Sorts lexicographically triples of indices
        inline void sortFaces(vector<pair<Triangle,TetraIndex>* >* faces){
            sortFaces(faces,0,faces->size());
        }

        // Sorts lexicographically subsets of a list of triples of indices
        void sortFaces(vector<pair<Triangle,TetraIndex>* >* faces, int left, int right);

        // Merges sorted sublists of triples of indices
        void merge(vector<pair<Triangle,TetraIndex>* >* faces, int left, int center, int right);

        // Comparing function for triples of indices
        inline int minusTriangle(Triangle t1, Triangle t2){
            if(t1.minindex()< t2.minindex()) return -1;
            else if(t1.minindex() > t2.minindex()) return 1;

            if(t1.middleindex()< t2.middleindex()) return -1;
            else if(t1.middleindex() > t2.middleindex()) return 1;

            if(t1.maxindex()< t2.maxindex()) return -1;
            else if(t1.maxindex() > t2.maxindex()) return 1;

            return 0;
        }

        // Star of an edge
        std::vector<TetraIndex> ET(Edge edge,TetraIndex oneAdjacentTetra);

    public:
        // Vertex reference getter
        inline Vertex<Coordinate,FieldValue>* getVertex(VertexIndex index){return (vertexes->at(index));}

        // Tetrahedron reference getter
        inline Tetra* getTetra(TetraIndex index){if(index==noAdjacentTetraIndex)return NULL; else return (tetras->at(index));}

        // Constructor
        inline Mesh(){
        }

        // Number of stored vertices
        inline VertexIndex vertexNumber(){return vertexes->size();}

        // Number of stored tetrrahedra
        inline TetraIndex tetraNumber(){return tetras->size();}

        //
        std::vector<TetraIndex> FT(Triangle face);
        //
        std::vector<Triangle> FF(Triangle face);

        // Star of an edge
        std::vector<TetraIndex> ET(Edge edge);
        //
        std::vector<Triangle> EF(Edge edge);//TODO
        //
        std::vector<Edge> EE(Edge edge);//TODO
        //

        // Star of a vertex
        std::vector<TetraIndex> VT(VertexIndex center,bool* foundBorder);
        //
        std::vector<Triangle> VF(VertexIndex center);
        //
        std::vector<Edge> VE(VertexIndex center);
        //
        std::vector<VertexIndex> VV(VertexIndex center);
        //

        // Link of a vertex
        std::vector<TetraIndex> Lnk(Edge edge,TetraIndex oneAdjacentTetra); //TODO

        // Loader for TS tetrahedral mesh files
        Exception loadTS(FILE* TSfile);
        // Matrix building function
        Exception build();
        // Consistent point insertion function
        Exception pointInsertion(Vertex<Coordinate,FieldValue>* vertex);
        // Mesh coherence test
        void testCoherence();
        // Standard output operator
        inline friend ostream & operator << (ostream & out, Mesh & mesh){
            out<<"Mesh["<<endl;
            for(VertexIndex i=0; i<mesh.vertexNumber(); i++){
                out<<i<<": "<<*mesh.getVertex(i)<<endl;
            }
            for(TetraIndex i=0; i<mesh.tetraNumber(); i++){
                out<<i<<": "<<*mesh.getTetra(i)<<endl;
            }
            out<<"]";
            return out;
        }


        bool fixIncoherence(Incoherence incoherence);

        //Utility Function
        template<class C> bool isIntoVector(C c, vector<C> &c_vect);

        inline void setLowerLeft(){

            Coordinate minimum[3] = {getVertex(0)->x(), getVertex(0)->y(), getVertex(0)->z()};

            for(VertexIndex VI = 1; VI < vertexNumber(); VI=(VertexIndex)(VI+1)){

                if(getVertex(VI)->x() < minimum[0])
                    minimum[0] = getVertex(VI)->x();
                if(getVertex(VI)->y() < minimum[1])
                    minimum[1] = getVertex(VI)->y();
                if(getVertex(VI)->z() < minimum[2])
                    minimum[2] = getVertex(VI)->z();
            }
            lowerLeft = new Vertex<Coordinate, FieldValue>(minimum[0], minimum[1], minimum[2], 0.0);
        }

        inline void setUpperRight(){

            Coordinate maximum[3] = {getVertex(0)->x(), getVertex(0)->y(), getVertex(0)->z()};

            for(VertexIndex VI = 1; VI < vertexNumber(); VI=(VertexIndex)(VI+1)){

                if(getVertex(VI)->x() > maximum[0])
                    maximum[0] = getVertex(VI)->x();
                if(getVertex(VI)->y() > maximum[1])
                    maximum[1] = getVertex(VI)->y();
                if(getVertex(VI)->z() > maximum[2])
                    maximum[2] = getVertex(VI)->z();
            }
            upperRight = new Vertex<Coordinate, FieldValue>(maximum[0], maximum[1], maximum[2], 0.0);
        }

        inline void getDiagonal(){

            BBDiagonal = sqrt((upperRight->x()-lowerLeft->x())*(upperRight->x()-lowerLeft->x()) + (upperRight->y()-lowerLeft->y())*(upperRight->y()-lowerLeft->y()) + (upperRight->z()-lowerLeft->z())*(upperRight->z()-lowerLeft->z()));
        }


        vector<Coordinate> readField(char* file_name);

        ///IO functions
        void writeTS(char* file_name, vector<Coordinate> fieldvalues); ///write a new mesh in .ts format with the original scalar field switched with fieldvalues
        void writeField(char* file_name, vector<Coordinate> fieldvalues); ///write a file with a column of values taken from fieldvalues (first value is the number of values)
        void writeVTK(char* file_name, vector<Coordinate> fieldvalues); ///write a new mesh in .vtk format with the original scalar field switched with fieldvalues (for paraview)
        void writeVTK_segmentation(char* filename, vector<int>* segm1, vector<int>* segm2);
        void writeSuperTetrasVTK(char* filename);
        void writeIndicesVTK(char* filename);
        //void writeSupertetrasBoundaries(char* filename);
        //void output_segmentation(const char* filename, int sz);
        void output_segAvgField(char* filename, int sz);

    };
}

#include "Mesh.cpp"

#endif // MESH_H
