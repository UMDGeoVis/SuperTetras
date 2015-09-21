#ifndef LIBTETRA_Mesh
#define LIBTETRA_Mesh
#include "Mesh.h"

namespace LibTetra{
    template<class Coordinate, class FieldValue>
    bool Mesh<Coordinate,FieldValue> :: fixIncoherence(Incoherence incoherence){
        int fixed=0;
        switch(incoherence){
        case missingVTstars:{
                std::vector<VertexIndex> toBeDeleted;
                for(VertexIndex vind=0; vind<vertexNumber(); vind++){
                    if(getVertex(vind)->VTstar()==uninitializedTetraIndex) toBeDeleted.push_back(vind);
                }
                if(toBeDeleted.size()==0) return true;
                sort(toBeDeleted.begin(),toBeDeleted.end());
                for(TetraIndex tind=0; tind<tetraNumber(); tind++){
                    Tetra* tp=getTetra(tind);
                    for(PositionIndex pind=a; pind<=d; pind=(PositionIndex)(pind+1)){
                        VertexIndex toScale=0;
                        VertexIndex tvi=tp->TV(pind);
                        while(tvi>toBeDeleted[toScale]&&toScale<toBeDeleted.size())toScale++;
                        if(tvi==toBeDeleted[toScale]){
                            //getVertex(toBeDeleted[toScale])->setVTstar(tind);
                            cerr<<"A VTstar was not set: VertexIndex "<<toBeDeleted[toScale]<<" into TetraIndex "<<tind<<endl;
                            //return false;
                        }
                        if(tvi<toScale){
                            cerr<<"Trying to set a negative vertexIndex while fixing coherence: this CANNOT happen...";
                        }
                        tp->setTV(pind,tvi-toScale);
                    }
                }
                std::vector<Vertex<Coordinate,FieldValue>* >* newVertexes=new std::vector<Vertex<Coordinate,FieldValue>* >();
                int del=0;
                for(VertexIndex vind=0; vind<vertexNumber(); vind++){
                    Vertex<Coordinate,FieldValue>* vp=getVertex(vind);
                    //cerr<<vind<<" "<<toBeDeleted[del]<<" "<<del<<endl;
                    if(vind==toBeDeleted[del]){
                        del++;
                        delete vp;
                    }
                    else{
                        newVertexes->push_back(vp);
                    }
                }
                std::vector<Vertex<Coordinate,FieldValue>* >* tmp=vertexes;
                vertexes=newVertexes;
                delete tmp;
            };
            break;
        }
        return true;
    }
    
    template<class Coordinate, class FieldValue>
    void Mesh<Coordinate,FieldValue> :: testCoherence(){
#define outOfBounds_TetraIndex(A) ( A !=noAdjacentTetraIndex && ( A<0 || A>=tetraNumber() ) )
#define outOfBounds_VertexIndex(A) ( A<0 || A>=vertexNumber() )
#define outOfBounds_PositionIndex(A) ( A !=unknown && ( A<a || A>d ) )
        int topologyIncoherences=0;
        int structureIncoherences=0;
        TetraIndex tnum=tetraNumber();
        VertexIndex vnum=vertexNumber();
        cerr<<"Beginning coherence check..."<<endl;
        for(VertexIndex vindex=0; vindex<vnum; vindex++){
            if(outOfBounds_TetraIndex(getVertex(vindex)->VTstar())){
                cerr<<"- Structural Incoherence: VTstar of vertex "<<vindex<<" is invalid because out of bounds: "<<getVertex(vindex)->VTstar()<<endl;
                structureIncoherences++;
                continue;    
            }
            if(!getTetra(getVertex(vindex)->VTstar())->contains(vindex)){
                cerr<<"- Topological Incoherence: tetra "<<getVertex(vindex)->VTstar()<<" does not contain vertex "<<vindex<<" although being its VTstar."<<endl;
                //cout<<*getTetra(getVertex(vindex)->VTstar());
                //cerr<<endl;
                topologyIncoherences++;
            }
        }
        for(TetraIndex tindex=0; tindex<tnum; tindex++){
            Tetra* tetra=getTetra(tindex);
            if(tetra==NULL){
                cerr<<"- Structural Incoherence: tetra "<<tindex<<" is NULL."<<endl;
                structureIncoherences++;
                continue;
            }
            for(PositionIndex pindex=a; pindex<=d; pindex=(PositionIndex)(pindex+1)){
                if(outOfBounds_TetraIndex(tetra->TT(pindex))){
                    cerr<<"- Structural Incoherence: tetra "<<tindex<<" TT at "<<pindex<<" is invalid because out of bounds."<<endl;
                    structureIncoherences++;
                    continue;
                }
                if(outOfBounds_VertexIndex(tetra->TV(pindex))){
                    cerr<<"- Structural Incoherence: tetra "<<tindex<<" TV at "<<pindex<<" is invalid because out of bounds."<<endl;
                    structureIncoherences++;
                    continue;
                }
                if(outOfBounds_PositionIndex(tetra->inverseTT(pindex))){
                    cerr<<"- Structural Incoherence: tetra "<<tindex<<" inverseTT at "<<pindex<<" is invalid because out of bounds."<<endl;
                    structureIncoherences++;
                    continue;
                }
                
                Tetra* TTtetra=getTetra(tetra->TT(pindex));
                if(TTtetra==NULL){
                    continue;
                }
                VertexIndex TTindex=tetra->TV(pindex);
                if(TTtetra->contains(TTindex)){
                    cerr<<"- Topological Incoherence: tetra "<<tetra->TT(pindex)<<" contains vertex "<<tetra->TV(pindex)<<" it's opposite to with respect to tetra "<<tindex<<"."<<endl;
                    topologyIncoherences++;
                }
                Triangle face=tetra->TF(pindex);
                if(!TTtetra->contains(face)){
                    cerr<<"- Topological Incoherence: tetra "<<tetra->TT(pindex)<<" does not contain face shared with tetra "<<tindex<<" at position "<<pindex<<" on the latter."<<endl;
                    topologyIncoherences++;
                }
                if(face.contains(tetra->TV(pindex))){
                    cerr<<"- Topological Incoherence: tetra "<<tetra->TT(pindex)<<" contains vertex "<<tetra->TV(pindex)<<" at least twice."<<endl;
                    topologyIncoherences++;
                }
                
                if(TTtetra->TT(tetra->inverseTT(pindex))!=tindex){
                    cerr<<"- Topological Incoherence: inverseTT relation at position "<<pindex<<" of tetra "<<tindex<<" is incoherent with respect to adjacent tetra "<<tetra->TT(pindex)<<"."<<endl;
                    topologyIncoherences++;
                }
                
            }
        }
        cerr<<"Total incoherences on structure: "<<structureIncoherences<<endl;
        cerr<<"Total incoherences on topology: "<<topologyIncoherences<<endl;
        
#undef outOfBounds_TetraIndex
#undef outOfBounds_VertexIndex
#undef outOfBounds_PositionIndex
    }
    
    template<class Coordinate, class FieldValue>
    std::vector<TetraIndex> Mesh<Coordinate,FieldValue> :: ET(Edge edge, TetraIndex oneAdjacentTetrahedron){
        std::vector<TetraIndex> rtets;
        std::vector<TetraIndex> ltets;
        queue<TetraIndex> rqu;
        queue<TetraIndex> lqu;
        
        Tetra* tetra=getTetra(oneAdjacentTetrahedron);
//#ifdef LIBTETRA_DEBUG
//        cerr<<"  Beginning ET generation for edge "<<edge<<" of tetra "<<*tetra<<endl;
//#endif
        tetra->setStatus(checked);
        
        PositionIndex indexes[2];
        int filled=0;
        
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    looking for east and west pole indexes"<<endl;
//#endif
        for(PositionIndex pind=a; pind<=d; pind=(PositionIndex)(pind+1)){
            if(tetra->TV(pind)!=edge.minindex() && tetra->TV(pind)!=edge.maxindex()) indexes[filled++]=pind;
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    initializing queues"<<endl;
//#endif
        rqu.push(tetra->TT(indexes[0]));
        lqu.push(tetra->TT(indexes[1]));
        
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    processing first queue"<<endl;
//#endif
        while(!rqu.empty()){
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      beginning queue cycle..."<<endl;
//#endif
            TetraIndex index=rqu.front();
            rqu.pop();
            if(index==noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      no adjacent tetra found..."<<endl;
//#endif
                continue;
            }
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      processing tetra "<<index<<endl;
//#endif
            tetra=getTetra(index);
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      pushing tetra into right chain"<<endl;
//#endif
            if(tetra->status()==unchecked)rtets.push_back(index);
            tetra->setStatus(checked);
            filled=0;
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      looking for east and west position indexes"<<endl;
//#endif
            for(PositionIndex pind=a; pind<=d; pind=(PositionIndex)(pind+1)){
                if(tetra->TV(pind)!=edge.minindex() && tetra->TV(pind)!=edge.maxindex()){
                    indexes[filled++]=pind;
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      "<<pind<<" is east or west position index"<<endl;
//#endif
                }
            }
            for(int i=0; i<2; i++){
                TetraIndex newtetra=tetra->TT(indexes[i]);
                if(newtetra==noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      no adjacent tetra to push"<<endl;
//#endif
                    continue;
                }
                if(getTetra(newtetra)->status()==unchecked){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      pushing tetra "<<newtetra<<" into right queue"<<endl;
//#endif
                    rqu.push(newtetra);
                }
            }
        }
        
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    processing second queue"<<endl;
//#endif
        while(!lqu.empty()){
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      beginning queue cycle..."<<endl;
//#endif
            TetraIndex index=lqu.front();
            lqu.pop();
            if(index==noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      no adjacent tetra found..."<<endl;
//#endif
                continue;
            }
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      processing tetra "<<index<<endl;
//#endif
            tetra=getTetra(index);
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      pushing tetra into left chain"<<endl;
//#endif
            if(tetra->status()==unchecked)ltets.push_back(index);
            tetra->setStatus(checked);
            filled=0;
//#ifdef LIBTETRA_DEBUG
//            cerr<<"      looking for east and west position indexes"<<endl;
//#endif
            for(PositionIndex pind=a; pind<=d; pind=(PositionIndex)(pind+1)){
                if(tetra->TV(pind)!=edge.minindex() && tetra->TV(pind)!=edge.maxindex()){
                    indexes[filled++]=pind;
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      "<<pind<<" is east or west position index"<<endl;
//#endif
                }
            }
            for(int i=0; i<2; i++){
                TetraIndex newtetra=tetra->TT(indexes[i]);
                if(newtetra==noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      no adjacent tetra to push"<<endl;
//#endif
                    continue;
                }
                if(getTetra(newtetra)->status()==unchecked){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      pushing tetra "<<newtetra<<" into right queue"<<endl;
//#endif
                    lqu.push(newtetra);
                }
            }
        }
        
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    creating edge star chain"<<endl;
//        cerr<<"    ";
//#endif
        std::vector<TetraIndex> tets;
        for(int i=ltets.size()-1; i>=0; i--){
            TetraIndex index=ltets[i];
//#ifdef LIBTETRA_DEBUG
//            cerr<<index<<" < ";
//#endif
            getTetra(index)->setStatus(unchecked);
            tets.push_back(index);
        }
        if(rtets.size()!=0){
//#ifdef LIBTETRA_DEBUG
//            cerr<<oneAdjacentTetrahedron<<" < ";
//#endif
        }
        else{ 
//#ifdef LIBTETRA_DEBUG
//            cerr<<oneAdjacentTetrahedron<<" |";
//#endif
        }
        tetra->setStatus(unchecked);
        tets.push_back(oneAdjacentTetrahedron);
        for(int i=0; i<rtets.size(); i++){
            TetraIndex index=rtets[i];
            if(i!=rtets.size()-1){ 
//#ifdef LIBTETRA_DEBUG
//                cerr<<index<<" < ";
//#endif
            }
            else{ 
//#ifdef LIBTETRA_DEBUG
//                cerr<<index<<" |";
//#endif
            }
            getTetra(index)->setStatus(unchecked);
            tets.push_back(index);
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    cleaning"<<endl;
//#endif
        
//#ifdef LIBTETRA_DEBUG
//        cerr<<"  ET generation ended."<<endl;
//#endif
        return tets;
    }
    
    template<class Coordinate, class FieldValue>
    void Mesh<Coordinate,FieldValue> :: sortIndexes(PositionIndex *indexes, int numberOfIndexes, Tetra* tetra){
        for(int i=0; i<numberOfIndexes; i++){
            int minimum=i;
            for(int j=i+1; j<numberOfIndexes; j++){
                if(tetra->TV(indexes[j])<tetra->TV(indexes[minimum])) minimum=j;
            }
            PositionIndex tmp=indexes[i];
            indexes[i]=indexes[minimum];
            indexes[minimum]=tmp;
        }
    }

    template<class Coordinate, class FieldValue>
    Exception Mesh<Coordinate,FieldValue> :: build(){
        unsigned int VTstarSet=0;
        unsigned int TTandInverseTTSet=0;
        std::vector<pair<Triangle,TetraIndex>* > faces(tetraNumber()*4);
        
        for(TetraIndex index=0; index<tetraNumber(); index=index+1){
            Tetra* tetra=getTetra(index);
//#ifdef LIBTETRA_DEBUG
//            cerr<<"Tetra number "<<index<<":"<<endl;
//#endif
            
            for(PositionIndex vp=a; vp<=d; vp=(PositionIndex)(vp+1)){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Vertex number "<<vp<<"["<<tetra->TV(vp)<<"]:"<<endl;
//#endif
                Vertex<Coordinate,FieldValue>* vertex=getVertex(tetra->TV(vp));
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    VT*: "<<vertex->VTstar()<<endl;
//#endif
                if(vertex->VTstar()==uninitializedTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    Setting VT*..."<<endl;
//#endif
                    VTstarSet++;
                    vertex->setVTstar(index);
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    VT* set:"<<vertex->VTstar()<<endl;
//#endif
                }
                faces[tetraNumber()*vp+index]=new pair<Triangle,TetraIndex>(tetra->TF(vp),index);
            }
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"Sorting..."<<endl;
//#endif
        sortFaces(&faces);
//#ifdef LIBTETRA_DEBUG
//        cerr<<"Beginning TT and InverseTT build ..."<<endl;
//#endif
        for(unsigned int i=0; i<faces.size()-1; i++){
            Triangle firstFace=faces[i]->first;
            Triangle secondFace=faces[i+1]->first;
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  Comparing "<<firstFace<<" with "<<secondFace<<endl;
//#endif
            if((firstFace)==(secondFace)){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"   They are equal."<<endl;
//#endif
                TetraIndex firstIndex=faces[i]->second;
                TetraIndex secondIndex=faces[i+1]->second;
//#ifdef LIBTETRA_DEBUG
//                cerr<<"   Indexes are respectively "<<firstIndex<<" and "<<secondIndex<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  SettingRelations..."<<endl;
//#endif
                Tetra* firstTetra=getTetra(firstIndex);
                Tetra* secondTetra=getTetra(secondIndex);
                firstTetra->setTT(firstFace.faceIndex(),secondIndex);
                firstTetra->setInverseTT(firstFace.faceIndex(),secondFace.faceIndex());
                secondTetra->setTT(secondFace.faceIndex(),firstIndex);
                secondTetra->setInverseTT(secondFace.faceIndex(),firstFace.faceIndex());
                TTandInverseTTSet+=2;
                
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Deleting pairs..."<<endl;
//#endif
                delete faces[i];
                delete faces[i+1];
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Deleting triangles..."<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  ...Cleanup complete."<<endl;
//#endif
                i++;
            }
            else{
//#ifdef LIBTETRA_DEBUG
//                cerr<<"   They are NOT equal."<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Deleting triangle..."<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Deleting pair..."<<endl;
//#endif
                delete faces[i];
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  ...Cleanup complete."<<endl;
//#endif
            }
        }
        
        if(VTstarSet!=vertexNumber()) return someVTstarCouldNotBeSet;
        if(TTandInverseTTSet!=4*tetraNumber()) return someTTorInverseTTCouldNotBeSet;
//#ifdef LIBTETRA_DEBUG
//        cerr<<"build complete"<<endl;
//#endif
        return noException;
    }
    
    template<class Coordinate, class FieldValue>
    void Mesh<Coordinate, FieldValue> :: merge(vector<pair<Triangle,TetraIndex>* >* faces, int left, int center, int right){

//#ifdef LIBTETRA_DEBUG
//        cerr<<"    Merge step step: left="<<left<<" right="<<right<<endl;
//#endif
        int size=right-left;
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    creating auxiliary std::vector..."<<endl;
//#endif
        std::vector<pair<Triangle,TetraIndex>* > aux(size);
        int i,j;
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    interleaving..."<<endl;
//#endif
        for(i=left,j=center; i<center&&j<right;){
            int compare=minusTriangle(faces->at(i)->first,faces->at(j)->first);
            if(compare>0){
                aux[i+j-center-left]=faces->at(j);
                j++;
            }
            else{
                aux[i+j-center-left]=faces->at(i);
                i++;
            }
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    adding remains..."<<endl;
//#endif
        for(;i<center;i++){
            aux[i+j-center-left]=faces->at(i);
        }
        for(;j<right;j++){
            aux[i+j-center-left]=faces->at(j);
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"    copying std::vector..."<<endl;
//#endif
        for(int k=0; k<size; k++){
            faces->at(k+left)=aux[k];
        }

//#ifdef LIBTETRA_DEBUG
//        cerr<<"    deleting auxiliary..."<<endl;
//#endif

//#ifdef LIBTETRA_DEBUG
//        cerr<<"    ended."<<endl;
//#endif
    }
    
    template<class Coordinate, class FieldValue>
    void Mesh<Coordinate, FieldValue> :: sortFaces(vector<pair<Triangle,TetraIndex>* >* faces, int left, int right){
//#ifdef LIBTETRA_DEBUG
//        cerr<<"  Sort step: left="<<left<<" right="<<right<<endl;
//#endif
        if (right-left>2) {
            int center = (left+right)/2;
            sortFaces(faces, left, center);
            sortFaces(faces, center, right);
            merge(faces, left, center, right);
        }
        else{
            int compare=minusTriangle(faces->at(left)->first,faces->at(right-1)->first);
            if(compare>0){
                pair<Triangle,TetraIndex>* tmp=faces->at(right-1);
                faces->at(right-1)=faces->at(left);
                faces->at(left)=tmp;
            }
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"  ended."<<endl;
//#endif
    }
    template<class Coordinate, class FieldValue>
    Exception Mesh<Coordinate,FieldValue> :: loadTS(FILE* TSfile){
        int vNumber;
        int tNumber;
        fscanf(TSfile,"%d %d ",&vNumber,&tNumber);
        vertexes=new std::vector<Vertex<Coordinate,FieldValue>* >(vNumber);
        tetras=new std::vector<Tetra*>(tNumber);
        for(VertexIndex i=0; i<vertexNumber(); i=i+1){
            Coordinate x;
            Coordinate y;
            Coordinate z;
            FieldValue field;
            fscanf(TSfile,"%lf %lf %lf %lf ",&x,&y,&z,&field);
            vertexes->at(i)=new Vertex<Coordinate,FieldValue>(x,y,z,field);
        }
        for(TetraIndex i=0; i<tetraNumber(); i++){
            VertexIndex a,b,c,d;
            fscanf(TSfile,"%ld %ld %ld %ld ",&a,&b,&c,&d);

            Tetra* t = new Tetra(a,b,c,d);
            tetras->at(i)= t;
        }
        fclose(TSfile);
        return noException;
    }
    
    template<class Coordinate, class FieldValue>
    Exception Mesh<Coordinate,FieldValue> :: pointInsertion(Vertex<Coordinate,FieldValue>* vertex){
        
        TVposition relativePosition=external;
        TetraIndex tetraIndex=0;
        void * faceOrEdgeVertexIsOn;
//#ifdef LIBTETRA_DEBUG
//        cerr<<"Point Insertion: "<<*vertex<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//        cerr<<"  Getting relative position..."<<endl;
//#endif
        while(relativePosition==external && tetraIndex<tetraNumber()){
            relativePosition=getRelativePosition(vertex,getTetra(tetraIndex++),&faceOrEdgeVertexIsOn);
        }
        tetraIndex=tetraIndex-1;
        if(relativePosition==internal){
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  Point is internal to tetrahedron "<<tetraIndex<<endl;
//#endif
            VertexIndex vnum=vertexNumber();
            VertexIndex vertexIndex=vnum;
            TetraIndex tnum=tetraNumber();
            
            Tetra* tetra=getTetra(tetraIndex);
            vertex->setVTstar(tetraIndex);
            pushVertex(vertex);
            
            Triangle face;
            Tetra* newTetras[4];
#define finalIndex(A) ( A==a ? tetraIndex : tnum+A-1 )
            std::vector<pair<Triangle,TetraIndex>* > faces(16);
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  Creating new Tetras..."<<endl;
//#endif
            for(PositionIndex vp=a; vp<=d; vp=(PositionIndex)(vp+1)){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    number"<<vp<<": "<<endl;
//#endif
                face=tetra->TF(vp);
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      creating..."<<endl;
//#endif
                newTetras[vp]=new Tetra(vertexIndex,face.minindex(),face.middleindex(),face.maxindex());
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      setting TT..."<<endl;
//#endif
                newTetras[vp]->setTT(a,tetra->TT(vp));
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      setting inverse TT..."<<endl;
//#endif
                newTetras[vp]->setInverseTT(a,tetra->inverseTT(vp));                
                if(tetra->TT(vp)!=noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      getting adjacent tetra and setting TT and inverse TT..."<<endl;
//#endif
                    Tetra* adjTetra=getTetra(tetra->TT(vp));
                    adjTetra->setTT(tetra->inverseTT(vp),finalIndex(vp));
                    adjTetra->setInverseTT(tetra->inverseTT(vp),a);
                }
//#ifdef LIBTETRA_DEBUG
//                cerr<<"      getting faces..."<<endl;
//#endif
                for(PositionIndex fp=a; fp<=d; fp=(PositionIndex)(fp+1)){
                    faces[4*vp+fp]=new pair<Triangle,TetraIndex>(newTetras[vp]->TF(fp),finalIndex(vp));
                }
            }
            tetras->at(tetraIndex)=newTetras[a];
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  ...new Tetras created and added"<<endl;
//#endif
            for(PositionIndex vp=b; vp<=d; vp=(PositionIndex)(vp+1)) pushTetra(newTetras[vp]);
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  sorting faces..."<<endl;
//#endif
            sortFaces(&faces);
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  setting relationships..."<<endl;
//#endif
            for(int i=0; i<faces.size()-1; i++){
                Triangle firstFace=faces[i]->first;
                Triangle secondFace=faces[i+1]->first;
                if((firstFace)==(secondFace)){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    coincidcent faces: "<<firstFace<<" == "<<secondFace<<endl;
//#endif
                    TetraIndex firstIndex=faces[i]->second;
                    TetraIndex secondIndex=faces[i+1]->second;
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    getting tetras..."<<endl;
//#endif
                    Tetra* firstTetra=getTetra(firstIndex);
                    Tetra* secondTetra=getTetra(secondIndex);
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    setting TTs and inverse TTs..."<<endl;
//#endif
                    firstTetra->setTT(firstFace.faceIndex(),secondIndex);
                    firstTetra->setInverseTT(firstFace.faceIndex(),secondFace.faceIndex());
                    secondTetra->setTT(secondFace.faceIndex(),firstIndex);
                    secondTetra->setInverseTT(secondFace.faceIndex(),firstFace.faceIndex());
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    cleaning up..."<<endl;
//#endif
                    delete faces[i];
                    delete faces[i+1];
                    i++;
                }
                else{
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    cleaning up..."<<endl;
//#endif
                    delete faces[i];
                }
            }
//#ifdef LIBTETRA_DEBUG
//            cerr<<"    fixing VTstar of old vertex..."<<endl;
//#endif
            getVertex(tetra->TV(a))->setVTstar(finalIndex(1));
#undef finalIndex
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  cleaning up..."<<endl;
//#endif
            delete tetra;
            return noException;
        }
        if(relativePosition==face){
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  Point is internal to face "<<*(Triangle*)faceOrEdgeVertexIsOn<<" of tetrahedron "<<tetraIndex<<endl;
//#endif
            VertexIndex vnum=vertexNumber();
            VertexIndex vertexIndex=vnum;
            TetraIndex tnum=tetraNumber();
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  retrieving incident tetra indexes..."<<endl;
//#endif
            
            PositionIndex incidentFaceIndexs[2];
            TetraIndex tetraIndexes[2];
            Tetra* tets[2];
            
            incidentFaceIndexs[0]=((Triangle*)faceOrEdgeVertexIsOn)->faceIndex();
            tets[0]=getTetra(tetraIndex);
            tetraIndexes[0]=tetraIndex;
            TetraIndex otherIndex=tets[0]->TT(incidentFaceIndexs[0]);
            vertex->setVTstar(tetraIndex);
            pushVertex(vertex);
            
            if(otherIndex!=noAdjacentTetraIndex){
                tets[1]=getTetra(otherIndex);
                tetraIndexes[1]=otherIndex;
                incidentFaceIndexs[1]=tets[0]->inverseTT(incidentFaceIndexs[0]);
                for(int i=0; i<4; i++) pushTetra(NULL);
                
                Triangle face;
                Tetra* newTetras[3][2];
                
#define finalIndex(A,B) ( A==0 ? tetraIndexes[B] : tnum+3*B+A-B-1 )

//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Creating new Tetras..."<<endl;
//#endif
                for(int tetr=0; tetr<2; tetr++){
                    int actualIndex=0;
                    PositionIndex sortedIndexes[3];
                    int filled=0;
                    for(PositionIndex ind=a; ind<=d; ind=(PositionIndex)(ind+1)){
                        if(ind!=incidentFaceIndexs[tetr]) sortedIndexes[filled++]=ind;
                    }
                    sortIndexes(sortedIndexes,3,tets[tetr]); //sorts indexes in crescent order of vertexIndexes at that position
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"  SortedIndexes:..."<<endl;
//                    for(int i=0;i<3;i++){
//                        cerr<<"    index: "<<sortedIndexes[i]<<", vertex position: "<<tets[tetr]->TV(sortedIndexes[i])<<"..."<<endl;
//                    }
//#endif
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"  Common face is at index "<<incidentFaceIndexs[tetr]<<"..."<<endl;
//#endif
                    for(int i=0; i<3; i++){
                        PositionIndex newt=sortedIndexes[i];
                        int actualIndex=i;
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"    number"<<finalIndex(actualIndex,tetr)<<" with index "<<actualIndex<<" and tetra "<<tetr<<": "<<endl;
//#endif
                        Triangle face=tets[tetr]->TF(newt);
                        
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      creating..."<<endl;
//#endif
                        
                        PositionIndex newVertexPosition=newt;
                        newTetras[actualIndex][tetr]=new Tetra(vertexIndex,newVertexPosition,tets[tetr]);
                        
                        tetras->at(finalIndex(actualIndex,tetr))=newTetras[actualIndex][tetr];
                        
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      setting TT..."<<endl;
//#endif
                        newTetras[actualIndex][tetr]->setTT(incidentFaceIndexs[tetr],finalIndex(actualIndex,1-tetr));
                        newTetras[actualIndex][tetr]->setTT(sortedIndexes[(i+1)%3],finalIndex((actualIndex+1)%3,tetr));
                        newTetras[actualIndex][tetr]->setTT(sortedIndexes[(i+2)%3],finalIndex((actualIndex+2)%3,tetr));
                        
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      setting inverse TT..."<<endl;
//#endif
                        newTetras[actualIndex][tetr]->setInverseTT(incidentFaceIndexs[tetr],incidentFaceIndexs[1-tetr]);
                        newTetras[actualIndex][tetr]->setInverseTT(sortedIndexes[(i+1)%3],newt);
                        newTetras[actualIndex][tetr]->setInverseTT(sortedIndexes[(i+2)%3],newt);
                        
                        if(i==0){
//#ifdef LIBTETRA_DEBUG
//                            cerr<<"      fixing old vertex VTstar..."<<endl;
//#endif
                            getVertex(tets[tetr]->TV(newVertexPosition))->setVTstar(finalIndex(1,tetr));
                        }
                    }
                }
#undef finalIndex
                delete tets[0];
                delete tets[1];
            }else{
                for(int i=0; i<2; i++) pushTetra(NULL);
                
                Triangle face;
                Tetra* newTetras[3][1];
                
#define finalIndex(A,B) ( A==0 ? tetraIndexes[B] : tnum+A-1 )
                
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Creating new Tetras..."<<endl;
//#endif
                int tetr=0;
                int actualIndex=0;
                PositionIndex sortedIndexes[3];
                int filled=0;
                for(PositionIndex ind=a; ind<=d; ind=(PositionIndex)(ind+1)){
                    if(ind!=incidentFaceIndexs[tetr]) sortedIndexes[filled++]=ind;
                }
                sortIndexes(sortedIndexes,3,tets[tetr]); //sorts indexes in crescent order of vertexIndexes at that position

//#ifdef LIBTETRA_DEBUGc
//                err<<"  SortedIndexes:..."<<endl;
//#endif
//                for(int i=0;i<3;i++){
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    index: "<<sortedIndexes[i]<<", vertex position: "<<tets[tetr]->TV(sortedIndexes[i])<<"..."<<endl;
//#endif
//                }
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  Common face is at index "<<incidentFaceIndexs[tetr]<<"..."<<endl;
//#endif
                for(int i=0; i<3; i++){
                    PositionIndex newt=sortedIndexes[i];
                    int actualIndex=i;
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    number"<<finalIndex(actualIndex,tetr)<<" with index "<<actualIndex<<" and tetra "<<tetr<<": "<<endl;
//#endif
                    Triangle face=tets[tetr]->TF(newt);

//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      creating..."<<endl;
//#endif

                    PositionIndex newVertexPosition=newt;
                    newTetras[actualIndex][tetr]=new Tetra(vertexIndex,newVertexPosition,tets[tetr]);

                    tetras->at(finalIndex(actualIndex,tetr))=newTetras[actualIndex][tetr];

//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      setting TT..."<<endl;
//#endif
                    newTetras[actualIndex][tetr]->setTT(sortedIndexes[(i+1)%3],finalIndex((actualIndex+1)%3,tetr));
                    newTetras[actualIndex][tetr]->setTT(sortedIndexes[(i+2)%3],finalIndex((actualIndex+2)%3,tetr));

//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      setting inverse TT..."<<endl;
//#endif
                    newTetras[actualIndex][tetr]->setInverseTT(sortedIndexes[(i+1)%3],newt);
                    newTetras[actualIndex][tetr]->setInverseTT(sortedIndexes[(i+2)%3],newt);

//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      fixing old vertex VTstar..."<<endl;
//#endif
                    getVertex(tets[tetr]->TV(newVertexPosition))->setVTstar(finalIndex(1,tetr));
                }
#undef finalIndex
                delete tets[0];
            }
            delete (Triangle*)faceOrEdgeVertexIsOn;
            return noException;
        }
        if(relativePosition==edge){
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  Point is on egde "<<*(Edge*)faceOrEdgeVertexIsOn<<" of tetrahedron "<<tetraIndex<<endl;
//#endif
            VertexIndex vnum=vertexNumber();
            VertexIndex vertexIndex=vnum;
            TetraIndex tnum=tetraNumber();
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  retrieving incident tetra indexes..."<<endl;
//#endif
            
            Edge* commonEdge=((Edge*)faceOrEdgeVertexIsOn);
            Tetra* oldTetra;
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  getting ET"<<endl;
//#endif
            std::vector<TetraIndex> indexes=ET(*commonEdge,tetraIndex); // already sorted
            int adjnum=indexes.size();
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  there are "<<adjnum<<" tetras in the star of given edge"<<endl;
//#endif
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  storing old tetras"<<endl;
//#endif
            std::vector<Tetra*> oldTetras;
            for(int i=0; i<adjnum;i++){
                oldTetras.push_back(getTetra(indexes[i]));
            }
            
            std::vector<Tetra*>* newTetras[2];
            for(int j=0; j<2; j++){ // j==0 -> northern, j==1 -> southern
                newTetras[j]=new std::vector<Tetra*>();
            }
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  adding vertex and setting VTstar"<<endl;
//#endif
            vertex->setVTstar(tetraIndex);
            pushVertex(vertex);
            
            for(int i=0; i<adjnum; i++) pushTetra(NULL);
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  getting poles"<<endl;
//#endif
            VertexIndex poles[2]; 
            poles[0]=commonEdge->minindex(); //North
            poles[1]=commonEdge->maxindex(); //South
            delete commonEdge;
            PositionIndex northPoleIndex=a; // pole belonging to tetra
            PositionIndex southPoleIndex=b; // pole replaced by new vertex
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  creating tetras"<<endl;
//#endif
            // creating Tetras
#define finalIndex(A,B) ( B==0 ? indexes[A] : (tnum+A) )
            for(int j=0; j<2; j++){
                for(int i=0; i<adjnum; i++){
                    oldTetra=oldTetras[i];
                    
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    looking for 'south pole' of tetra "<<i<<" of chain "<<j<<endl;
//#endif
                    for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                        if(oldTetra->TV(pi)==poles[j]){
                            northPoleIndex=pi;
                            continue;
                        }
                        if(oldTetra->TV(pi)==poles[1-j]){
                            southPoleIndex=pi;
                            continue;
                        }
                    }
                    
                    // "north" relations are set during creation"
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    creating tetra"<<endl;
//#endif
                    Tetra* newTetra=new Tetra(vertexIndex,southPoleIndex,oldTetra);
                    newTetras[j]->push_back(newTetra);
                    tetras->at(finalIndex(i,j))=newTetra;
                }
            }

            // the north and south poles used are the ones of newTetras[1]
            getVertex(newTetras[1]->at(0)->TV(northPoleIndex))->setVTstar(finalIndex(0,1));
            
            // relationships setting variables
            PositionIndex eastPoleIndex; // adjacent to following tetra
            PositionIndex westPoleIndex; // adjacent to preceding tetra
            PositionIndex oldWestPoleIndex;
            PositionIndex oldEastPoleIndex;
            int filled=0;
            
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  presetting initial values for east, west and oldWest poles"<<endl;
//#endif
            //presetting westPoleIndex and eastPoleIndex
            oldTetra=oldTetras[0];
            for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                int j=0;
                if(oldTetra->TV(pi)==poles[1-j]){
                    continue;
                }
                if(oldTetra->TV(pi)==poles[j]){
                    continue;
                }
                if(filled++==0){
                    westPoleIndex=pi;
                    continue;
                }
                eastPoleIndex=pi;
            }
            
            //special case, the north and south poles used are the ones of newTetras[1]
            if(adjnum==1){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"  [SPECIAL CASE]: there is only one tetra in the star of given edge"<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    setting relations"<<endl;
//#endif
                newTetras[0]->at(0)->setTT(southPoleIndex,finalIndex(0,1));
                newTetras[0]->at(0)->setInverseTT(southPoleIndex,northPoleIndex);
                newTetras[1]->at(0)->setTT(northPoleIndex,finalIndex(0,0));
                newTetras[1]->at(0)->setInverseTT(northPoleIndex,southPoleIndex);
                
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    cleaning..."<<endl;
//#endif
                delete oldTetras[0];
                delete newTetras[0];
                delete newTetras[1];
                return noException;
            }
            
            bool continuous=true;
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  checking correct order for first tetra east and west poles"<<endl;
//#endif
            // correct order, just discontinuous
            if(oldTetra->TT(westPoleIndex)==noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    they were in the correct order"<<endl;
//#endif
                continuous=false;
            }
            // this means west and east poles were inverted
            if(oldTetra->TT(eastPoleIndex)==noAdjacentTetraIndex){
                PositionIndex tmp=eastPoleIndex;
                eastPoleIndex=westPoleIndex;
                westPoleIndex=tmp;
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    they were in the wrong order and have been switched"<<endl;
//#endif
                continuous=false;
            }
            oldEastPoleIndex=westPoleIndex;
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  checking if adjacency graph of edge* tetras is homeomorphic to a circle..."<<endl;
//#endif
            if(continuous){
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    it is homeomorphic to a circle"<<endl;
//#endif
                oldTetra=oldTetras[adjnum-1];
                for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
                    int j=0;
                    if(oldTetra->TV(pi)==poles[1-j]){
                        continue;
                    }
                    if(oldTetra->TV(pi)==poles[j]){
                        continue;
                    }
                    if(oldTetra->TV(pi)==oldTetras[0]->TV(eastPoleIndex)){
                        PositionIndex tmp=eastPoleIndex;
                        eastPoleIndex=westPoleIndex;
                        westPoleIndex=tmp;
                        pi=a;
                        continue;
                    }
                    if(oldTetra->TV(pi)==oldTetras[0]->TV(westPoleIndex)){
                        continue;
                    }
                    oldWestPoleIndex=pi;
                }
            }
            else{ 
//#ifdef LIBTETRA_DEBUG
//                cerr<<"    it is NOT homeomorphic to a circle"<<endl;
//#endif
            }
            
            PositionIndex initial_eastPoleIndex=eastPoleIndex; // adjacent to following tetra
            PositionIndex initial_westPoleIndex=westPoleIndex; // adjacent to preceding tetra
            PositionIndex initial_oldWestPoleIndex=oldWestPoleIndex;
            PositionIndex initial_oldEastPoleIndex=oldEastPoleIndex;
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  setting relations"<<endl;
//#endif
            // setting relations
            for(int j=0; j<2; j++){
                
                eastPoleIndex=initial_eastPoleIndex;
                westPoleIndex=initial_westPoleIndex; 
                oldWestPoleIndex=initial_oldWestPoleIndex;
                oldEastPoleIndex=initial_oldEastPoleIndex;
                
                for(int i=0; i<adjnum; i++){
                    oldTetra=oldTetras[i];
                    PositionIndex northPoleIndex=a; // pole belonging to tetra
                    PositionIndex southPoleIndex=b; // pole replaced by new vertex
                    int poleCovering[]={0,0,0,0};
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    Tetra "<<i<<" of chain "<<j<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"    getting pole indexes"<<endl;
//#endif
                    if(i>0){
                        oldWestPoleIndex=westPoleIndex;
                        oldEastPoleIndex=eastPoleIndex;
                    }
                    for(PositionIndex pi=a; pi<=d; pi=(PositionIndex)(pi+1)){
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      trying south..."<<endl;
//#endif
                        if(oldTetra->TV(pi)==poles[1-j]){
//#ifdef LIBTETRA_DEBUG
//                            cerr<<"      ...it's south"<<endl;
//#endif
                            poleCovering[1]++;
                            southPoleIndex=pi;
                            continue;
                        }
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      trying north..."<<endl;
//#endif
                        if(oldTetra->TV(pi)==poles[j]){
//#ifdef LIBTETRA_DEBUG
//                            cerr<<"      ...it's north"<<endl;
//#endif
                            poleCovering[0]++;
                            northPoleIndex=pi;
                            continue;
                        }
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      trying west..."<<endl;
//#endif
                        if(oldTetra->TT(pi)==noAdjacentTetraIndex){
                            if(i==adjnum-1){
//#ifdef LIBTETRA_DEBUG
//                                cerr<<"      ...it's west"<<endl;
//#endif
                                poleCovering[3]++;
                                westPoleIndex=pi;
                                continue;
                            }
                            if(i==0){
//#ifdef LIBTETRA_DEBUG
//                                cerr<<"      ...it's east"<<endl;
//#endif
                                poleCovering[2]++;
                                eastPoleIndex=pi;
                                continue;
                            }
                            return disconnectedEdgeStar;
                        }
                        if(oldTetra->TV(pi)==oldTetras[(i+adjnum-1)%adjnum]->TV(oldEastPoleIndex)){
//#ifdef LIBTETRA_DEBUG
//                            cerr<<"      ...finally it's west"<<endl;
//#endif
                            poleCovering[3]++;
                            westPoleIndex=pi;
                            continue;
                        }
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      trying east..."<<endl;
//#endif
                        if(!continuous && oldTetra->TT(pi)!=noAdjacentTetraIndex){
                            if(i==0){
//#ifdef LIBTETRA_DEBUG
//                                cerr<<"      ...it's west"<<endl;
//#endif
                                poleCovering[3]++;
                                westPoleIndex=pi;
                                continue;
                            }
                            if(i==adjnum-1){
//#ifdef LIBTETRA_DEBUG
//                                cerr<<"      ...it's east"<<endl;
//#endif
                                poleCovering[2]++;
                                eastPoleIndex=pi;
                                continue;
                            }
                        }
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      ...finally it's east"<<endl;
//#endif
                        poleCovering[2]++;
                        eastPoleIndex=pi;
                    }
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      pole covering: N"<<poleCovering[0]<<" S"<<poleCovering[1]<<
//                            " E"<<poleCovering[2]<<" W"<<poleCovering[3]<<endl;
//#endif
                    Tetra* newTetra=newTetras[j]->at(i);
//#ifdef LIBTETRA_DEBUG
//                    if(newTetra==NULL) cerr<<"    AAHH! a NULL pointer!!!!"<<endl;
//#endif
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      setting south relations"<<endl;
//#endif
                    // "south" relations setting...
                    newTetra->setTT(northPoleIndex,finalIndex(i,1-j));
                    newTetra->setInverseTT(northPoleIndex,southPoleIndex);
                    
                    // "east" and "west" relation setting
                    if(oldTetra->TT(westPoleIndex)!=noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      setting west TT"<<endl;
//#endif
                        newTetra->setTT(westPoleIndex,finalIndex((i+1)%adjnum,j));
                    }
                    if(oldTetra->TT(eastPoleIndex)!=noAdjacentTetraIndex){
//#ifdef LIBTETRA_DEBUG
//                        cerr<<"      setting east TT, inverseTT and previous tetra west InverseTT"<<endl;
//#endif
                        newTetra->setTT(eastPoleIndex,finalIndex((i+adjnum-1)%adjnum,j));
                        newTetra->setInverseTT(eastPoleIndex,oldWestPoleIndex);
                        newTetras[j]->at((i+adjnum-1)%adjnum)->setInverseTT(oldWestPoleIndex,eastPoleIndex);
                    }
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      fixing TT and inverseTT of adjacent tetra opposite to old south pole"<<endl;
//#endif
                    Tetra* adjTetra=getTetra(oldTetra->TT(southPoleIndex));
                    if(adjTetra!=NULL){
                        adjTetra->setTT(oldTetra->inverseTT(southPoleIndex),finalIndex(i,j));
                        adjTetra->setInverseTT(oldTetra->inverseTT(southPoleIndex),southPoleIndex);
                    }
//#ifdef LIBTETRA_DEBUG
//                    cerr<<"      fixing VTstar of old south pole"<<endl;
//#endif
                    getVertex(oldTetra->TV(northPoleIndex))->setVTstar(finalIndex(i,j));
                }
            }
#undef finalIndex
            
//#ifdef LIBTETRA_DEBUG
//            cerr<<"  cleaning"<<endl;
//#endif
            for(int i=0; i<adjnum; i++)delete oldTetras[i];
            delete newTetras[0];
            delete newTetras[1];
            return noException;
        }
        if(relativePosition==vertexPoint){
            return vertexAlreadyPresent;
        }
        if(relativePosition==external){
            return vertexExternalFromMesh;
        }
        return noException;
    }
    
    template<class Coordinate, class FieldValue>
    TVposition Mesh<Coordinate,FieldValue> :: getRelativePosition(Vertex<Coordinate,FieldValue>* vertex, Tetra* tetra, void** storeForIncidentFaceOrEdge){
//#ifdef LIBTETRA_DEBUG
//        cerr<<"Entering getRelativePosition..."<<endl;
//#endif
        Turn dTurn=collinear; // used as uninitialized value
        bool coll=false;
        Triangle collinearFace;
        for(PositionIndex i=a; i<=d; i=(PositionIndex)(i+1)){
            Triangle face=tetra->TF(i);
            Turn turn=getTurn(vertex,face);
            if(turn==collinear){
                coll=true;
                collinearFace=face;
                break;
            }
//#ifdef LIBTETRA_DEBUG
//            cerr<<"...external."<<endl;
//#endif
            if(turn!=dTurn && dTurn!=collinear) return external;
            dTurn=turn;
        }
        if(coll){
            dTurn=collinear;
            PositionIndex edgeIndex=unknown;
            Turn turn[3];
            int sum=0;
            for(PositionIndex i=a; i<d; i=(PositionIndex)(i+1)){
                turn[i]=getTurn(vertex,collinearFace,i);
                if(turn[i]==collinear){
                    edgeIndex=i;
                }
                sum=sum+turn[i];
            }

            if(sum==3 || sum==-3){
                *storeForIncidentFaceOrEdge=collinearFace.clone();
//#ifdef LIBTETRA_DEBUG
//                cerr<<"...face."<<endl;
//#endif
                return face; //(can't be all three collinear!)
            }
            if(sum==2 || sum==-2){
                *storeForIncidentFaceOrEdge=collinearFace.FE(edgeIndex);
//#ifdef LIBTETRA_DEBUG
//                cerr<<"...edge."<<endl;
//#endif
                return edge;
            }

//#ifdef LIBTETRA_DEBUG
//            cerr<<"...external(2)."<<endl;
//#endif
            return external;
        }
//#ifdef LIBTETRA_DEBUG
//        cerr<<"...internal."<<endl;
//#endif
        return internal;
    }

    template<class Coordinate, class FieldValue>
    std::vector<Triangle> Mesh<Coordinate,FieldValue> :: FF(Triangle face)
    {
        vector<Triangle> faces;
        vector<Triangle> tetra0 = this->EF(*face.FE(a));
        vector<Triangle> tetra1 = this->EF(*face.FE(b));
        vector<Triangle> tetra2 = this->EF(*face.FE(c));
        for(unsigned int i=0;i<tetra0.size();i++)
        {
            if(face != tetra0.at(i))
            {
                faces.push_back(tetra0.at(i));
            }
        }
        for(unsigned int i=0;i<tetra1.size();i++)
        {
            if(face != tetra1.at(i))
            {
                faces.push_back(tetra1.at(i));
            }
        }
        for(unsigned int i=0;i<tetra2.size();i++)
        {
            if(face != tetra2.at(i))
            {
                faces.push_back(tetra2.at(i));
            }
        }
        return faces;
    }

    template<class Coordinate, class FieldValue>
    std::vector<TetraIndex> Mesh<Coordinate,FieldValue> :: FT(Triangle face)
    {
        vector<TetraIndex> tetra;
        VertexIndex v0 = face.FV(a);
        bool* isBorder = new bool();
        vector<TetraIndex> all_tetra = this->VT(v0,isBorder);
        for(unsigned int i=0;i<all_tetra.size();i++)
        {
            for(PositionIndex j=a;j<=d;j=(PositionIndex)(j+1))
            {
                if(face == this->getTetra(all_tetra.at(i))->TF(j))
                {
                    tetra.push_back(all_tetra.at(i));
                    if(this->getTetra(all_tetra.at(i))->TT(j)!=noAdjacentTetraIndex)
                        tetra.push_back(this->getTetra(all_tetra.at(i))->TT(j));
                    return tetra;
                }
            }
        }
        //non dovrebbe mai arrivarci salvo se si inserisce una faccia che non appartiene a nessun tetra
        return tetra;
    }

    template<class Coordinate, class FieldValue>
    std::vector<Edge> Mesh<Coordinate,FieldValue> :: EE(Edge edge)
    {
        vector<Edge> edges;
        vector<Edge> edge0 = this->VE(edge.EV(a));
        vector<Edge> edge1 = this->VE(edge.EV(b));
        for(unsigned int i=0;i<edge0.size();i++)
        {
            if(edge != edge0.at(i))
                edges.push_back(edge0.at(i));
        }
        for(unsigned int i=0;i<edge1.size();i++)
        {
            if(edge != edge1.at(i))
                edges.push_back(edge1.at(i));
        }
        return edges;
    }

    template<class Coordinate, class FieldValue>
    std::vector<Triangle> Mesh<Coordinate,FieldValue> :: EF(Edge edge)
    {
        vector<Triangle> faces;
        vector<Triangle> all_faces = this->VF(edge.EV(a));
        for(unsigned int i=0;i<all_faces.size();i++)
        {
            for(PositionIndex pos = a; pos < d; pos = (PositionIndex)(pos+1))
            {
                if(edge == *(all_faces.at(i).FE(pos)))
                {
                    faces.push_back(all_faces.at(i));
                    break;
                }
            }
        }

        return faces;
    }

    template<class Coordinate, class FieldValue>
    std::vector<TetraIndex> Mesh<Coordinate,FieldValue> :: ET(Edge edge)
    {
        vector<TetraIndex> tetra;
        VertexIndex v0 = edge.EV(a);
        bool* isBorder = new bool();
        vector<TetraIndex> vt = this->VT(v0,isBorder);
        for(unsigned int w=0; w<vt.size(); w++)
        {
            bool isFound = false;
            for(PositionIndex i=a;i<=d;i=(PositionIndex)(i+1))
            {
                for(PositionIndex j=PositionIndex(i+b);j<=d;j=(PositionIndex)(j+1))
                {
                    Edge e = Edge(this->getTetra(vt.at(w))->TV(i),this->getTetra(vt.at(w))->TV(j));
                    if(edge == e)
                    {
                        tetra.push_back(vt.at(w));
                        isFound = true;
                        break;
                    }
                }
                if(isFound)
                    break;
            }
        }
        return tetra;
    }

    template<class Coordinate, class FieldValue>
    std::vector<VertexIndex> Mesh<Coordinate,FieldValue> :: VV(VertexIndex center)
    {
        vector<VertexIndex> vertices;

        vector<TetraIndex> tetra;//contiene i tetra visitati
        TetraIndex tet = this->getVertex(center)->VTstar();
        queue<TetraIndex> queue;
        queue.push(tet);
        tetra.push_back(tet);
        this->getTetra(tet)->setStatus(checked);

        PositionIndex k=(PositionIndex)-1;

        while(!queue.empty())
        {
            k=(PositionIndex)-1;
            TetraIndex current = queue.front();

            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(PositionIndex pos = a; pos <= d; pos = (PositionIndex)(pos+1))
            {
                if(this->getTetra(current)->TV(pos) == center)
                {
                    k = pos;
                    break;
                }
            }

            for(PositionIndex pos = b; pos <= d; pos = (PositionIndex)(pos+1))
            {
                VertexIndex v = this->getTetra(current)->TV(((PositionIndex)(k+pos)%4));
                if(!this->isIntoVector(v,vertices))
                    vertices.push_back(v);

                TetraIndex adj = this->getTetra(current)->TT((PositionIndex)((k+pos)%4));
                if( adj != noAdjacentTetraIndex && this->getTetra(adj)->status() != checked)
                {
                    this->getTetra(adj)->setStatus(checked);
                    queue.push(adj);
                    tetra.push_back(adj);
                }
            }
            queue.pop();
        }

        for(unsigned int i=0; i<tetra.size(); i++)
            this->getTetra(tetra.at(i))->setStatus(unchecked);

        return vertices;
    }

    template<class Coordinate, class FieldValue>
    std::vector<Edge> Mesh<Coordinate,FieldValue> :: VE(VertexIndex center)
    {
        vector<Edge> edges;
        vector<TetraIndex> tetra;//contiene i tetra visitati
        TetraIndex tet = this->getVertex(center)->VTstar();
        queue<TetraIndex> queue;
        queue.push(tet);
        tetra.push_back(tet);
        this->getTetra(tet)->setStatus(checked);

        PositionIndex k=(PositionIndex)-1;

        while(!queue.empty())
        {
            k=(PositionIndex)-1;
            TetraIndex current = queue.front();

            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(PositionIndex pos = a; pos <= d; pos = (PositionIndex)(pos+1))
            {
                if(this->getTetra(current)->TV(pos) == center)
                {
                    k = pos;
                    break;
                }
            }

            for(PositionIndex pos = b; pos <= d; pos = (PositionIndex)(pos+1))
            {
                Edge e = Edge(this->getTetra(current)->TV(k),this->getTetra(current)->TV(((PositionIndex)(k+pos)%4)));
                if(!this->isIntoVector(e,edges))
                    edges.push_back(e);

                TetraIndex adj = this->getTetra(current)->TT((PositionIndex)((k+pos)%4));
                if( adj != noAdjacentTetraIndex && this->getTetra(adj)->status() != checked)
                {
                    this->getTetra(adj)->setStatus(checked);
                    queue.push(adj);
                    tetra.push_back(adj);
                }
            }
            queue.pop();
        }

        for(unsigned int i=0; i<tetra.size(); i++)
            this->getTetra(tetra.at(i))->setStatus(unchecked);

        return edges;
    }

    template<class Coordinate, class FieldValue>
    std::vector<Triangle> Mesh<Coordinate,FieldValue> :: VF(VertexIndex center)
    {
        vector<Triangle> faces;
        vector<TetraIndex> tetra;//contiene i tetra visitati
        TetraIndex tet = this->getVertex(center)->VTstar();
        queue<TetraIndex> queue;
        queue.push(tet);
        tetra.push_back(tet);
        this->getTetra(tet)->setStatus(checked);

        PositionIndex k=(PositionIndex)-1;

        while(!queue.empty())
        {
            k=(PositionIndex)-1;

            TetraIndex current = queue.front();

            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(PositionIndex pos = a; pos <= d; pos = (PositionIndex)(pos+1))
            {
                if(this->getTetra(current)->TV(pos) == center)
                {
                    k = pos;
                    break;
                }
            }

            for(PositionIndex pos = b; pos <= d; pos = (PositionIndex)(pos+1))
            {
                Triangle faceIn = this->getTetra(current)->TF((PositionIndex)((k+pos)%4));
                if(!this->isIntoVector(faceIn,faces))
                    faces.push_back(faceIn);

                TetraIndex adj = this->getTetra(current)->TT((PositionIndex)((k+pos)%4));
                if( adj != noAdjacentTetraIndex && this->getTetra(adj)->status() != checked)
                {
                    this->getTetra(adj)->setStatus(checked);
                    queue.push(adj);
                    tetra.push_back(adj);
                }
            }
            queue.pop();
        }

        for(unsigned int i=0; i<tetra.size(); i++)
            this->getTetra(tetra.at(i))->setStatus(unchecked);

        return faces;
    }

    template<class Coordinate, class FieldValue>
    std::vector<TetraIndex> Mesh<Coordinate,FieldValue> :: VT(VertexIndex center,bool* foundBorder)
    {
        *foundBorder = false;
        vector<TetraIndex> tetras;
        TetraIndex tet = this->getVertex(center)->VTstar();
        queue<TetraIndex> coda;
        coda.push(tet);
        tetras.push_back(tet);
        this->getTetra(tet)->setStatus(checked);

        PositionIndex k=(PositionIndex)-1;

        while(!coda.empty()){            
            k=(PositionIndex)-1;
            TetraIndex attuale = coda.front();

            //cerco la posizione del vertice nell'array dei vertici del triangolo
            for(PositionIndex pos = a; pos <= d; pos = (PositionIndex)(pos+1))
            {
                if(this->getTetra(attuale)->TV(pos) == center)
                {
                    k = pos;
                    break;
                }
            }

            for(PositionIndex pos = b; pos <= d; pos = (PositionIndex)(pos+1))
            {
                TetraIndex adj = this->getTetra(attuale)->TT((PositionIndex)((k+pos)%4));
                if( adj != noAdjacentTetraIndex && this->getTetra(adj)->status() != checked)
                {
                    this->getTetra(adj)->setStatus(checked);
                    coda.push(adj);
                    tetras.push_back(adj);
                }
                else if(adj == noAdjacentTetraIndex && this->getTetra(tet)->TV((PositionIndex)((k+pos)%4)) != center)
                    *foundBorder = true;
            }
            coda.pop();
        }

        for(unsigned int i=0; i<tetras.size(); i++)
            this->getTetra(tetras.at(i))->setStatus(unchecked);

        return tetras;
    }

    //UTILITY FUNCTION
    template<class Coordinate, class FieldValue>
    template<class C> bool Mesh<Coordinate,FieldValue> :: isIntoVector(C c, vector<C> &c_vect)
    {
        for(unsigned int i=0; i<c_vect.size();i++)
            if(c == c_vect.at(i))
                return true;
        return false;
    }

}

#endif //LIBTETRA_Mesh
