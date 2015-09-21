#ifndef TETRA_H
#define TETRA_H

//#include "Define.h"
#include "Triangle.h"

#define coutTInd(C,T) if(T==noAdjacentTetraIndex)C<<"None";else C<<T;
#define coutPosInd(C,P) if(P==unknown)C<<"Unknown";else C<<P;

namespace LibTetra{
    // Class storing complete topological information for tetrahedra (top simplices) with reference to geometrical information of their vertices
    class Tetra{
    private:
        VertexIndex tetraTV[4]; // Indices of geometrical information of the vertices
        TetraIndex tetraTT[4]; // Adjacency relations
        PositionIndex tetraInverseTT[4]; // PositionIndex of face shared with each adjacent tetrahedron, with respect to the other tetrahedron
        TetraStatus tetraStatus; // Checked status of tetrahedron
        //TODO : this must be optimized

    public:

        // Constructor
        inline Tetra(VertexIndex a, VertexIndex b, VertexIndex c, VertexIndex d){
            tetraTV[0]=a;
            tetraTV[1]=b;
            tetraTV[2]=c;
            tetraTV[3]=d;
            setTT(noAdjacentTetraIndex,noAdjacentTetraIndex,noAdjacentTetraIndex,noAdjacentTetraIndex);
            setInverseTT(unknown,unknown,unknown,unknown);
            setStatus(unchecked);
        }

        // Vertex replacement constructor
        inline Tetra(VertexIndex Avertex, PositionIndex vertexAFinalPosition, Tetra* oldTetra){
            for(PositionIndex ind=a; ind<=d;ind=(PositionIndex)(ind+1)){
                tetraTV[ind]=oldTetra->TV(ind);
            }
            tetraTV[vertexAFinalPosition]=Avertex;
            setTT(noAdjacentTetraIndex,noAdjacentTetraIndex,noAdjacentTetraIndex,noAdjacentTetraIndex);
            setTT(vertexAFinalPosition,oldTetra->TT(vertexAFinalPosition));
            setInverseTT(unknown,unknown,unknown,unknown);
            setInverseTT(vertexAFinalPosition,oldTetra->inverseTT(vertexAFinalPosition));
            setStatus(unchecked);
        }

        // Single vertex index setter
        inline void setTV(PositionIndex vertexPosition, VertexIndex vertex){
            tetraTV[vertexPosition]=vertex;
        }

        // Complete adjacency relations setter
        inline void setTT(TetraIndex a, TetraIndex b, TetraIndex c, TetraIndex d){
            tetraTT[0]=a;
            tetraTT[1]=b;
            tetraTT[2]=c;
            tetraTT[3]=d;
        }

        // Single adjacency relation setter
        inline void setTT(PositionIndex oppositeVertexPosition, TetraIndex tetra){
            tetraTT[oppositeVertexPosition]=tetra;
        }

        // Complete inverse adjacency relations setter
        inline void setInverseTT(PositionIndex a, PositionIndex b, PositionIndex c, PositionIndex d){
            tetraInverseTT[0]=a;
            tetraInverseTT[1]=b;
            tetraInverseTT[2]=c;
            tetraInverseTT[3]=d;
        }

        // Single inverse adjacency relation setter
        inline void setInverseTT(PositionIndex oppositeVertexPosition, PositionIndex position){
            tetraInverseTT[oppositeVertexPosition]=position;
        }

        // Check status setter
        inline void setStatus(TetraStatus status){
            tetraStatus=status;
        }

        // Vertices getter
        inline VertexIndex TV(PositionIndex vertexPosition){
            return tetraTV[vertexPosition];
        }

        // Adjacent tetrahedra getter
        inline TetraIndex TT(PositionIndex TT){
            return tetraTT[TT];
        }

        // Inverse adjacency relations getter
        inline PositionIndex inverseTT(PositionIndex TT){
            return tetraInverseTT[TT];
        }

        // Belonging of index to tetrahedron
        inline bool contains(VertexIndex vertex){
            for(PositionIndex pindex=a; pindex<=d; pindex=(PositionIndex)(pindex+1)){
                if(TV(pindex)==vertex) return true;
            }
            return false;
        }

        // Belonging of a triple of vertex indices to tetrahedron
        inline bool contains(Triangle face){
//            int vertexesContained=0;
            for(PositionIndex pindex=a; pindex<=d; pindex=(PositionIndex)(pindex+1)){
                if(face.contains(TV(pindex))) return true;
            }
            return false;
        }

        // Face getter
        inline Triangle TF(PositionIndex faceIndex){
            PositionIndex indexes[3];
            int filled=0;
            for(PositionIndex index=a; index<=d; index=(PositionIndex)(index+1)){
                if(index!=faceIndex) indexes[filled++]=index;
            }
            Triangle t= Triangle(TV(indexes[0]),TV(indexes[1]),TV(indexes[2]),faceIndex);
            return t;
        }

         // Check status getter
        inline TetraStatus status(){
            return tetraStatus;
        }

        // Standard output operator getter
        inline friend ostream & operator << (ostream & out, Tetra & tetra){
            out
                    << "Tetra[" <<endl
                    << " TV("
                    << tetra.TV(a) << " " << tetra.TV(b) << " " << tetra.TV(c) << " " << tetra.TV(d) << ") " <<endl
                    << " TT(";
            coutTInd(cout,tetra.TT(a))
                    cout<< " ";
            coutTInd(cout,tetra.TT(b))
                    cout << " ";
            coutTInd(cout,tetra.TT(c))
                    cout << " ";
            coutTInd(cout,tetra.TT(d))
                    cout << ") " <<endl;
            cout << " InverseTT(";
            coutPosInd(cout,tetra.inverseTT(a))
                    cout << " ";
            coutPosInd(cout,tetra.inverseTT(b))
                    cout << " ";
            coutPosInd(cout,tetra.inverseTT(c))
                    cout << " ";
            coutPosInd(cout,tetra.inverseTT(d))
                    cout << ") " <<endl;
            cout <<"]";
            return out;
        }
    };
}

#endif // TETRA_H
