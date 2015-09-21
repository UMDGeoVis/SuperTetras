#ifndef TRIANGLE_H
#define TRIANGLE_H

//#include "Define.h"
#include "Edge.h"

namespace LibTetra{
    // Class for temporary storage of vertices 2-adjacency
    class Triangle{
    private:
        VertexIndex vertexes[3]; // Vertex indices
        PositionIndex fIndex; // Position index of the face in the original tetrahedron (if appliable)

    public:

        // Constructor
        inline Triangle(VertexIndex vi1, VertexIndex vi2, VertexIndex vi3, PositionIndex index=unknown){
            vertexes[0]=min(min(vi1,vi2),vi3);
            vertexes[2]=max(max(vi1,vi2),vi3);
            if(vi1>minindex()&&vi1<maxindex()) vertexes[1]=vi1;
            else if(vi2>minindex()&&vi2<maxindex()) vertexes[1]=vi2;
            else vertexes[1]=vi3;
            fIndex=index;
        }

        // Null constructor
        inline Triangle(){
            vertexes[0]=vertexes[1]=vertexes[2]=0;
            fIndex=unknown;
        }

        // Swallow copy
        inline Triangle* clone(){
            return new Triangle(minindex(),middleindex(),maxindex(),faceIndex());
        }

        // Smaller-index vertex getter
        inline VertexIndex minindex(){return vertexes[0];}

        // Middle-index vertex getter
        inline VertexIndex middleindex(){return vertexes[1];}

        // Bigger-index vertex getter
        inline VertexIndex maxindex(){return vertexes[2];}

        // PositionIndex getter
        inline PositionIndex faceIndex(){return fIndex;}

        // Vertices getter
        inline VertexIndex FV(PositionIndex index){return vertexes[index%3];}

        // Bounding edges getter
        inline Edge* FE(PositionIndex index){return new Edge(vertexes[(index+1)%3],vertexes[(index+2)%3]);}

        // Index equality
        inline bool operator==(Triangle triangle){
            if(minindex()!=triangle.minindex()) return false;
            if(middleindex()!=triangle.middleindex()) return false;
            if(maxindex()!=triangle.maxindex()) return false;
            return true;
        }

        // Index inequality
        inline bool operator!=(Triangle triangle){
            return !((*this)==(triangle));
        }

        // Belonging of an index to the face
        inline bool contains(VertexIndex vertex){
            for(PositionIndex pindex=a; pindex<d; pindex=(PositionIndex)(pindex+1)){
                if(FV(pindex)==vertex) return true;
            }
            return false;
        }

        // Standard output operator
        inline friend ostream & operator << (ostream & out, Triangle & triangle){
            out << "Triangle[(" << triangle.minindex() << " " << triangle.middleindex() << " " << triangle.maxindex() << ")]";
            return out;
        }

    };
}

#endif // TRIANGLE_H
