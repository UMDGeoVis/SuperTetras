#ifndef EDGE_H
#define EDGE_H
#include "Vertex.h"

namespace LibTetra{
    // Class for temporary storage of vertices 1-adjacency topological relation
    class Edge{
    private:
        VertexIndex vertexes[2]; // Indices of the endpoints
    public:

        // Constructor
        inline Edge(VertexIndex vi1, VertexIndex vi2){
            vertexes[0]=min(vi1,vi2);
            vertexes[1]=max(vi1,vi2);
        }

        // Smaller-index endpoint getter
        inline VertexIndex minindex(){return vertexes[0];}

        // Bigger-index endpoint getter
        inline VertexIndex maxindex(){return vertexes[1];}

        inline VertexIndex EV(PositionIndex vPos)
        {
            return vertexes[vPos];
        }

        // Index equality
        inline bool operator==(Edge edge){
            if(minindex()!=edge.minindex()) return false;
            if(maxindex()!=edge.maxindex()) return false;
            return true;
        }

        // Index inequality
        inline bool operator!=(Edge edge){
            return !((*this)==(edge));
        }

        // Standard output operator
        inline friend ostream & operator << (ostream & out, Edge & edge){
            out << "Edge[(" << edge.minindex() << " " << edge.maxindex() << ")]"<<endl;
            return out;
        }

    };
}

#endif // EDGE_H
