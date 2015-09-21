#ifndef DEFINE_H
#define DEFINE_H

namespace LibTetra{
    typedef unsigned long VertexIndex;
    typedef unsigned long TetraIndex;
    static TetraIndex reservedTetraIndexes=2;
    static TetraIndex maxTetraNumber=1000000000-reservedTetraIndexes; // TODO: change 1000000000 in MAX_LONG_VALUE (or whatever name it has....)
//    static long maxTetraNumber=-1;
    static TetraIndex noAdjacentTetraIndex=maxTetraNumber+1;
    static TetraIndex uninitializedTetraIndex=maxTetraNumber+2;

    // Position of a face or a vertex in a tetrahedron
    typedef enum{a,b,c,d,unknown} PositionIndex;

    // Check status of a tetrahedron (for adjacency graph spanning operations)
    typedef enum{unchecked, checked} TetraStatus;

    // Turn value
    typedef enum{collinear,left,right} Turn;

    // Incoherences
    typedef enum{missingVTstars} Incoherence;

    // Exception codes
    typedef enum{
          noException,
          someVTstarCouldNotBeSet,
          someTTorInverseTTCouldNotBeSet,
          vertexAlreadyPresent,
          inverseTTrelationBroken,
          vertexExternalFromMesh,
          disconnectedEdgeStar
          } Exception;

    // Relative position of a vertex with respect to a tetrahedron
    typedef enum{vertexPoint,edge,face,internal,external} TVposition;
}

#endif // DEFINE_H
