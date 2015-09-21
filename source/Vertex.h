#ifndef VERTEX_H
#define VERTEX_H
#include "Matrix.cpp"
#include "Define.h"

namespace LibTetra{
    template<class Coordinate, class FieldValue>
    // Class storing geometrical information of vertices
    class Vertex{
    private:
        Coordinate vertexCoordinates[3]; // Coordinates of vertices
        FieldValue vertexFieldValue; // Scalar Field value over vertex
        TetraIndex vertexVTstar; // VT* relationship

        float vertexDistortion;
//        bool isvertexBoundary;

    public:

        // Constructor
        inline Vertex(Coordinate x, Coordinate y, Coordinate z,FieldValue value){
            vertexCoordinates[0]=x;
            vertexCoordinates[1]=y;
            vertexCoordinates[2]=z;
            setFieldValue(value);
            setVTstar(uninitializedTetraIndex);
            ///------------///
            //ADDED IN 20 May 2011
            vertexDistortion = -1;
//            isvertexBoundary = false;
        }

        // Field value setter
        inline void setFieldValue(FieldValue value){
            vertexFieldValue=value;
        }

        // VT* setter
        inline void setVTstar(TetraIndex VTstar){
            vertexVTstar=VTstar;
        }

        // First coordinate getter
        inline Coordinate x(){
            return vertexCoordinates[0];
        }

        // Second coordinate getter
        inline Coordinate y(){
            return vertexCoordinates[1];
        }

        // Third coordinate getter
        inline Coordinate z(){
            return vertexCoordinates[2];
        }

        // Field value getter
        inline FieldValue fieldValue(){
            return vertexFieldValue;
        }

        // VT* getter
        inline TetraIndex VTstar(){
            return vertexVTstar;
        }

        inline Coordinate EuclideanDistance3D(Vertex other){

            return sqrt((this->x()-other.x())*(this->x()-other.x()) + (this->y()-other.y())*(this->y()-other.y()) + (this->z()-other.z())*(this->z()-other.z()));
        }

        ///------------///
        //ADDED IN 20 May 2011
        inline void setDistortion(float dist)
        {
            this->vertexDistortion = dist;
        }

//        inline void setIsBoundary(bool isB)
//        {
//            this->isvertexBoundary = isB;
//        }

        inline float distortion()
        {
            return this->vertexDistortion;
        }

//        inline bool isBoundary()
//        {
//            return this->isvertexBoundary;
//        }
        ///------------///

        // Coordinate equality
        inline bool operator==(Vertex vertex){
            return x()==vertex.x() && y()==vertex.y() && z()==vertex.z();
        }

        // Standard output operator
        inline friend ostream & operator << (ostream & out, Vertex & vertex){
            out << "Vertex[(" << vertex.x() << " " << vertex.y() << " " << vertex.z() << " " << vertex.fieldValue() << ")"<<endl;
            out << "VT*:" << vertex.VTstar() <<endl;
            out << "]";
            return out;
        }

        // Difference operator
        inline Vertex* operator-(Vertex other){
            return new Vertex(x()-other.x(),y()-other.y(),z()-other.z(),fieldValue()-other.fieldValue());
        }

    };
}

#endif // VERTEX_H
