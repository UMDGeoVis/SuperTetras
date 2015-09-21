#include "Mesh.h"

namespace LibTetra{

template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeTS(char* file_name, vector<Coordinate> fieldvalues){


    FILE* file;
    file = fopen(file_name, "w");

    fprintf(file,"%d %d\n", vertexNumber(), tetraNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f %f %f %f\n", getVertex(i)->x(), getVertex(i)->y(), getVertex(i)->z(), fieldvalues[i]);


    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%ld %ld %ld %ld\n", getTetra(i)->TV(a), getTetra(i)->TV(b), getTetra(i)->TV(c), getTetra(i)->TV(d));

    fclose(file);
}


template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeField(char* file_name, vector<Coordinate> fieldvalues){

    FILE* file;
    file = fopen(file_name, "w");

    fprintf(file, "%d\n", fieldvalues.size());

    for(int i=0; i<fieldvalues.size(); i++)
        fprintf(file, "%f\n", fieldvalues[i]);

    fclose(file);
}


template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeVTK(char* file_name, vector<Coordinate> fieldvalues){

    FILE* file;
    file = fopen(file_name, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i)->x(), getVertex(i)->y(), getVertex(i)->z());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", tetraNumber(), tetraNumber()*5);

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "4 %d %d %d %d\n", getTetra(i)->TV(a), getTetra(i)->TV(b), getTetra(i)->TV(c), getTetra(i)->TV(d));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%d ", 10);
    fprintf(file, "\n\n");


    fprintf(file, "POINT_DATA %d \n", vertexNumber());
    fprintf(file, "FIELD FieldData 2\n");
    fprintf(file, "original_field 1 %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f ", getVertex(i)->fieldValue());

    fprintf(file, "\n\n");
    fprintf(file, "distortion_field 1 %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%lf ", fieldvalues[i]);

    fprintf(file, "\n");

    fclose(file);
}

//template<class Coordinate, class FieldValue>
//void Mesh<Coordinate,FieldValue>::output_segmentation(const char* filename, int sz){

//    FILE* file;
//    file = fopen(filename, "w");

//    fprintf(file, "%d \n", this->tetraNumber());
//    fprintf(file, "%d \n", sz);
//    for(int i=0; i<this->tetraNumber(); i++)
//        fprintf(file, "%d \n", this->clusterIndex[i]);

//    fprintf(file, "\n");

//    //adjacencies
//    vector< vector<int> > V=this->buildAdjacencies();

//    for(unsigned int ii; ii<V.size(); ii++){
//        vector<int> row = V.at(ii);
//        for(unsigned int jj=0;jj<row.size();jj++)
//            fprintf(file, "%d ", row.at(jj));
//        fprintf(file, "\n");
//    }

//    fclose(file);
//}

template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeSuperTetrasVTK(char *filename){

    FILE* file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i)->x(), getVertex(i)->y(), getVertex(i)->z());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", tetraNumber(), tetraNumber()*5);

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "4 %d %d %d %d\n", getTetra(i)->TV(a), getTetra(i)->TV(b), getTetra(i)->TV(c), getTetra(i)->TV(d));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%d\n ", 10);
    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", vertexNumber());
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "fieldvalue 1 %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f ", getVertex(i)->fieldValue());

    fprintf(file, "\n\n");

    fprintf(file, "CELL_DATA %d \n", tetraNumber());
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "segmentation 1 %d float\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%f ", this->clusterFieldValue[i]);
//    fprintf(file, "segmentation 1 %d int\n", tetraNumber());

//    for(int i=0; i<tetraNumber(); i++)
//        fprintf(file, "%d ", this->clusterIndex[i]);

    fprintf(file, "\n\n");

    fclose(file);
}

template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeIndicesVTK(char *filename){

    FILE* file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII\n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i)->x(), getVertex(i)->y(), getVertex(i)->z());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", tetraNumber(), tetraNumber()*5);

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "4 %d %d %d %d\n", getTetra(i)->TV(a), getTetra(i)->TV(b), getTetra(i)->TV(c), getTetra(i)->TV(d));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%d\n ", 10);
    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", vertexNumber());
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "fieldvalue 1 %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f ", getVertex(i)->fieldValue());

    fprintf(file, "\n\n");

    fprintf(file, "CELL_DATA %d \n", tetraNumber());
    fprintf(file, "FIELD FieldData 1\n");

//    fprintf(file, "segmentation 1 %d float\n", tetraNumber());

//    for(int i=0; i<tetraNumber(); i++)
//        fprintf(file, "%f ", this->clusterFieldValue[i]);
    fprintf(file, "segmentation 1 %d int\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%d ", this->clusterIndex[i]);

    fprintf(file, "\n\n");

    fclose(file);
}

template<class Coordinate, class FieldValue>
void Mesh<Coordinate,FieldValue>::writeVTK_segmentation(char* filename, vector<int>* ascending, vector<int>* descending){
    FILE* file;
    file = fopen(filename, "w");

    fprintf(file, "# vtk DataFile Version 2.0\n\n");
    fprintf(file, "ASCII \n");
    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
    fprintf(file, "POINTS %d float\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%f %f %f\n", getVertex(i)->x(), getVertex(i)->y(), getVertex(i)->z());
    fprintf(file, "\n\n");

    fprintf(file, "CELLS %d %d\n", tetraNumber(), tetraNumber()*5);

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "4 %d %d %d %d\n", getTetra(i)->TV(a), getTetra(i)->TV(b), getTetra(i)->TV(c), getTetra(i)->TV(d));
    fprintf(file, "\n");

    fprintf(file, "CELL_TYPES %d\n", tetraNumber());

    for(int i=0; i<tetraNumber(); i++)
        fprintf(file, "%d ", 10);
    fprintf(file, "\n\n");

    fprintf(file, "POINT_DATA %d \n", vertexNumber());
    fprintf(file, "FIELD FieldData 1\n");

    fprintf(file, "critici 1 %d int\n", vertexNumber());

    for(int i=0; i<vertexNumber(); i++)
        fprintf(file, "%d ", critical_vertices[i]);

    fprintf(file, "\n\n");

    int segm=0;
    if(ascending!=NULL)segm+=1;
    if(descending!=NULL)segm+=1;

    fprintf(file, "CELL_DATA %d \n", tetraNumber());
    fprintf(file, "FIELD FieldData %d\n", segm);

    if(ascending!=NULL){
        fprintf(file, "ascending 1 %d int\n", tetraNumber());

        for(int i=0; i<tetraNumber(); i++)
            fprintf(file, "%d ", (*ascending)[i]);

        fprintf(file, "\n\n");
    }

    if(descending!=NULL){
        fprintf(file, "descending 1 %d int\n", tetraNumber());

        for(int i=0; i<tetraNumber(); i++)
            fprintf(file, "%d ", (*descending)[i]);

        fprintf(file, "\n\n");
    }

    fclose(file);
}

template<class Coordinate, class FieldValue>
vector<Coordinate> Mesh<Coordinate,FieldValue>::readField(char* file_name){

    vector<Coordinate> fieldvalues;
    FILE* file;
    file = fopen(file_name, "r");

    int num_values;
    fscanf(file, "%d", &num_values);
    float val;
    for(int i=0; i<num_values; i++){
        fscanf(file, "%f", &val);
        fieldvalues.push_back(val);
    }

    fclose(file);

    return fieldvalues;
}

}
