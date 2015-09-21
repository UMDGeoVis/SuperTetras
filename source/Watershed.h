#ifndef WATERSHED_H
#define WATERSHED_H
#include "Mesh.h"
#include <cmath>
#include "assert.h"
#include <algorithm>
#include <set>
#include <list>

#define PIGRECO 3.1415927




namespace LibTetra
{


    template<class Coordinate, class FieldValue>
    class Watershed: public Mesh<Coordinate,FieldValue>
    {
    private:
        map<int,int> region_to_critic;

    public:
        inline Watershed() {this->critical_vertices=vector<int>(0);}

        inline int choose_best_vertex_label(TetraIndex t, vector<int>* v_label, bool ascending){
            int label =-1;
            float fieldValue=0;

            for(int i=0; i<4; i++){
                VertexIndex v = getTetra(t)->TV(i);
                if((*v_label)[v] <= -1) continue;

                if(ascending){
                    if(label==-1){
                        label = (*v_label)[v];
                        fieldValue = getVertex(v)->fieldValue();
                    }
                    else if(fieldValue > getVertex(v)->fieldValue()){
                        label = (*v_label)[v];
                        fieldValue = getVertex(v)->fieldValue();
                    }
                }
                else{
                    if(label==-1){
                        label = (*v_label)[v];
                        fieldValue = getVertex(v)->fieldValue();
                    }
                    else if(fieldValue < getVertex(v)->fieldValue()){
                        label = (*v_label)[v];
                        fieldValue = getVertex(v)->fieldValue();
                    }
                }
            }

            return label;
        }

        inline int choose_best_tetra_label(map<int,int>* tetra_labels, bool ascending){
            int label=-1;
            int n_labels=0;

            for(map<int,int>::iterator it = tetra_labels->begin(); it != tetra_labels->end(); it++){
                if(it->second > n_labels){
                    label = it->first;
                    n_labels = it->second;
                }
                else if(it->second == n_labels){
                    if(ascending && getVertex(region_to_critic[it->first])->fieldValue() < getVertex(region_to_critic[label])->fieldValue()){
                        label = it->first;
                        n_labels = it->second;
                    }
                    else if(!ascending && getVertex(region_to_critic[it->first])->fieldValue() > getVertex(region_to_critic[label])->fieldValue()){
                        label = it->first;
                        n_labels = it->second;
                    }
                }
            }
            assert(label != -1);
            return label;
        }

        vector<int>* compute_watershed_segmentation(bool); ///distortion computation at vertex v
    };

    struct Comparer{
        inline operator()(pair<int,float> v1, pair<int, float> v2){
            return v1.second > v2.second;
        }
    };


    template<class Coordinate, class FieldValue>
    vector<int>* Watershed<Coordinate,FieldValue> :: compute_watershed_segmentation(bool ascending){

        vector<int>* segmentation = new vector<int>(tetraNumber(), -1);
        vector<int> v_labels = vector<int>(vertexNumber(), -1);
        vector<bool> visited = vector<bool>(vertexNumber(), false);

        if(this->critical_vertices.size() == 0) this->critical_vertices = vector<int>(vertexNumber(),0);

        priority_queue<pair<int,float>, vector<pair<int,float> >, Comparer> sorted_vertices;

        if(ascending)
            for(int i=0; i<vertexNumber(); i++)
                sorted_vertices.push(pair<int,FieldValue>(i,getVertex(i)->fieldValue()));
        else
            for(int i=0; i<vertexNumber(); i++)
                sorted_vertices.push(pair<int,FieldValue>(i,-getVertex(i)->fieldValue()));


        vector<VertexIndex> vv;
        set<int> local_labels;
        queue<VertexIndex> flat_region;
        list<VertexIndex> touched_vertexes;
        bool watershed_in_area;
        while(!sorted_vertices.empty()){
            watershed_in_area = false;

            vv.clear();
            local_labels.clear();
            touched_vertexes.clear();

            VertexIndex vertice = sorted_vertices.top().first;
            sorted_vertices.pop();


            if(visited[vertice]) continue;

            flat_region.push(vertice);
            visited[vertice] = true;
            touched_vertexes.push_back(vertice);

            while(!flat_region.empty()){
                VertexIndex v = flat_region.front();
                flat_region.pop();
                vv = VV(v);


                for(int j=0; j<vv.size(); j++){
                    if(getVertex(v)->fieldValue() == getVertex(vv[j])->fieldValue() ){
                        if(!visited[vv[j]]){
                            visited[vv[j]] = true;
                            flat_region.push(vv[j]);
                            touched_vertexes.push_back(vv[j]);
                        }
                    }
                    else{
                        if(v_labels[vv[j]] > -1) local_labels.insert(v_labels[vv[j]]);
                        if(v_labels[vv[j]] == -2) watershed_in_area = true;
                    }
                }

            }

            int label=0;
            if(local_labels.size() == 0 && !watershed_in_area){
                label = vertice;
                if(ascending) this->critical_vertices[vertice] = 1;
                else          this->critical_vertices[vertice] = 2;
            }
            else if(local_labels.size() == 1) label = *(local_labels.begin());
            else label = -2;

            for(list<VertexIndex>::iterator it = touched_vertexes.begin(); it != touched_vertexes.end(); it++)
                v_labels[*it] = label;

        }

        assert(v_labels.size() == this->vertexNumber());
        list<TetraIndex> watershed_tetra;
        for(int i=0; i<tetraNumber(); i++){
            local_labels.clear();

            for(int j=0; j<4; j++){
                VertexIndex v = getTetra(i)->TV(j);
                assert(v_labels[v] != -1);
                if(v_labels[v] != -2)
                    local_labels.insert(v_labels[v]);
            }

            if(local_labels.empty()){
                (*segmentation)[i] = -1; //watershed tetrahedron
                watershed_tetra.push_back(i);
            }
            else if(local_labels.size() == 1)
                (*segmentation)[i] = (*local_labels.begin());

            else
                (*segmentation)[i] = choose_best_vertex_label(i, &v_labels, ascending);
        }


        map<int,int> regions;
        while(!watershed_tetra.empty()){
            TetraIndex t = watershed_tetra.front();
            watershed_tetra.pop_front();

            for(int j=0; j<4; j++){
                TetraIndex t1 = getTetra(t)->TT(j);
                if(t1 == noAdjacentTetraIndex || (*segmentation)[t1] <= -1) continue;

                if(regions.find((*segmentation)[t1]) == regions.end()){
                    regions[(*segmentation)[t1]] = 1;
                }
                else{
                    int reg = regions[(*segmentation)[t1]];
                    regions[(*segmentation)[t1]] = ++reg;
                }

            }


            if(regions.size() == 0) watershed_tetra.push_back(t);
            else if(regions.size() == 1) (*segmentation)[t] = (regions.begin())->first;
            else (*segmentation)[t] = choose_best_tetra_label(&regions,ascending);

            regions.clear();
        }

        return segmentation;

    }
}
#endif // WATERSHED_H
