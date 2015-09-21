#ifndef REGIONSTRUCTURE_H
#define REGIONSTRUCTURE_H

#include "Vertex.h"
#include <tr1/unordered_set>
#include <vector>

using namespace std;

namespace LibTetra {

template<class Coordinate, class FieldValue> class RegionStructure{

private:
    std::tr1::unordered_set<TetraIndex> tetras;
    int regIndex;
    Vertex<Coordinate, FieldValue> regCenter;
    FieldValue regColor;


public:
    inline std::tr1::unordered_set<TetraIndex> getTetras(){
        return tetras;
    }

    inline TetraIndex getNthTetra(int index){
        return tetras[index];
    }

    inline bool isInRegion(TetraIndex TI){

        if(tetras.count(TI)>0)
            return true;
        return false;
    }

    inline void setCentroid(Vertex<Coordinate, FieldValue> V){
        regCenter=V;
    }

    inline void setIndex(int R){
        regIndex=R;
    }
};
}

#endif // REGIONSTRUCTURE_H
