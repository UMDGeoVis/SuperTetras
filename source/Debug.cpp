/*
 * File:   debug.cpp
 * Author: Riccardo Fellegara
 *
 * Created on 31 agosto 2009, 16.27
 */

#include <iostream>
#include "Debug.h"
#include <stdio.h>


//Debug::Debug() {
//}
//
//Debug::Debug(const Debug& orig) {
//}
//
//Debug::~Debug() {
//}

timespec Debug::time;
clock_t Debug::start;
clock_t Debug::finish;
string Debug::operation;

void Debug::Notice(string text) {
    cerr << text << endl;
}

void Debug::Begin(string op)
{
    operation = op;
    cerr << operation << flush;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time);
    start = clock();
}

void Debug::End()
{
    timespec time2;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
    finish = (double) clock();
    double delta = finish - start;
    double result = delta / static_cast<double> (CLOCKS_PER_SEC);
    fprintf(stderr, "%lf ", result);
    cerr << double(time2.tv_sec - time.tv_sec) + (time2.tv_nsec - time.tv_nsec)/1000000000.0 << endl;
}

bool Debug::checkVT(Mesh<float, float> &mesh, vector<vector<TetraIndex> > &allVTres)
{
    cout<<"check VT results.."<<endl;
    int vertexId = mesh.vertexNumber()/100;
    int start = vertexId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<TetraIndex> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(vertexId == mesh.vertexNumber())
            vertexId--;
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            for(int j = 0; j<4; j++)
            {
                if(vertexId == t->TV(j))
                {
                    ret.push_back(i);
                    break;
                }
            }
        }

        if(allVTres.at(count).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<allVTres.at(count).size()<<" "<<ret.size()<<endl;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allVTres.at(count).size();i++)
            {
                if(ret.at(j)==allVTres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<"tetra with id "<<ret.at(j)<<" not found"<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        vertexId += start;
        count++;
    }
    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkVTALL(Mesh<float, float> &mesh, vector<vector<TetraIndex> > &allVTres)
{
    cout<<"check VT results.."<<endl;
    bool allOk = true;
    for(unsigned long int vid = 0; vid<mesh.vertexNumber(); vid++)
    {
        vector<TetraIndex> ret;
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            for(int j = 0; j<4; j++)
            {
                if(vid == t->TV(j))
                {
                    ret.push_back(i);
                    break;
                }
            }
        }

        if(allVTres.at(vid).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<vid<<endl;
            cout<<allVTres.at(vid).size()<<" "<<ret.size()<<endl;
            for(unsigned int i=0;i<allVTres.at(vid).size();i++)
//                cout<<*mesh.getTetra(allVTres.at(vid).at(i))<<endl;
                cout<<allVTres.at(vid).at(i)<<" ";
            cout<<endl;
            for(unsigned int i=0;i<ret.size();i++)
//                cout<<*mesh.getTetra(ret.at(i))<<endl;
                cout<<ret.at(i)<<" ";
            cout<<endl;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allVTres.at(vid).size();i++)
            {
                if(ret.at(j)==allVTres.at(vid).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<"tetra with id "<<ret.at(j)<<" not found"<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
    }
    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkVF(Mesh<float, float> &mesh, vector<vector<Triangle> > &allVFres)
{
    cout<<"check VF results.."<<endl;
    int vertexId = mesh.vertexNumber()/100;
    int start = vertexId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<Triangle> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(vertexId == mesh.vertexNumber())
            vertexId--;
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);            
            for(int j = 0; j<4; j++)
            {
                for(int w=0;w<3;w++)
                {
                    if(vertexId == t->TF(j).FV(w) && !mesh.isIntoVector(t->TF(j),ret))
                    {
                        ret.push_back(t->TF(j));
                        break;
                    }
                }
            }
        }

        if(allVFres.at(count).size()!=ret.size())
        {
            cout<<count<<" "<<vertexId<<endl;
            cout<<"Error: the two list have different size"<<endl;
            cout<<allVFres.at(count).size()<<" "<<ret.size()<<endl;
            for(unsigned int i=0;i<allVFres.at(count).size();i++)
            {
                for(int j=0;j<3;j++)
                    cout<<allVFres.at(count).at(i).FV(j)<<" ";
                cout<<endl;
            }
            cout<<endl;
            for(unsigned int i=0;i<ret.size();i++)
            {
                for(int j=0;j<3;j++)
                    cout<<ret.at(i).FV(j)<<" ";
                cout<<endl;
            }
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allVFres.at(count).size();i++)
            {
                if(ret.at(j)==allVFres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<"tetra with id "<<ret.at(j)<<" not found"<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        vertexId += start;
        count++;
    }
    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkVE(Mesh<float, float> &mesh, vector<vector<Edge> > &allVEres)
{
    cout<<"check VE results.."<<endl;
    int vertexId = mesh.vertexNumber()/100;
    int start = vertexId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<Edge> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(vertexId == mesh.vertexNumber())
            vertexId--;
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            bool isFound = false;
            for(int i=0;i<4;i++)
            {
                for(int j=i+1;j<4;j++)
                {
                    Edge e = Edge(t->TV(i),t->TV(j));
                    if((vertexId == t->TV(i) || vertexId == t->TV(j)) && !mesh.isIntoVector(e,ret))
                    {
                        ret.push_back(e);
                        isFound = true;
                    }
                }
            }
        }

        if(allVEres.at(count).size()!=ret.size())
        {
            cout<<count<<" "<<vertexId<<endl;
            cout<<"Error: the two list have different size"<<endl;
            cout<<allVEres.at(count).size()<<" "<<ret.size()<<endl;
            for(unsigned int i=0;i<allVEres.at(count).size();i++)
            {
                for(int j=0;j<2;j++)
                    cout<<allVEres.at(count).at(i).EV(j)<<" ";
                cout<<endl;
            }
            cout<<endl;
            for(unsigned int i=0;i<ret.size();i++)
            {
                for(int j=0;j<2;j++)
                    cout<<ret.at(i).EV(j)<<" ";
                cout<<endl;
            }
            int a; cin>>a;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allVEres.at(count).size();i++)
            {
                if(ret.at(j)==allVEres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<"tetra with id "<<found<<" not found"<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        vertexId += start;
        count++;
    }
    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkVV(Mesh<float, float> &mesh, vector<vector<VertexIndex> > &allVVres)
{
    cout<<"check VV results.."<<endl;
    int vertexId = mesh.vertexNumber()/100;
    int start = vertexId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<VertexIndex> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(vertexId == mesh.vertexNumber())
            vertexId--;
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            for(int i=0;i<4;i++)
            {
                for(int j=i+1;j<4;j++)
                {
                    if(vertexId == t->TV(i) && !mesh.isIntoVector(t->TV(j),ret))
                        ret.push_back(t->TV(j));
                    if( vertexId == t->TV(j) && !mesh.isIntoVector(t->TV(i),ret))
                        ret.push_back(t->TV(i));
                }
            }
        }

        if(allVVres.at(count).size()!=ret.size())
        {
            cout<<count<<endl;
            cout<<"Error: the two list have different size"<<endl;
            cout<<allVVres.at(count).size()<<" "<<ret.size()<<endl;
            for(unsigned int i=0;i<allVVres.at(count).size();i++)
            {
                cout<<allVVres.at(count).at(i)<<" ";
            }
            cout<<endl<<endl;
            for(unsigned int i=0;i<ret.size();i++)
            {
                cout<<ret.at(i)<<" ";
            }
            cout<<endl;
            int a; cin>>a;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allVVres.at(count).size();i++)
            {
                if(ret.at(j)==allVVres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        vertexId += start;
        count++;
    }
    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkET(Mesh<float, float> &mesh, vector<vector<TetraIndex> > &allETres)
{
    cout<<"check ET results.."<<endl;
    int tetraId = mesh.tetraNumber()/100;
    int start = tetraId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<TetraIndex> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(tetraId == mesh.tetraNumber())
            tetraId--;
        Tetra* tmp = mesh.getTetra(tetraId);
        Edge* edge = tmp->TF(0).FE(1);
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            bool isFound = false;
            for(int w=0;w<4;w++)
            {
                for(int j=w+1;j<4;j++)
                {
                    Edge e = Edge(t->TV(w),t->TV(j));
                    if(*edge == e)
                    {
                        ret.push_back(i);
                        isFound = true;
                        break;
                    }
                }
                if(isFound)
                    break;
            }
        }

        if(allETres.at(count).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<allETres.at(count).size()<<" "<<ret.size()<<endl;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allETres.at(count).size();i++)
            {
                if(ret.at(j)==allETres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<ret.size()<<" "<<allETres.at(count).size()<<endl;
                for(unsigned int i=0;i<ret.size();i++)
                    cout<<ret.at(i)<<" ";
                cout<<endl;
                for(unsigned int i=0;i<allETres.at(count).size();i++)
                    cout<<allETres.at(count).at(i)<<" ";
                cout<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        tetraId += start;
        count++;
    }

    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkEF(Mesh<float, float> &mesh, vector<vector<Triangle> > &allEFres)
{
    cout<<"check EF results.."<<endl;
    int tetraId = mesh.tetraNumber()/100;
    int start = tetraId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<Triangle> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(tetraId == mesh.tetraNumber())
            tetraId--;
        Tetra* tmp = mesh.getTetra(tetraId);
        Edge* e = tmp->TF(0).FE(1);
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            for(int w=0;w<4;w++)
            {
                for(int j=0;j<3;j++)
                {
                    if(*t->TF(w).FE(j)==*e && !mesh.isIntoVector(t->TF(w),ret))
                    {
                        ret.push_back(t->TF(w));
                        break;
                    }
                }
            }
        }

        if(allEFres.at(count).size()!=ret.size())
        {
            cout<<count<<endl;
            cout<<"Error: the two list have different size"<<endl;
            cout<<allEFres.at(count).size()<<" "<<ret.size()<<endl;
            cout<<"edge: "<<*e<<endl;
            cout<<"ret:"<<endl;
            for(unsigned int i=0;i<ret.size();i++)
                cout<<ret.at(i)<<endl;
            cout<<endl<<"allEFres:";
            for(unsigned int i=0;i<allEFres.at(count).size();i++)
                cout<<allEFres.at(count).at(i)<<endl;
            cout<<endl<<"VF debug:";
            vector<Triangle> all_faces = mesh.VF(e->EV(0));
            for(unsigned int i=0;i<all_faces.size();i++)
                cout<<all_faces.at(i)<<endl;
            allOk = false;
            int a; cin>>a;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allEFres.at(count).size();i++)
            {
                if(ret.at(j)==allEFres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<ret.size()<<" "<<allEFres.at(count).size()<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        tetraId += start;
        count++;
    }

    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkEE(Mesh<float, float> &mesh, vector<vector<Edge> > &allEEres)
{
    cout<<"check EE results.."<<endl;
    int tetraId = mesh.tetraNumber()/100;
    int start = tetraId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<Edge> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(tetraId == mesh.tetraNumber())
            tetraId--;
        Tetra* tmp = mesh.getTetra(tetraId);
        Edge* e = tmp->TF(0).FE(1);
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            bool isFound = false;
            for(int i=0;i<4;i++)
            {
                for(int j=i+1;j<4;j++)
                {
                    Edge edge = Edge(t->TV(i),t->TV(j));
                    if(*e != edge)
                    {
                        for(int w=0;w<2;w++)
                        {
                            if(edge.EV(w)==e->EV(0) && !mesh.isIntoVector(edge,ret))
                            {
                                ret.push_back(edge);
                                isFound = true;
                            }
                            if(edge.EV(w)==e->EV(1) && !mesh.isIntoVector(edge,ret))
                            {
                                ret.push_back(edge);
                                isFound = true;
                            }
                        }
                    }
                }
            }
        }

        if(allEEres.at(count).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<allEEres.at(count).size()<<" "<<ret.size()<<endl;
            allOk = false;
            int a; cin>>a;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allEEres.at(count).size();i++)
            {
                if(ret.at(j)==allEEres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<ret.size()<<" "<<allEEres.at(count).size()<<endl;
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        tetraId += start;
        count++;
    }

    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkFT(Mesh<float, float> &mesh, vector<vector<TetraIndex> > &allFTres)
{
    cout<<"check FT results.."<<endl;
    int tetraId = mesh.tetraNumber()/100;
    int start = tetraId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<TetraIndex> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(tetraId == mesh.tetraNumber())
            tetraId--;
        Tetra* tmp = mesh.getTetra(tetraId);
        Triangle face = tmp->TF(1);
        for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
        {
            Tetra* t = mesh.getTetra(i);
            for(int w=0;w<4;w++)
            {
                if(t->TF(w)==face)
                {
                    ret.push_back(i);
                    break;
                }
            }
            if(ret.size()==2)
                break;
        }

        if(allFTres.at(count).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<allFTres.at(count).size()<<" "<<ret.size()<<endl;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allFTres.at(count).size();i++)
            {
                if(ret.at(j)==allFTres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                cout<<"face: "<<face.FV(0)<<" "<<face.FV(1)<<" "<<face.FV(2)<<endl;
                cout<<"ret:"<<endl;
                for(unsigned int i=0;i<ret.size();i++)
                    cout<<ret.at(i)<<" "<<*mesh.getTetra(ret.at(i))<<" ";
                cout<<endl<<"allFTres:"<<endl;
                for(unsigned int i=0;i<allFTres.at(count).size();i++)
                    cout<<allFTres.at(count).at(i)<<" "<<*mesh.getTetra(allFTres.at(count).at(i))<<" ";
                cout<<endl<<"vt di debug: ";
                VertexIndex v0 = face.FV(0);
                cout<<"id vertice: "<<v0<<endl;
                bool* isBorder = new bool();
                for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
                    mesh.getTetra(i)->setStatus(unchecked);
                vector<TetraIndex> all_tetra = mesh.VT(v0,isBorder);
                for(unsigned int i=0;i<all_tetra.size();i++)
                {
                    cout<<all_tetra.at(i)<<" ";
                }
                cout<<endl<<"vt (versione bruteforce):"<<endl;
                vector<TetraIndex> ret2;
                for(unsigned int i = 0; i < mesh.tetraNumber(); i++)
                {
                    Tetra* t = mesh.getTetra(i);
                    for(int j = 0; j<4; j++)
                    {
                        if(v0 == t->TV(j))
                        {
                            ret2.push_back(i);
                            break;
                        }
                    }
                }
                for(unsigned int i=0;i<ret2.size();i++)
                {
                    cout<<ret2.at(i)<<" ";
                }
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        tetraId += start;
        count++;
    }

    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}

bool Debug::checkFF(Mesh<float, float> &mesh, vector<vector<Triangle> > &allFFres)
{
    cout<<"check FF results.."<<endl;
    int tetraId = mesh.tetraNumber()/100;
    int start = tetraId;
    int count = 0;
    bool allOk = true;
    while(count < 100)
    {
        vector<Triangle> ret;
        //ATTENZIONE SE TROVO UN NUMERO DI VERTICI/TETRA DIVISIBILE PER 100
        //ALL'ULTIMO GIRO SI PIANTA PERCHE' SFORO IL NUMERO DI VERTICI/TETRA DI UNO
        //TAPULLO X ORA
        if(tetraId == mesh.tetraNumber())
            tetraId--;
        Tetra* tmp = mesh.getTetra(tetraId);
        Triangle face = tmp->TF(b);
        for(TetraIndex i = a; i < mesh.tetraNumber(); i=(TetraIndex)(i+1))
        {
            Tetra* t = mesh.getTetra(i);
            for(PositionIndex w=a;w<=d;w=(PositionIndex)(w+1))
            {
                Triangle f = t->TF(w);
                if(f != face && !mesh.isIntoVector(f,ret))
                {
                    bool isFound =false;
                    for(PositionIndex j=a;j<d;j=(PositionIndex)(j+1))
                    {
                        for(PositionIndex x=a;x<d;x=(PositionIndex)(x+1))
                        {
                            if( *f.FE(j) == *face.FE(x) )
                            {
                                ret.push_back(f);
                                isFound = true;
                                break;
                            }
                        }
                        if(isFound)
                            break;
                    }
                }
            }
        }

        if(allFFres.at(count).size()!=ret.size())
        {
            cout<<"Error: the two list have different size"<<endl;
            cout<<allFFres.at(count).size()<<" "<<ret.size()<<endl;
            cout<<"ret:"<<endl;
            for(unsigned int i=0;i<ret.size();i++)
            {
                cout<<ret.at(i)<<" ";
                if(i%3==0)
                    cout<<endl;
            }
            cout<<endl<<"allFTres:"<<endl;
            for(unsigned int i=0;i<allFFres.at(count).size();i++)
            {
                cout<<allFFres.at(count).at(i)<<" ";
                if(i%3==0)
                    cout<<endl;
            }
            int a; cin>>a;
        }

        for(unsigned int j=0;j<ret.size();j++)
        {
            bool found = false;
            for(unsigned int i=0;i<allFFres.at(count).size();i++)
            {
                if(ret.at(j)==allFFres.at(count).at(i))
                {
                    found=true;
                    break;
                }
            }
            if(!found)
            {
                allOk = false;
                int a; cin>>a;
            }
            found=false;
        }
        ret.clear();
        tetraId += start;
        count++;
    }

    if(allOk)
    {
        cout<<"   OK"<<endl;
        return true;
    }
    else
    {
        cout<<"   Error"<<endl;
        return false;
    }
}
