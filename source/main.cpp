#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <QString>
#include <QStringList>
#include "LibraryHeader.h"
#include "Debug.h"

#define BOLD  "\033[1m\033[33m" //for dark background shell
//#define BOLD "\033[1m\033[31m"  //for white background shell
#define RESET   "\033[0m"


void print_help();
void print_paragraph(string stringa, int cols);

int main (int argc, char* argv[]){

    if(argc == 1){
        print_help();
        return 0;
    }

    FILE* file_mesh;
    /// load mesh in .ts format

    cout<<"Before initialize class"<<endl;
    SuperTetras<double, double>* mesh = new SuperTetras<double, double>();
    cout<<"After"<<endl;

    if(argc < 3) printf("Input Error [not enough argument]: please type ./distortion for help\n");
    int regions = 100; //default values
    double weightF = 0.2;
    int normStrategy = 0; //default
    char* outfile = "dummy.vtk";
    char* boundariesFile = "dummy_boundaries.vtk";
    bool addIter = false;

    for(int i=1; i<(argc); i += 2){
        char* tag = argv[i];
        if(strcmp(tag, "-i") == 0){
            file_mesh = fopen(argv[i+1], "r");
            mesh->loadTS(file_mesh);
            mesh->build();
            printf("%d %d\n", mesh->vertexNumber(), mesh->tetraNumber());
        }

        else if(strcmp(tag, "-nSeg")==0){
            regions = atoi(argv[i+1]);
            cout<<"Regions: "<<regions<<endl;
        }

        else if(strcmp(tag, "-weight")==0){
            weightF = atof(argv[i+1]);
            cout<<"Weight: "<<weightF<<endl;
        }

        else if(strcmp(tag, "-norm")==0){
            normStrategy = atoi(argv[i+1]);
        }

        else if(strcmp(tag, "-vtk")==0){
            outfile = argv[i+1];
            cout<<"File output "<<outfile<<endl;
        }

        else if(strcmp(tag, "-vtkBound")==0){
            boundariesFile = argv[i+1];
            cout<<"File boundaries "<<boundariesFile<<endl;
        }

        else if(strcmp(tag, "-iter")==0){
            addIter=true;
            cout<<"Iteratively adding mode"<<endl;
        }
    }

    mesh->setNDesiredRegions(regions);
    mesh->setWeightOfField(weightF);
    if(normStrategy==0 || normStrategy==1)
        mesh->setNormStrategy(normStrategy);
    else
        mesh->setNormStrategy(0);
    mesh->initStructs();

    cout<<"Parameters: regions="<<mesh->getNDesiredRegions()<<" and weight="<<mesh->getWeightOfField()<<endl;
    //if(addIter)
      //  cout<<"Iterative region adding mode on"<<endl;
    //mesh->setAddIterative(addIter);

    mesh->triggerSegmentation();

    mesh->writeSuperTetrasVTK(outfile);
    mesh->writeSupertetrasBoundaries(boundariesFile);

    QString qs = QString::fromUtf8(outfile);
    QStringList qsl = qs.split(".");
    qs = qsl.at(0);
    qs.append(".dat");
    cout<<"Out "<<qs.toStdString().c_str()<<endl;
    QString qsegm=qsl.at(0);
    qsegm.append("_segm.vtk");

    mesh->output_segmentation(qs.toStdString().c_str(), mesh->countReg());
    qs=qsl.at(0);
    qs.append("_field.dat");
    mesh->output_segAvgField(qs.toStdString().c_str());
    //mesh->writeGraphCutResult(qsegm.toStdString().c_str());

    return 0;
}

void print_help(){

    //annoying stuff to get the dimension of the output shell (!!! not sure it works on Mac,
    //everything based on the command tput. If it doesn't work the dimension si setted to 80 by default)
    FILE* fp;
    char path[1035];

    int cols;
    fp = popen("tput cols", "r");
    if(fp != NULL){
      fgets(path, sizeof(path)-1, fp);
      cols = atoi(path);
    }
    else{
      cols = 80;
    }

    //printf start
      printf(BOLD "\n  NAME:\n\n" RESET);
      printf("\tSimulated Immersion for 3D mesh \n\n" RESET);

      printf(BOLD "  USAGE: \n\n" RESET);
      printf(BOLD "    -i [input file]\n\n" RESET);
      print_paragraph("the input file in .ts format", cols);

      printf(BOLD "    -f [field input file]\n\n" RESET);
      print_paragraph("the input file for the field values on which compute the descending or ascending Morse complexes. First value of the [field input file] indicates the number of vertexes and the remaining are the field values for each vertex", cols);

      printf(BOLD "    -descending [output file]\n\n" RESET);
      print_paragraph("compute the segmentation coresponding to the descending Morse complex outputting the result in the [output file]. First entry of the file is the set of faces, remaining entries are integers corresponding to the belonging region of each face",cols);

      printf(BOLD "    -ascending [output file]\n\n" RESET);
      print_paragraph("compute the segmentation coresponding to the ascending Morse complex outputting the result in the [output file]. First entry of the file is the set of faces, remaining entries are integers corresponding to the belonging region of each face",cols);

      printf(BOLD, "   -numS [number]\n\n" RESET);
      print_paragraph("set the number of desired regions for the segmentation", cols);

      printf(BOLD, "   -weight [w]\n\n" RESET);
      print_paragraph("set the parameter which states the importance of the field value respect to the spatial distance", cols);

      printf(BOLD, "   -norm [ns]\n\n" RESET);
      print_paragraph("Use range (ns=0) or standard-deviation-based (ns=1) strategy to normalize field value of the mesh", cols);

      printf(BOLD "    -vtk [output file]\n\n" RESET);
      print_paragraph("write a file [output file] in .vtk formant with all the segmentations computed", cols);

      printf(BOLD "  EXAMPLE: \n\n" RESET);
      printf("          .\\SimulatedImmersion -i mesh.ts -descending mesh_desc.dat -ascending mesh_asc.dat -vtk mesh.vtk\n\n");
      print_paragraph("read as input file the mesh [mesh.ts] computes both the descending and ascending Morse complexes outputting the results in a vtk file", cols);


      printf(BOLD "  IMPLEMENTATION:\n\n" RESET);
      printf("          Author: Ciccio\n");
      printf("          Group: G3 Geometry and Graphics Group\n");
      printf("          Last Update: 08/10/2012\n\n");

      printf(BOLD "  DESCRIPTION: \n\n" RESET);
      print_paragraph("Implementation of the watershed based algorithm for simulated immersion (3D)", cols);
}

void print_paragraph(string stringa, int cols){
    if(stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}


