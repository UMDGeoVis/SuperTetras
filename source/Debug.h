/*
 * File:   debug.h
 * Author: Riccardo Fellegara
 *
 * Created on 31 agosto 2009, 16.27
 */

#ifndef _DEBUG_H
#define	_DEBUG_H
#include <sys/time.h>
#include <time.h>
#include <string>
#include "LibraryHeader.h"
using namespace std;
using namespace LibTetra;

///A static class used for debugging
class Debug {
public:
    //Debug();
    //Debug(const Debug& orig);
    //virtual ~Debug();
    ///A public static method that writes to standard output the string passed as argument
    /*!
     * \param text a string argument, represents the string to write
     */
    static void Notice(string text);
    ///A public static method that writes to standard output the string passed as argument, with time information
    /*!
     * \param operation a string argument, represents the string to write
     */
    static void Begin(string operation);
    ///A public static method that writes to standard output the time information of an ended operation
    static void End();
    ///
    static bool checkVT(Mesh<float,float>& mesh, vector< vector<TetraIndex> >& allVTres);
     static bool checkVTALL(Mesh<float,float>& mesh, vector< vector<TetraIndex> >& allVTres);
    ///
    static bool checkVF(Mesh<float,float>& mesh, vector< vector<Triangle> >& allVFres);
    ///
    static bool checkVE(Mesh<float,float>& mesh, vector< vector<Edge> >& allVEres);
    ///
    static bool checkVV(Mesh<float,float>& mesh, vector< vector<VertexIndex> >& allVVres);
    ///
    static bool checkET(Mesh<float,float>& mesh, vector< vector<TetraIndex> >& allETres);
    ///
    static bool checkEF(Mesh<float,float>& mesh, vector< vector<Triangle> >& allEFres);
    ///
    static bool checkEE(Mesh<float,float>& mesh, vector< vector<Edge> >& allEEres);
    ///
    static bool checkFT(Mesh<float,float>& mesh, vector< vector<TetraIndex> >& allFTres);
    ///
    static bool checkFF(Mesh<float,float>& mesh, vector< vector<Triangle> >& allFFres);
private:
    ///A static timeval variable
    static timespec time;
    ///A static clock_t variable
    static clock_t start;
    ///A static clock_t variable
    static clock_t finish;
    ///A static string variable
    static string operation;

};

#endif	/* _DEBUG_H */

