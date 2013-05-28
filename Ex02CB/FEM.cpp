
#include "FEM.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include "Colsamm\Colsamm.h"

using namespace std;
using namespace ::_COLSAMM_;

double FEM::K2(double x,double y){
    return (100+delta)* exp(-50*(x*x+y*y))-100;
}

void FEM::ReadData(string filename){
    ifstream infile(filename.c_str());
    string line;
    int lineNum = 0;
    istringstream ss;

    //Read data of vertices
    getline(infile,line); lineNum = 1;
    ss.clear(); ss.str(line);
    int numVertices;
    ss>>numVertices;
    if (ss.fail()){
        cerr<<"error in reading number of vertices at line "<<lineNum<<endl;
        return;
    }
    getline(infile,line); ++lineNum;  //this line will be ignored
    for (int i=0;i<numVertices; ++i){
        getline(infile,line); ++lineNum;
        ss.clear(); ss.str(line);
        unsigned int index = 0;
        double x = 0.0;
        double y = 0.0;
        ss>>index>>x>>y;
        if (infile.fail()){
            cerr<<"error in reading line "<<lineNum<<endl;
            return;
        }
        Coordinates point(x,y);
        vertices.insert(pair<unsigned int, Coordinates>(index,point));
    }

    //Read data of faces
    getline(infile,line); ++lineNum;
    ss.clear();ss.str(line);
    int numFaces;
    ss>>numFaces;
    if (ss.fail()){
        cerr<<"error in reading number of faces at line "<<lineNum<<endl;
        return;
    }
    getline(infile,line); ++lineNum;//this line will be ignored
    for (int i=0;i<numFaces; ++i){
        getline(infile,line); ++lineNum;
        ss.clear();ss.str(line);
        unsigned int vertex1 = 0;
        unsigned int vertex2 = 0;
        unsigned int vertex3 = 0;
        ss>>vertex1>>vertex2>>vertex3;
        if (ss.fail()){
            cerr<<"error in reading line "<<lineNum<<endl;
            return;
        }
        vector<unsigned int> face;
        face.push_back(vertex1);
        face.push_back(vertex2);
        face.push_back(vertex3);
        faces.push_back(face);
    }

    infile.close();
}

void FEM::WriteK2(string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (map<unsigned int, Coordinates>::iterator it=vertices.begin(); it!=vertices.end(); ++it){
            double val = K2(it->second.x, it->second.y);
            outfile<<it->second.x<<" "<<it->second.y<<" "<<val<<endl;
        }
        outfile.close();
    }
}
