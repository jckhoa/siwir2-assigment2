#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct Coordinates{
    Coordinates(double _x, double _y):x(_x),y(_y){};
    double x,y;
};

class FEM{
public:
    FEM(double _delta, double _epsilon):delta(_delta),epsilon(_epsilon){};
    ~FEM(){};

    double K2(double x,double y); //calculate square of k
    void SetUpNeighborLists();    //???
    void ReadData(string filename); //read vertices and faces from a text file
    void WriteK2(string filename);  //write values of k square to a text file
    void WriteAh(string filename);
    void WriteMh(string filename);

private:
    map<unsigned int, Coordinates> vertices;
    vector< vector<unsigned int> > faces;

    double delta;
    double epsilon;
};

#endif // FEM_H
