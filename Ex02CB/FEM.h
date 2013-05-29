#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <map>
#include <vector>

using namespace std;

struct Coordinates{
    Coordinates():x(0.0),y(0.0){};
    Coordinates(double _x, double _y):x(_x),y(_y){};
    double x, y;
};

struct GlobalMatrices{
    GlobalMatrices():A(0.0),M(0.0){};
    GlobalMatrices(double _A, double _M):A(_A),M(_M){};
    double A, M;
};

class FEM{
public:
    //FEM(double _delta, double _epsilon):delta(_delta),epsilon(_epsilon){};
    FEM(double _epsilon):epsilon(_epsilon){};
    ~FEM(){};

    static double delta;
    static double K2(double x,double y); //calculate square of k

    void SetUpNeighborLists();
    void SetDelta(double _delta){delta = _delta;}
    void ReadData(string filename); //read vertices and faces from a text file
    void WriteK2(string filename);  //write values of k square to a text file
    void WriteAh(string filename);
    void WriteMh(string filename);
    void Solve();                   //calculate
    void SolveEigenProblem();
private:
    map<int, Coordinates> vertices;
    vector< vector<int> > faces;
    map<int, map<int,GlobalMatrices> > AM; //Ah[i][j]

    double epsilon;
};

#endif // FEM_H
