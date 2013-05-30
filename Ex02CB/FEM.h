#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <map>
#include <vector>

using namespace std;

typedef map<int, map<int,double> > FemMat ;
typedef map<int, vector<double> > FemVertices ;
class FEM{
public:
    //FEM(double _delta, double _epsilon):delta(_delta),epsilon(_epsilon){};
    FEM(double _epsilon):epsilon(_epsilon){};
    ~FEM(){};

    static double delta;
    static double K2(double x,double y); //calculate square of k

    void SetDelta(double _delta){delta = _delta;}
    void ReadData(string filename); //read vertices and faces from a text file
    void SolveAll();
private:
    FemVertices vertices;
    vector< vector<int> > faces;
    FemMat A, M; //Global matrices Ah and Mh
    double epsilon;

    void WriteK2(string filename);  //write values of k square to a text file
    void WriteMatrix(FemMat &mat, string filename);  //write values of global matrix Ah to a text file
    void WriteMh(string filename);  //write values of global matrix Mh to a text file
    void WriteU(string filename);   //write values of U to a text file
    void SolveLocalMatrix();                   //calculate
    void SolveEigenProblem();
    double DotProduct(const vector<double> &x, const vector<double> &y);
    void CalculateResidual(vector<double> &result,const vector<double> &f,FemMat &B, const vector<double> &x); //calculate result = f - B.x
    void CalculateBx(vector<double> &result, FemMat &B, const vector<double> &x); //result = B . x
    double CalculateNorm2(const vector<double> &u);  //calculate square of norm of u
    void Normalize(vector<double> &u);
    double CalculateLambda(const vector<double> &u);

    //These methods are for solving linear system with Conjugate Gradient method
    void UpdateVector(vector<double> &result, double coefficient, vector<double> &x);
    void UpdateSearchVector(vector<double> &result, double coefficient, vector<double> &x);
    void ConjugateGradientSolve(vector<double> &result, const vector<double> &f);  //solve linear system of equation

    //This method is for solving linear system with GaussSeidel method in 1 dimension
    void GaussSeidelSolve(vector<double> &result, const vector<double> &f);

    vector<double> F, U, Au; //Note that F = M . u
    double lambda;
};

#endif // FEM_H
