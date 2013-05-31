#ifndef FEM_H
#define FEM_H

#include <iostream>
#include <map>
#include <vector>

using namespace std;

typedef vector<double> FemVec;
typedef map<int, map<int,double> > FemMat;
typedef map<int, FemVec > FemVertices ;

class FEM{
public:
    //FEM(double _delta, double _epsilon):delta(_delta),epsilon(_epsilon){};
    explicit FEM(double _epsilon, int _nLevel):epsilon(_epsilon),nLevel(_nLevel),lambda(0){}
    ~FEM(){};

    static double delta;
    static double K2(double x,double y); //calculate square of k

    void ReadData(string filename); //read vertices and faces from a text file
    void SetDelta(double _delta){delta = _delta;}
    void SolveAll();

private:
    FemVertices vertices;
    vector< vector<int> > faces;
    FemMat A, M; //Global matrices Ah and Mh
    double epsilon;
    int nLevel;

    //These methods write results to text files
    void WriteK2(string filename);  //write values of k square to a text file
    void WriteMatrix(FemMat &mat, string filename);  //write values of global matrix Ah to a text file
    void WriteMh(string filename);  //write values of global matrix Mh to a text file
    void WriteU(string filename);   //write values of U to a text file
    void WriteLambda(string filename);   //write values of lambda to a text file

    void SplitFaces();
    int InsertVertex(FemVec &coords, FemVertices &vertexList);

    void SolveLocalMatrix();  //Solve local matrices and transfer them to global matrices A and M
    void SolveEigenProblem();

    double CalculateLambda(FemVec &u);
    double DotProduct(FemVec &x, FemVec &y);
    void CalculateResidual(FemVec &result,FemVec &f,FemMat &B, FemVec &x); //calculate result = f - B.x
    void CalculateBx(FemVec &result, FemMat &B, FemVec &x); //result = B . x
    double CalculateNorm2(FemVec &u);  //calculate square of norm of u
    void Normalize(FemVec &u);

    //These methods are for solving linear system with Conjugate Gradient method
    void UpdateVector(FemVec &result, double coefficient, FemVec &x);
    void UpdateSearchVector(FemVec &result, double coefficient, FemVec &x);
    void ConjugateGradientSolve(FemVec &result, FemVec &f);  //solve linear system of equation

    //This method is for solving linear system with GaussSeidel method in 1 dimension
    void GaussSeidelSolve(FemVec &result, FemVec &f);

    FemVec F, U, Au; //Note that F = M . u
    double lambda;
};

#endif // FEM_H
