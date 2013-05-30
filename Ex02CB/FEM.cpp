#include "FEM.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Colsamm/Colsamm.h"

using namespace std;
using namespace ::_COLSAMM_;

double FEM::delta = 0;

double FEM::DotProduct(const vector<double> &x, const vector<double> &y){
    if (x.size() != y.size()){
        cerr<<"Matrix size mismatch when calculating dot product"<<endl;
        return 0.0;
    }
    double result = 0.0;
    for (unsigned int i=0; i< x.size(); ++i)
        result += x[i]*y[i];
    return result;
}

double FEM::CalculateNorm2(const vector<double> &u){
    return DotProduct(u,u); //square of norm of u
}

void FEM::CalculateBx(vector<double> &result,  FemMat &B, const vector<double> &x){
    for (unsigned int row=0; row<x.size(); ++row){
        FemMat::iterator i = B.find(row); //check row is present in Mh
        if (i != B.end()){  //yes row is present in Mh
            double sum = 0.0;
            for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
                sum += j->second * x[j->first];  //multiply row of Mh with u
            result[row] = sum;
        }
        else{ //row k is not present in Mh, i.e. this row of Mh is zero
            result[row] = 0.0;
        }
    }
}

void FEM::CalculateResidual(vector<double> &result, const vector<double> &f, FemMat &B, const vector<double> &x){
    for (unsigned int row=0; row<x.size(); ++row){
        FemMat::iterator i = B.find(row); //check row is present in Mh
        if (i != B.end()){  //yes row is present in Ah
            double sum = 0.0;
            for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
                sum += j->second * x[j->first];  //multiply row of Mh with u
            result[row] = f[row] - sum;
        }
        else{ //row k is not present in Mh, i.e. this row of Mh is zero
            result[row] = f[row];
        }
    }
}

void FEM::UpdateVector(vector<double> &result, double coefficient, vector<double> &x){
    for (unsigned int i=0; i<x.size(); ++i){
        result[i] += coefficient * x[i];
    }
}

void FEM::UpdateSearchVector(vector<double> &result, double coefficient, vector<double> &x){
    for (unsigned int i=0; i<x.size(); ++i){
        result[i] = x[i] + coefficient * result[i];
    }
}

void FEM::ConjugateGradientSolve(vector<double> &result, const vector<double> &f){
    vector<double> Ad = vector<double>(result.size(),0.0);
    vector<double> r = vector<double>(result.size(),0.0);
    CalculateResidual(r,f,A,result);   //first residual
    vector<double> d = r;   //first search vector
    double normr2 = CalculateNorm2(r);
    double epsilon2 = epsilon*epsilon;
    while (normr2>epsilon2){
        CalculateBx(Ad,A,d); //find Ad from d
        double alpha = normr2/DotProduct(d,Ad);
        UpdateVector(result,alpha,d); //update x from alpha and d
        UpdateVector(r,-alpha,Ad); //update r from alpha and Ad
        double normr2next = CalculateNorm2(r);
        UpdateSearchVector(d,normr2next/normr2,r);
        normr2 = normr2next;
        //cout<<normr2<<endl;
    }
}

void FEM::GaussSeidelSolve(vector<double> &result, const vector<double> &f){
    vector<double> Ad = vector<double>(result.size(),0.0);
    vector<double> r = vector<double>(result.size(),0.0);
    CalculateResidual(r,f,A,result);   //first residual
    vector<double> d = r;   //first search vector
    double normr2 = CalculateNorm2(r);
    double epsilon2 = epsilon*epsilon;
    while (normr2>epsilon2){
        for (unsigned int row=0; row<result.size(); ++row){
            FemMat::iterator i = A.find(row); //check row is present in Mh
            if (i != A.end()){  //yes row is present in Ah
                double sum = 0.0;
                double aii = 0.0;
                for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
                    if (j->first != row)
                        sum += j->second * result[j->first];  //multiply row of Mh with u
                    else
                        aii = j->second;
                result[row] = (f[row] - sum) / aii;
            }
            else{ //row k is not present in Mh, i.e. this row of Mh is zero
                result[row] = 0;
            }
        }
        CalculateResidual(r,f,A,result);
        normr2 = CalculateNorm2(r);
    }
}

void FEM::Normalize(vector<double> &u){
    double normu = sqrt(CalculateNorm2(u));
    for (unsigned int i=0; i<u.size(); ++i)
        u[i] /= normu;
}

double FEM::CalculateLambda(const vector<double> &u){
    CalculateBx(Au,A,u);
    double numerator = DotProduct(u,Au); //Calculate numerator
    double denominator = DotProduct(u,F); //Calculate denominator
    return numerator/denominator;
}

void FEM::SolveEigenProblem(){
    lambda = 0.0;
    double lambdaOld = 1.0;
    U = vector<double>(vertices.size(),1.0);
    F = vector<double>(vertices.size(),0.0);
    Au = vector<double>(vertices.size(),0.0);
    int iter = 0;
    while (fabs((lambda-lambdaOld)/lambdaOld)>pow((double)10.00,-10)){
        lambdaOld = lambda;
        CalculateBx(F,M,U);
        ConjugateGradientSolve(U,F);   //solve linear system of equations to find U
        //GaussSeidelSolve(U,F);   //solve linear system of equations to find U
        Normalize(U);
        lambda = CalculateLambda(U);
        cout<<"lambda at iter "<<iter<<" = "<<lambda<<endl;
        ++iter;
    }
}

double FEM::K2(double x,double y){
    return (100+delta)* exp(-50*(x*x+y*y))-100;
}

void FEM::SolveLocalMatrix(){
   //  function<double<const FEM&, double,double)> GetK2 = &FEM::K2;
    for (unsigned int i=0; i< faces.size(); ++i){
        //find local matrix for each face
        ELEMENTS::Triangle my_element;
        vector< vector< double > > my_local_matrix;
        vector<double> corners(6, 0.0);
        // array corners contains the x- and y-coordinates of the
        // triangle corners in the order x0, y0, x1, y1, x2, y2
        int vertex0 = faces[i][0];
        int vertex1 = faces[i][1];
        int vertex2 = faces[i][2];

        corners[0] = vertices[vertex0][0];
        corners[1] = vertices[vertex0][1];
        corners[2] = vertices[vertex1][0];
        corners[3] = vertices[vertex1][1];
        corners[4] = vertices[vertex2][0];
        corners[5] = vertices[vertex2][1];
        // pass the corners to the finite element
        my_element(corners);

        //Computer local matrix A
        my_local_matrix =
        my_element.integrate(grad(v_()) * grad(w_())-func<double>(K2)*v_()*w_());
        //transfering value from local matrix A to global matrix Ah
        A[vertex0][vertex0] += my_local_matrix[0][0];
        A[vertex0][vertex1] += my_local_matrix[0][1];
        A[vertex0][vertex2] += my_local_matrix[0][2];
        A[vertex1][vertex0] += my_local_matrix[1][0];
        A[vertex1][vertex1] += my_local_matrix[1][1];
        A[vertex1][vertex2] += my_local_matrix[1][2];
        A[vertex2][vertex0] += my_local_matrix[2][0];
        A[vertex2][vertex1] += my_local_matrix[2][1];
        A[vertex2][vertex2] += my_local_matrix[2][2];

        //compute local matrix M
        my_local_matrix =
        my_element.integrate(v_()*w_());
        //transfering value from local matrix M to global matrix Mh
        M[vertex0][vertex0] += my_local_matrix[0][0];
        M[vertex0][vertex1] += my_local_matrix[0][1];
        M[vertex0][vertex2] += my_local_matrix[0][2];
        M[vertex1][vertex0] += my_local_matrix[1][0];
        M[vertex1][vertex1] += my_local_matrix[1][1];
        M[vertex1][vertex2] += my_local_matrix[1][2];
        M[vertex2][vertex0] += my_local_matrix[2][0];
        M[vertex2][vertex1] += my_local_matrix[2][1];
        M[vertex2][vertex2] += my_local_matrix[2][2];
    }
}

void FEM::WriteMatrix(FemMat &mat, string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (FemMat::iterator i = mat.begin(); i != mat.end(); ++i){
            for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j){
                outfile<<i->first<<" "<<j->first<<" "<<j->second<<endl;
            }
        }
        outfile.close();
    }
}

void FEM::SolveAll(){
    WriteK2("ksq.txt");
    SolveLocalMatrix();
    WriteMatrix(A,"A.txt");
    WriteMatrix(M,"M.txt");
    SolveEigenProblem();
    WriteU("eigenmode.txt");
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
        int index = 0;
        double x = 0.0;
        double y = 0.0;
        ss>>index>>x>>y;
        if (infile.fail()){
            cerr<<"error in reading line "<<lineNum<<endl;
            return;
        }
        vector<double> point;
        point.push_back(x);
        point.push_back(y);
        vertices[index] = point;
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
        int vertex0 = 0;
        int vertex1 = 0;
        int vertex2 = 0;
        ss>>vertex0>>vertex1>>vertex2;
        if (ss.fail()){
            cerr<<"error in reading line "<<lineNum<<endl;
            return;
        }
        vector<int> face;
        face.push_back(vertex0);
        face.push_back(vertex1);
        face.push_back(vertex2);
        faces.push_back(face);
    }

    infile.close();
}

void FEM::WriteK2(string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (FemVertices::iterator it = vertices.begin(); it != vertices.end(); ++it){
            double val = K2(it->second[0], it->second[1]);
            outfile<<it->second[0]<<" "<<it->second[1]<<" "<<val<<endl;
        }
        outfile.close();
    }
}

void FEM::WriteU(string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (unsigned int i=0; i< U.size(); ++i){
            FemVertices::iterator vertexIt = vertices.find(i);
            if (vertexIt != vertices.end())
                outfile<<vertexIt->second[0]<<" "<<vertexIt->second[1]<<" "<<U[i]<<endl;
        }
        outfile.close();
    }
}
