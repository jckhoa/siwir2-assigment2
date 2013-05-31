#include "FEM.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Colsamm/Colsamm.h"
#include <iomanip>

using namespace std;
using namespace ::_COLSAMM_;

double FEM::delta = 0;

int FEM::InsertVertex(FemVec &coords, FemVertices &vertexList){
    //this assumes that the maximum index in the vertexList is always (vertexList.size()-1)
    int numVertices = vertexList.size();
    FemVertices::iterator it = vertexList.begin();
    for (; it != vertexList.end(); ++it){
        if (it->second[0] == coords[0] && it->second[1] == coords[1])
            break;
    }
    if (it == vertexList.end()){ //this pair of coordinates is not in the vertexList
        vertexList[numVertices] = coords;  //add it to the vertexList with index = numVertices
        return numVertices;
    }
    else return it->first;
}

void FEM::SplitFaces(){
    for (int lvl = 0; lvl < nLevel; ++lvl){ //Split at level lvl
        int numFaces = faces.size();
        for (int iface=0; iface<numFaces; ++iface){
            int vertex0 = faces[iface][0];
            int vertex1 = faces[iface][1];
            int vertex2 = faces[iface][2];

            FemVec &coords0 = vertices[vertex0];
            FemVec &coords1 = vertices[vertex1];
            FemVec &coords2 = vertices[vertex2];

            FemVec coords01, coords02, coords12; //coordinates of middle points of edges (0,1), (0,2), (1,2)
            coords01.push_back((coords0[0]+coords1[0])/2);
            coords01.push_back((coords0[1]+coords1[1])/2);

            coords02.push_back((coords0[0]+coords2[0])/2);
            coords02.push_back((coords0[1]+coords2[1])/2);

            coords12.push_back((coords1[0]+coords2[0])/2);
            coords12.push_back((coords1[1]+coords2[1])/2);

            //insert coord01, coord02, coord12 into vertices
            int vertex01 = InsertVertex(coords01,vertices);
            int vertex02 = InsertVertex(coords02,vertices);
            int vertex12 = InsertVertex(coords12,vertices);

            //replace the current face by the top face
            faces[iface][1] = vertex01;
            faces[iface][2] = vertex02;

            //add the other 3 faces
            vector<int> face1, face2, face3;
            face1.push_back(vertex01);
            face1.push_back(vertex1);
            face1.push_back(vertex12);

            face2.push_back(vertex01);
            face2.push_back(vertex02);
            face2.push_back(vertex12);

            face3.push_back(vertex02);
            face3.push_back(vertex2);
            face3.push_back(vertex12);

            faces.push_back(face1);
            faces.push_back(face2);
            faces.push_back(face3);
        }
    }
}
double FEM::DotProduct(FemVec &x, FemVec &y){
    if (x.size() != y.size()){
        cerr<<"Matrix size mismatch when calculating dot product"<<endl;
        return 0.0;
    }
    double result = 0.0;
    for (unsigned int i=0; i< x.size(); ++i)
        result += x[i]*y[i];
    return result;
}

double FEM::CalculateNorm2(FemVec &u){
    return DotProduct(u,u); //square of norm of u
}

void FEM::CalculateBx(FemVec &result,  FemMat &B, FemVec &x){
    for (FemMat::iterator i = B.begin(); i != B.end(); ++i){
        double sum = 0.0;
        for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
            sum += j->second * x[j->first];  //multiply row of Mh with u
        result[i->first] = sum;
    }
}

void FEM::CalculateResidual(FemVec &result, FemVec &f, FemMat &B, FemVec &x){
    for (FemMat::iterator i = B.begin(); i != B.end(); ++i){
        double sum = 0.0;
        for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
            sum += j->second * x[j->first];  //multiply row of Mh with u
        result[i->first] = f[i->first] - sum;
    }
}

void FEM::UpdateVector(FemVec &result, double coefficient, FemVec &x){
    for (unsigned int i=0; i<x.size(); ++i){
        result[i] += coefficient * x[i];
    }
}

void FEM::UpdateSearchVector(FemVec &result, double coefficient, FemVec &x){
    for (unsigned int i=0; i<x.size(); ++i){
        result[i] = x[i] + coefficient * result[i];
    }
}

void FEM::ConjugateGradientSolve(FemVec &result, FemVec &f){
    FemVec Ad = result;
    FemVec r = result;
    CalculateResidual(r,f,A,result);   //first residual
    FemVec d = r;   //first search vector
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
    }
}

void FEM::GaussSeidelSolve(FemVec &result, FemVec &f){
    FemVec Ad = result;
    FemVec r = result;
    CalculateResidual(r,f,A,result);   //first residual
    FemVec d = r;   //first search vector
    double normr2 = CalculateNorm2(r);
    double epsilon2 = epsilon*epsilon;
    while (normr2>epsilon2){
        for (FemMat::iterator i = A.begin(); i != A.end(); ++i){
                double sum = 0.0;
                double aii = 1.0;
                for (map<int,double>::iterator j = i->second.begin(); j != i->second.end(); ++j)
                    if (j->first != i->first)
                        sum += j->second * result[j->first];  //multiply row of Mh with u
                    else
                        aii = j->second;
                result[i->first] = (f[i->first] - sum) / aii;
        }
        CalculateResidual(r,f,A,result);
        normr2 = CalculateNorm2(r);
    }
}

void FEM::Normalize(FemVec &u){
    double normu = sqrt(CalculateNorm2(u));
    for (unsigned int i=0; i<u.size(); ++i)
        u[i] /= normu;
}

double FEM::CalculateLambda(FemVec &u){
    CalculateBx(Au,A,u);
    double numerator = DotProduct(u,Au); //Calculate numerator
    double denominator = DotProduct(u,F); //Calculate denominator
    return numerator/denominator;
}

void FEM::SolveEigenProblem(){
    lambda = 0.0;
    double lambdaOld = 1.0;
    U = FemVec(A.size(),1.0); //set U to 1 to prevent trivial solution
    F = FemVec(A.size(),0.0);
    Au = FemVec(A.size(),0.0);

    int iter = 1;
    while (fabs((lambda-lambdaOld)/lambdaOld)>pow((double)10.00,-10)){
        lambdaOld = lambda;
        CalculateBx(F,M,U);
        ConjugateGradientSolve(U,F);   //solve linear system of equations to find U
        //GaussSeidelSolve(U,F);   //solve linear system of equations to find U
        Normalize(U);
        lambda = CalculateLambda(U);
        cout<<"lambda after iter "<<iter<<" = "<<setprecision(10)<<lambda<<endl;
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
        vector< double > corners(6, 0.0);
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
        my_local_matrix = my_element.integrate(grad(v_()) * grad(w_())-func<double>(K2)*v_()*w_());
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
        my_local_matrix = my_element.integrate(v_()*w_());
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
    cout<<"delta = "<<delta<<", epsilon = "<<epsilon<<", refLvl = "<<nLevel<<endl;
    cout<<endl;
    if (nLevel>0)
        cout<<"Splitting faces. Please be patient ..."<<endl<<endl<<"\"Breathe. Smile. and Go slowly.\""<<endl<<endl;
    SplitFaces();
    WriteK2("ksq.txt");
    SolveLocalMatrix();
    WriteMatrix(A,"A.txt");
    WriteMatrix(M,"M.txt");
    SolveEigenProblem();
    WriteU("eigenmode.txt");
    WriteLambda("lambda.txt");
    if (nLevel>0){
        cout<<endl<<"See! It's really simple. "<<endl;
        cout<<endl<<"\"Just sit quietly, doing nothing, spring comes and the grass grows by itself.\""<<endl;
        cout<<endl<<"Thank you for your patience. Have a wonderful day :^)"<<endl;
    }
    cout<<endl;
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
        FemVec point;
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

void FEM::WriteLambda(string filename){
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        outfile<<setprecision(10)<<lambda;
        outfile.close();
    }
}
