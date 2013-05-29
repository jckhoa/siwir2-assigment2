
#include "FEM.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include "Colsamm/Colsamm.h"

using namespace std;
using namespace ::_COLSAMM_;

double FEM::delta = 0;

void FEM::SolveEigenProblem(){
    double lambda = 0;
    double lambdaOld = 1;
    vector<double> uh(vertices.size(),0.0);
    vector<double> f(vertices.size(),0.0);
    while (fabs((lambda-lambdaOld)/lambdaOld)>epsilon){
        lambdaOld = lambda;
//        f = MatrixMultiplication(Mh,uh);
 //       uh = LinearSolve(Ah,f);
  //      Normalize(uh);
   //     lambda = CalculateLambda(uh)
    }
}

double FEM::K2(double x,double y){
    return (100+delta)* exp(-50*(x*x+y*y))-100;
}

void FEM::SetUpNeighborLists(){
    for (unsigned int i=0; i< faces.size(); ++i){
        int vertex0 = faces[i][0];
        int vertex1 = faces[i][1];
        int vertex2 = faces[i][2];
        AM[vertex0][vertex0] = GlobalMatrices(0.0,0.0);
        AM[vertex0][vertex1] = GlobalMatrices(0.0,0.0);
        AM[vertex0][vertex2] = GlobalMatrices(0.0,0.0);
        AM[vertex1][vertex0] = GlobalMatrices(0.0,0.0);
        AM[vertex1][vertex1] = GlobalMatrices(0.0,0.0);
        AM[vertex1][vertex2] = GlobalMatrices(0.0,0.0);
        AM[vertex2][vertex0] = GlobalMatrices(0.0,0.0);
        AM[vertex2][vertex1] = GlobalMatrices(0.0,0.0);
        AM[vertex2][vertex2] = GlobalMatrices(0.0,0.0);

    }
}
double myfunc(double x, double y){ return 0.0;}
void FEM::Solve(){
   //  function<double<const FEM&, double,double)> GetK2 = &FEM::K2;
    for (unsigned int i=0; i< faces.size(); ++i){
        //find local matrix for each face
        ELEMENTS::Triangle my_elementA, my_elementM;
        vector< vector< double > > my_local_matrixA, my_local_matrixM;
        vector<double> corners(6, 0.0);
        // array corners contains the x- and y-coordinates of the
        // triangle corners in the order x0, y0, x1, y1, x2, y2
        int vertex0 = faces[i][0];
        int vertex1 = faces[i][1];
        int vertex2 = faces[i][2];

        corners[0] = vertices[vertex0].x;
        corners[1] = vertices[vertex0].y;
        corners[2] = vertices[vertex1].x;
        corners[3] = vertices[vertex1].y;
        corners[4] = vertices[vertex2].x;
        corners[5] = vertices[vertex2].y;
        // pass the corners to the finite element
        my_elementA(corners);
        my_local_matrixA =
        my_elementA.integrate(grad(v_()) * grad(w_())-func<double>(K2)*v_()*w_());

        //transfering value from local matrix to global matrix Ah
        AM[vertex0][vertex0].A += my_local_matrixA[0][0];
        AM[vertex0][vertex1].A += my_local_matrixA[0][1];
        AM[vertex0][vertex2].A += my_local_matrixA[0][2];
        AM[vertex1][vertex0].A += my_local_matrixA[1][0];
        AM[vertex1][vertex1].A += my_local_matrixA[1][1];
        AM[vertex1][vertex2].A += my_local_matrixA[1][2];
        AM[vertex2][vertex0].A += my_local_matrixA[2][0];
        AM[vertex2][vertex1].A += my_local_matrixA[2][1];
        AM[vertex2][vertex2].A += my_local_matrixA[2][2];



        my_elementM(corners);
        my_local_matrixM =
        my_elementM.integrate(v_()*w_());

        //transfering value from local matrix to global matrix Ah
        AM[vertex0][vertex0].M += my_local_matrixM[0][0];
        AM[vertex0][vertex1].M += my_local_matrixM[0][1];
        AM[vertex0][vertex2].M += my_local_matrixM[0][2];
        AM[vertex1][vertex0].M += my_local_matrixM[1][0];
        AM[vertex1][vertex1].M += my_local_matrixM[1][1];
        AM[vertex1][vertex2].M += my_local_matrixM[1][2];
        AM[vertex2][vertex0].M += my_local_matrixM[2][0];
        AM[vertex2][vertex1].M += my_local_matrixM[2][1];
        AM[vertex2][vertex2].M += my_local_matrixM[2][2];
    }
}

void FEM::WriteAh(string filename){
    /* This output is for view in Matlab
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (int i=0;i<vertices.size();++i){
            for (int j=0;j<vertices.size();++j){
                map<int,GlobalMatrices>::iterator it = AM[i].find(j);
                if (it == AM[i].end())
                    outfile<<0.000000<<" ";
                else
                    outfile<<it->second.A<<" ";
            }
            outfile<<endl;
        }

        outfile.close();
    }
    */

    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (int i=0;i<vertices.size();++i){
            for (int j=0;j<vertices.size();++j){
                map<int,GlobalMatrices>::iterator it = AM[i].find(j);
                if (it != AM[i].end())
                    outfile<<i<<" "<<j<<" "<<it->second.A<<endl;
            }
        }

        outfile.close();
    }
}

void FEM::WriteMh(string filename){
    /* This output is for view in Matlab
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (int i=0;i<vertices.size();++i){
            for (int j=0;j<vertices.size();++j){
                map<int,GlobalMatrices>::iterator it = AM[i].find(j);
                if (it == AM[i].end())
                    outfile<<0.000000<<" ";
                else
                    outfile<<it->second.M<<" ";
            }
            outfile<<endl;
        }

        outfile.close();
    }
    */
    ofstream outfile(filename.c_str());
    if (outfile.good()){
        for (int i=0;i<vertices.size();++i){
            for (int j=0;j<vertices.size();++j){
                map<int,GlobalMatrices>::iterator it = AM[i].find(j);
                if (it != AM[i].end())
                    outfile<<i<<" "<<j<<" "<<it->second.M<<endl;
            }
        }
        outfile.close();
    }
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
        Coordinates point(x,y);
        vertices.insert(pair<int, Coordinates>(index,point));
    }
    //check if all vertices were read correctly
    /*
    for (map<int,Coordinates>::iterator it= vertices.begin(); it !=vertices.end();++it)
        cout<<it->second.x<<" "<<it->second.y<<" "<<endl;
    */

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
        for (map<int, Coordinates>::iterator it=vertices.begin(); it!=vertices.end(); ++it){
            double val = K2(it->second.x, it->second.y);
            outfile<<it->second.x<<" "<<it->second.y<<" "<<val<<endl;
        }
        outfile.close();
    }
}

void SetUpNeighborLists();
