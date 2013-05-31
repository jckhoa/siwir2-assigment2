#include "FEM.h"
#include <sstream>

using namespace std;
int main(int argc, char** argv)
{
    double epsilon = 0.0000001;
    double delta = 0.01;
    double refLvl = 0;
    stringstream ss;
    if (argc >3){
        ss<<argv[1]<<" "<<argv[2]<<" "<<argv[3];
        ss>>delta>>epsilon>>refLvl;
    }
    else if (argc>2){
        ss<<argv[1]<<" "<<argv[2];
        ss>>delta>>epsilon;
    }
    else if (argc>1){
        ss<<argv[1];
        ss>>delta;
    }
    if (ss.fail()){
        cerr<<"One of the parameters is illegal."<<endl;
        return -1;
    }
    //FEM myfem(0.01,0.01);
    FEM myfem(epsilon, refLvl);
    myfem.SetDelta(delta);
    myfem.ReadData("unit_circle.txt");
    myfem.SolveAll();
    return 0;
}
