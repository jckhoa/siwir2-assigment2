#include "FEM.h"
#include <math.h>
int main(int argc, char** argv)
{
    //FEM myfem(0.01,0.01);
    FEM myfem(pow((double)10.00,-10.00));
    myfem.SetDelta(0.01);
    myfem.ReadData("unit_circle.txt");
    myfem.WriteK2("ksq.txt");
    myfem.SetUpNeighborLists();
    myfem.Solve();
    myfem.WriteAh("A.txt");
    myfem.WriteMh("M.txt");
    return 0;
}
