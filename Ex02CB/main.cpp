#include "FEM.h"
#include <math.h>
int main(int argc, char** argv)
{
    //FEM myfem(0.01,0.01);
    FEM myfem(0.00001);
    myfem.SetDelta(0.1);
    myfem.ReadData("unit_circle.txt");
    myfem.SolveAll();
    return 0;
}
