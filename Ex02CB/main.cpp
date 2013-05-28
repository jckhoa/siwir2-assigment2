

#include "FEM.h"

int main(int argc, char** argv)
{
    FEM myfem(0.01,0.01);
    myfem.ReadData("unit_circle.txt");
    myfem.WriteK2("ksq.txt");
    return 0;
}
