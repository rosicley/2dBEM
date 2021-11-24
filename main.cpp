static char help[] = "Code to solve potential problem using boundary element method";

#include "src/Problem.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);

// #include "examples/exPaper.h"
#include "examples/barra.h"

    PetscFinalize();

    return 0;
}