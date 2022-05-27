static char help[] = "Code to solve potential and elasticity problem using boundary element method";

#include "src/Problem.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char *)0, help);

    // #include "examples/exPaper.h"
    // #include "examples/chapa2.h"
    // #include "examples/viga.h"
    // #include "examples/viga_subregioes.h"
    // #include "examples/teste.h"
    

#include "examples/chapa.h"

    PetscFinalize();
 // teste
    return 0;
}