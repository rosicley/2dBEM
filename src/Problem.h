#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <unordered_set>
#include <time.h>
#include <iomanip>
#include "BoundaryConditions.h"
#include "BoundaryElement.h"
#include "mesh/Geometry.h"
#include <metis.h>
#include <petscksp.h>

#ifdef _WIN32
#include <direct.h>
#define getCurrentDir _getcwd
#define remove "del "
#else
#include <unistd.h>
#define getCurrentDir getcwd
#define remove "rm "
#endif

class Problem
{
public:
    Problem();

    Problem(const int &nIP, const double &collocParam, const double &offsetSourceP);

    ~Problem();

    void addNode(const int &index, const vecDouble &coordinate);

    void addBoundaryElement(const int &index, const std::string &name, const std::vector<Node *> &connection);

    void addMaterial(const double &young, const double &poisson = 0.0, const double &density = 0.0);
    
    void exportToParaviewGeometricMesh(const int &index);

    void exportToParaviewCollocationMesh_Potential(const int &index);

    void exportToParaviewCollocationMesh_Elasticity(const int &index);

    void exportToParaviewSourcePoints(const int &index);

    void createSourcePoints();

    void applyPotentialBoundaryConditions();

    void generateMesh(Geometry *geometry, const int &order, const std::string &algorithm = "AUTO", const bool &plotMesh = true, const bool &deleteFiles = false, const bool &showInfo = false);

    void gmshReading(const std::string &meshFile);

    int solvePotentialProblem();

    int solveElasticityProblem(const std::string &planeState);

    void applyElasticityBoundaryConditions();

    int computeInternalPoints();

    void addInternalPoints(std::vector<std::vector<double>> coord);

    void teste();

private:
    // Solid variables
    std::unordered_map<std::string, std::vector<BoundaryElement *>> lineElements_;
    std::vector<BoundaryElement *> elements_;
    std::vector<BoundaryElement *> partElements_;

    std::vector<Node *> nodes_;
    std::unordered_set<Node *> duplicatedNodes_;

    std::vector<CollocationPoint *> collocPoints_;
    std::vector<int> collocCondition_; //<collocation index, 0 for known flux or 1 for known potential>
    std::vector<SourcePoint *> sourcePoints_;

    std::vector<SourcePoint *> internalPoints_;
    std::vector<double> internalPotential_;
    std::vector<std::vector<double>> internalFlux_; //[potential, flux1, flux2]

    std::vector<Material *> materials_;

    // Analysis parameters
    IntegQuadrature *quadrature_;

    // Metis
    idx_t *elementPartition_;
    idx_t *nodePartition_;

    // PETSC variables
    Mat A, C;
    Vec b, x, All;
    PetscErrorCode ierr;
    PetscInt Istart, Iend, Idof, Ione, iterations, *dof;
    KSP ksp;
    PC pc;
    VecScatter ctx;
    // PetscScalar val, value;

    ///
    Geometry *geometry_;
    std::string current_working_dir_;
    int order_;

    Geometry *internalGeometry_;

    double collocParam_;
    bool sourceOut_;
    double offsetSourceOut_;
};