#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <unordered_set>
#include <time.h>
#include <iomanip>
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

    int computeInternalPointsPotentialProblem();

    int computeInternalPointsElasticityProblem();

    void addInternalPoints(std::vector<std::vector<double>> coord);

    void addInternalPoints(Surface *surface, const std::vector<std::vector<double>> &coord);

    void teste();

    void teste2();

    void coupleLines(std::vector<std::vector<Line *>> coupledLines);

    void addBodyForceToSurface(Surface *surface, const vecDouble &values);

private:
    // Elements
    std::unordered_map<std::string, std::vector<BoundaryElement *>> lineElements_;
    std::unordered_map<std::string, std::vector<BoundaryElement *>> subElements_;
    std::vector<BoundaryElement *> elements_;

    // Geometric nodes
    std::vector<Node *> nodes_;
    std::unordered_set<Node *> discontinuousNodes_;

    // Collocation points
    std::vector<CollocationPoint *> collocPoints_;
    std::vector<int> collocCondition_; //<collocation index, 0 for known neumman condition or 1 for known dirichelet condition, 2 or 3 for coupled lines>
    
    // Source points
    std::unordered_map<std::string, std::vector<SourcePoint *>> subSourcePoints_;
    std::vector<SourcePoint *> sourcePoints_;

    // Coupling lines
    std::unordered_map<int, int> coupledCollocFirst_;
    std::unordered_map<int, int> coupledCollocSecond_;

    // Internal Points
    std::vector<SourcePoint *> internalPoints_;
    std::vector<double> internalPotential_;
    std::vector<std::vector<double>> internalFlux_;
    std::unordered_map<std::string, std::vector<SourcePoint *>> subInternalPoints_;
    std::vector<std::vector<double>> internalDisplacements_;
    std::vector<std::vector<double>> internalStresses_;

    // Materials
    std::vector<Material *> materials_;

    // Body force
    std::unordered_map<std::string, std::vector<double>> bodyForces_;

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

    // Geometry
    Geometry *geometry_;
    std::string current_working_dir_;
    int order_;

    // Parameters
    double collocParam_;
    bool sourceOut_;
    double offsetSourceOut_;
};