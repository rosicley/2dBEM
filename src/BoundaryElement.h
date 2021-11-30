#pragma once

#include "Node.h"
#include "IntegrationQuadrature.h"
#include "Material.h"

class BoundaryElement
{
public:
    BoundaryElement(const int &index, const std::vector<Node *> &geoConnection, const std::vector<CollocationPoint *> &collocConnection, IntegQuadrature *quadrature);

    BoundaryElement(const int &order);

    ~BoundaryElement();

    std::vector<Node *> getConnection();

    Node *getNode(const int &index);

    void crossSectionProperties(double &area, double &perim, double &Ix, double &Iy, double &Mex, double &Mey);

    void getShapeFunction(const double &xsi, vecDouble &phi);

    void getShapeFunctionAndDerivate(const double &xsi, vecDouble &phi, vecDouble &dphi_dxsi);

    int getIndex();

    std::vector<CollocationPoint *> getCollocationConnection();

    void interpolateGeometricalCoordinate(const double &xsi, vecDouble &coord);

    void coordOfSourcePoints(const double &xsi, const double &offset, vecDouble &coord);

    virtual void potentialContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH);

    virtual void calculateInternalFlux(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH);

    virtual void elasticityContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material);

    virtual void calculateInternalStress(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material);

    void elasticityBodyForceContribution(const std::vector<SourcePoint *> &sourcePoints, vecDouble &vecB, Material *material, const vecDouble &force);

protected:
    int index_;

    std::vector<Node *> geoConnection_;

    std::vector<CollocationPoint *> collocConnection_;

    IntegQuadrature *quadrature_;

    int order_;

    vecDouble xsis_;
};

class DiscontBoundaryElement : public BoundaryElement
{
public:
    DiscontBoundaryElement(const int &index, const std::vector<Node *> &geoConnection, const std::vector<CollocationPoint *> &collocConnection, IntegQuadrature *quadrature, const std::string &discont, const double &paramColloc);

    DiscontBoundaryElement(const int &order, const std::string &discont, const double &param);

    ~DiscontBoundaryElement();

    virtual void potentialContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH);

    void getCollocationShapeFunction(const double &xsi, vecDouble &phi);

    virtual void calculateInternalFlux(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH);

    virtual void elasticityContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material);

    virtual void calculateInternalStress(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material);

private:
    std::string discont_; //lef, right, both
    double collocParam_;
    vecDouble collocXsis_;
};