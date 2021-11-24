#pragma once
#include "General.h"
#include <unordered_set>

// -------------- //
// GEOMETRIC NODE //
// -------------- //

class Node
{
public:
    /// Constructor - Defines a new node
    /// @param int index @param vecDouble coordinates @param int indexFE (optional)
    Node(const int &index, const vecDouble &initialCoordinate);

    Node();

    /// Destructor - Delete the node
    virtual ~Node();

    /// Returns the node or control point index
    /// @return node or control point index
    int getIndex();

    /// Modifies a vector with the initial coordinates
    /// @param vecDouble vector that will be modified
    void getCoordinate(vecDouble &coord);

    vecDouble getCoordinate();

    double getCoordinate(const int &direction);

    void setCoordinate(const vecDouble &newCoord);

protected:
    /// Index
    int index_;

    /// Initial nodal coordinate vector
    vecDouble coordinate_;
};

// ---------------- //
// COLOCATION POINT //
// ---------------- //

class CollocationPoint : public Node
{
public:
    CollocationPoint(Node *node);

    ~CollocationPoint();

    double getPotential();

    double getFlux();

    double getForce(const int &direction);

    double getDisplacement(const int &direction);

    double getCurrentCoordinate(const int &direction);

    void setPotential(const double &potential);

    void setFlux(const double &flux);

    void setForce(const vecDouble &force);

    void setDisplacement(const vecDouble &displacement);

    void setForce(const double &value, const int &direction);

    void setDisplacement(const double &value, const int &direction);

private:
    double potential_;
    double flux_;
    vecDouble force_;
    vecDouble displacement_;
};

// ------------ //
// SOURCE POINT //
// ------------ //

class SourcePoint : public Node
{
public:
    SourcePoint(CollocationPoint *coloc);

    SourcePoint(const int &index, const vecDouble &initialCoordinate);

    ~SourcePoint();

    bool checkElement(const int &indexElement);

    void addInsideElement(const int &indexElement);

private:
    std::unordered_set<int> insideElements_; //índices dos elementos que possuem esse ponto fonte
};