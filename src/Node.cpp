#include "Node.h"

// -------------- //
// GEOMETRIC NODE //
// -------------- //

Node::Node(const int &index, const vecDouble &initialCoordinate)
{
    index_ = index;
    coordinate_.resize(3, 0.0);
    coordinate_[0] = initialCoordinate[0];
    coordinate_[1] = initialCoordinate[1];
    coordinate_[2] = initialCoordinate[2];
}

Node::Node()
{
}

Node::~Node()
{
}

int Node::getIndex()
{
    return index_;
}

void Node::getCoordinate(vecDouble &coord)
{
    coord[0] = coordinate_[0];
    coord[1] = coordinate_[1];
    coord[2] = coordinate_[2];
}

vecDouble Node::getCoordinate()
{
    return coordinate_;
}

double Node::getCoordinate(const int &direction)
{
    return coordinate_[direction];
}

void Node::setCoordinate(const vecDouble &newCoord)
{
    coordinate_[0] = newCoord[0];
    coordinate_[1] = newCoord[1];
    coordinate_[2] = newCoord[2];
}

// ---------------- //
// COLOCATION POINT //
// ---------------- //

CollocationPoint::CollocationPoint(Node *node) : Node(node->getIndex(), node->getCoordinate())
{
    flux_ = 0.0;
    potential_ = 0.0;
    force_.resize(3, 0.0);
    displacement_.resize(3, 0.0);
}

CollocationPoint::~CollocationPoint()
{
}

double CollocationPoint::getPotential()
{
    return potential_;
}

double CollocationPoint::getFlux()
{
    return flux_;
}

double CollocationPoint::getForce(const int &direction)
{
    return force_[direction];
}

double CollocationPoint::getDisplacement(const int &direction)
{
    return displacement_[direction];
}

double CollocationPoint::getCurrentCoordinate(const int &direction)
{
    return (coordinate_[direction] + displacement_[direction]);
}

void CollocationPoint::setPotential(const double &potential)
{
    potential_ = potential;
}

void CollocationPoint::setFlux(const double &flux)
{
    flux_ = flux;
}

void CollocationPoint::setForce(const vecDouble &force)
{
    force_ = force;
}

void CollocationPoint::setDisplacement(const vecDouble &displacement)
{
    displacement_ = displacement;
}

void CollocationPoint::setForce(const double &value, const int &direction)
{
    force_[direction] = value;
}

void CollocationPoint::setDisplacement(const double &value, const int &direction)
{
    displacement_[direction] = value;
}

// ------------ //
// SOURCE POINT //
// ------------ //

SourcePoint::SourcePoint(CollocationPoint *coloc) : Node(coloc->getIndex(), coloc->getCoordinate())
{
}

SourcePoint::SourcePoint(const int &index, const vecDouble &initialCoordinate) : Node(index, initialCoordinate)
{
}

SourcePoint::~SourcePoint()
{
}

bool SourcePoint::checkElement(const int &indexElement)
{
    bool check;
    if (insideElements_.count(indexElement) == 0)
    {
        check = false;
    }
    else
    {
        check = true;
    }
    return check;
}

void SourcePoint::addInsideElement(const int &indexElement)
{
    insideElements_.insert(indexElement);
}
