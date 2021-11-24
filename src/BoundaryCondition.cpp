#include "BoundaryConditions.h"

/// DIRICHLET
DirichletCondition::DirichletCondition(Node *node,
                                       const int &direction,
                                       const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

DirichletCondition::~DirichletCondition()
{
}

Node *DirichletCondition::getNode()
{
    return node_;
}

int DirichletCondition::getDirection()
{
    return direction_;
}

double DirichletCondition::getValue()
{
    return value_;
}

///NEUMANN
NeumannCondition::NeumannCondition(Node *node,
                                   const int &direction,
                                   const double &value)
{
    node_ = node;
    direction_ = direction;
    value_ = value;
}

NeumannCondition::~NeumannCondition() {}

Node *NeumannCondition::getNode()
{
    return node_;
}

int NeumannCondition::getDirection()
{
    return direction_;
}

double NeumannCondition::getValue()
{
    return value_;
}

void NeumannCondition::incrementValue(const double &value)
{
    value_ += value;
}