#pragma once

#include "Node.h"

class DirichletCondition
{
public:
    DirichletCondition(Node *node,
                       const int &direction,
                       const double &value);

    ~DirichletCondition();

    Node *getNode();

    int getDirection();

    double getValue();

    void applyDirichletCondition(const double &rate);

private:
    Node *node_;

    int direction_;

    double value_;
};

class NeumannCondition
{
public:
    NeumannCondition(Node *node,
                     const int &direction,
                     const double &value);

    ~NeumannCondition();

    Node *getNode();

    int getDirection();

    double getValue();

    void incrementValue(const double &value);

private:
    Node *node_;

    int direction_;

    double value_;
};