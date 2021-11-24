#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "../Node.h"

class Point
{
public:
  Point(const std::string &name,
        const vecDouble &coordinates,
        const double &lcar,
        const bool &discretization);

  ~Point();

  std::string getName();

  std::string getGmshCode();

  void addNodeToPoint(Node *node);

  void getCoordinates(vecDouble &coord);

  void setCoordinates(const vecDouble &newCoordinates);

  double getlcar();

  Node *getPointNode();

  void setlcar(const double &newlcar);

private:
  std::string name_;      // Gmsh Physical entity name
  vecDouble coordinates_; // Coordinates vector (x,y)
  double lcar_;           // Characteristic length Gmsh parameter
  bool discretization_;   // Choose to discretize a point with a mesh node
  Node *pointNode_;       // Defines the mesh node discretizing the point
};