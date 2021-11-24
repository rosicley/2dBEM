// #pragma once

// #include "Line.h"
// #include "Surface.h"

// class CrackLine
// {
// public:
//   CrackLine();

//   CrackLine(const std::string &name, std::vector<Point *> points, Surface *surface, const std::string &openBoundary);

//   ~CrackLine();

//   std::string getName();

//   std::vector<Point *> getPoints();

//   Surface *getSurface();

//   std::string getGmshCodeCrackPlugin();

//   Point *getLastPoint();

//   Point *getFirstPoint();

//   void addLastPoint(Point *newPoint);

//     void addFirstPoint(Point *newPoint);

//   double getLastAngle();

//   std::pair<std::vector<Point *>, std::vector<Point *>> getPointsOfLocalGeometry();

//   void addLastPointsOfLocalGeometry(Point *boundary1, Point *boundary2);

//   void setSurface(Surface *newSurface);

//   void setAuxPoint(Point *auxPoint);

//   Point *getAuxPoint();

//   std::string getOpenBoundary();

// private:
//   std::string name_;
//   std::string openBoundary_;
//   std::vector<Point *> points_;

//   Surface *surface_;
//   std::vector<Point *> geometryBoundary1_; //only local
//   std::vector<Point *> geometryBoundary2_; //only local
//   Point *auxPoint_;                        //transition point between geometryBoundary1 and geometryBoundary2

// };