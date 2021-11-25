#pragma once

#include "Surface.h"
#include <math.h>
#include <unordered_map>

class Geometry
{
public:
    Geometry(const std::string &name);

    ~Geometry();

    Point *addPoint(std::vector<double> coordinates, const double &lcar = 1.0, const bool &discretization = true);

    Line *addLine(std::vector<Point *> points, const bool &discretization = true);

    Line *addCircle(std::vector<Point *> points, const bool &discretization = true); //{initial point, center point, end point}

    Line *addSpline(std::vector<Point *> points, const bool &discretization = true);

    LineLoop *addLineLoop(std::vector<Line *> lines, const bool &verify = false);

    // Surface *addSurface(std::vector<LineLoop *> lineLoop, const int &indexMaterial = 0, const double &thickness = 1.0);

    Surface *addSurface(std::vector<Line *> lineLoop, const int &indexMaterial = 0, const double &thickness = 1.0);

    SurfaceLoop *addSurfaceLoop(std::vector<Surface *> surfaces, const bool &verify = true);

    void addNeumannCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2 = {}, const std::vector<double> &directionX3 = {});

    void addNeumannCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2 = {}, const std::vector<double> &directionX3 = {});

    void addDirichletCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2 = {}, const std::vector<double> &directionX3 = {});

    void addDirichletCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2 = {}, const std::vector<double> &directionX3 = {});

    std::unordered_map<std::string, vecDouble> getDirichletCondition();

    std::unordered_map<std::string, vecDouble> getNeumannCondition();

    std::string getName();

    Point *getPoint(const std::string &name);

    Surface *getSurface(const std::string &name);

    std::unordered_map<std::string, Point *> getPoints();

    bool localGeometry();

    std::string createGmshCode();

    std::unordered_map<std::string, Surface *> getSurfaces();

    std::unordered_map<std::string, std::vector<std::string>> getRegions();

    bool verifyDuplicatedPoint(const std::string &namePoint);

    // void addCrack(std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary, const double &jradius, const double &lcarOfJintegral);

    // void addCrackOnGlobal(std::vector<Point *> points, const std::string &openBoundary, const double &jradius, const double &lcarOfJintegral, const double &lengthOffset, const double &lcarOfLocalBoundary);

    void transfiniteLine(std::vector<Line *> lines, const int &divisions, const double &progression = 1);

    // // void transfiniteSurface(std::vector<PlaneSurface *> surfaces, std::string oientation = "Left", std::vector<Point *> points = std::vector<Point *>());

    // double getRadiusJintegral();

    // void createGeometryFromCrack();

    // void addCrackPoint(const std::string &name, const vecDouble &coordinates);

    int getIndexMaterial(const std::string &name);

private:
    //Geometry
    std::unordered_map<std::string, Point *> points_;
    std::unordered_map<std::string, Line *> lines_;
    std::unordered_map<std::string, LineLoop *> lineLoops_;
    std::unordered_map<std::string, Surface *> surfaces_;
    std::unordered_map<std::string, SurfaceLoop *> surfaceLoops_;

    std::vector<std::vector<Line *>> transfiniteLines_;
    std::vector<int> divisionLines_;
    std::vector<int> progressionLines_;

    // std::unordered_map<std::string, Crack *> crackes_;

    ///BoundaryConditions
    std::unordered_map<std::string, vecDouble> neumannConditions_;
    std::unordered_map<std::string, vecDouble> diricheletConditions_;

    std::string name_;

    // std::string domain_; //local ou global
    // double Jradius_;
    // double lengthOffset_;
    // double lcarJ_;
    // double lcarOfLocalBoundary_;
    // int remeshNumber_;
};
