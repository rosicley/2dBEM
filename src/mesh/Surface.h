#pragma once

#include "Line.h"
#include "../Material.h"

class Surface
{
public:
    Surface();

    Surface(const std::string &name, std::vector<LineLoop *> lineLoop, const int &indexMaterial = 0, const double &thickness = 1.0);

    ~Surface();

    std::string getName();

    double getThickness();

    int getIndexMaterial();

    std::string getGmshCode();

    LineLoop *getLineLoop(const int &index);

    void getIndexMaterialAndThickness(int &indexMaterial, double &thickness);

    std::vector<LineLoop *> getLineLoops();

private:
    std::string name_;
    double thickness_;
    std::vector<LineLoop *> lineLoop_;
    int indexMaterial_;
};

class SurfaceLoop
{
public:
    SurfaceLoop();

    SurfaceLoop(const std::string &name, std::vector<Surface *> surfaces);

    ~SurfaceLoop();

    std::string getName();

    Surface *getSurface(const int &index);

    std::vector<Surface *> getSurfaces();

    std::string getGmshCode();

    void verification();

private:
    std::string name_;
    std::vector<Surface *> surfaces_;
};
