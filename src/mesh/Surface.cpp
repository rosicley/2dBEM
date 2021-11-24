#include "Surface.h"

/////////////
// Surface //
/////////////

Surface::Surface() {}

Surface::Surface(const std::string &name, std::vector<LineLoop *> lineLoop, const int &indexMaterial, const double &thickness)
{
    name_ = name;
    indexMaterial_ = indexMaterial;
    thickness_ = thickness;
    for (LineLoop *ll : lineLoop)
    {
        lineLoop_.push_back(ll);
    }
}

Surface::~Surface() {}

std::string Surface::getName()
{
    return name_;
}

double Surface::getThickness()
{
    return thickness_;
}

int Surface::getIndexMaterial()
{
    return indexMaterial_;
}

std::string Surface::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = news; Plane Surface(" << name_ << ") = {";
    for (int i = 0, nLineLoop = lineLoop_.size(); i < nLineLoop; i++)
    {
        text << lineLoop_[i]->getName();
        if (nLineLoop - i != 1)
        {
            text << ", ";
        }
    }
    text << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
    return text.str();
}

LineLoop *Surface::getLineLoop(const int &index)
{
    return lineLoop_[index];
}

void Surface::getIndexMaterialAndThickness(int &indexMaterial, double &thickness)
{
    indexMaterial = indexMaterial_;
    thickness = thickness_;
}

/////////////////
// SurfaceLoop //
/////////////////

SurfaceLoop::SurfaceLoop() {}

SurfaceLoop::SurfaceLoop(const std::string &name, std::vector<Surface *> surfaces)
{
    name_ = name;
    int nSurface = surfaces.size();
    surfaces_.resize(nSurface);
    for (int i = 0; i < nSurface; i++)
    {
        surfaces_[i] = surfaces[i];
    }
}

SurfaceLoop::~SurfaceLoop() {}

std::string SurfaceLoop::getName()
{
    return name_;
}

Surface *SurfaceLoop::getSurface(const int &index)
{
    return surfaces_[index];
}

std::vector<Surface *> SurfaceLoop::getSurfaces()
{
    return surfaces_;
}

std::string SurfaceLoop::getGmshCode()
{
    std::stringstream text;
    text << name_ << " = newsl; Surface Loop(" << name_ << ") = {";
    for (size_t i = 0; i < surfaces_.size(); i++)
    {
        text << surfaces_[i]->getName();
        if (i != (surfaces_.size() - 1))
            text << ", ";
    }
    text << "};\n//\n";
    return text.str();
}

void SurfaceLoop::verification()
{
    if (surfaces_.size() < 4)
    {
        std::stringstream text;
        text << "The surfaces ";
        for (size_t j = 0; j < surfaces_.size(); j++)
        {
            text << surfaces_[j]->getName();
            if (j != (surfaces_.size() - 1))
                text << ", ";
        }
        text << " do not form a closed volume." << std::endl;
        std::cout << text.str();
    }
}