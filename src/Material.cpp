#include "Material.h"

Material::Material(const double &young, const double &poisson, const double &density)
{
    young_ = young;
    poisson_ = poisson;
    density_ = density;
}

Material::~Material() {}

double Material::getYoung()
{
    return young_;
}

double Material::getPoisson()
{
    return poisson_;
}

double Material::getDensity()
{
    return density_;
}

void Material::getProperties(double &young, double &poisson, double &density)
{
    young = young_;
    poisson = poisson_;
    density = density_;
}