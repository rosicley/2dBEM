#pragma once
#include <iostream>
#include <string>

class Material
{
public:
    /// Constructor - Defines a Material
    /// @param double young @param double poisson @param double density
    Material(const double &young, const double &poisson, const double &density);

    /// Destructor - delete the material
    ~Material();

    /// Returns the young's modulus
    /// @return young's modulus
    double getYoung();

    /// Returns the poisson's ratio
    /// @return poisson's ratio
    double getPoisson();

    /// Returns the density
    /// @return density
    double getDensity();

    /// Modifies the arguments with the material properties
    /// @param double young @param double poison @param double density
    void getProperties(double &young, double &poisson, double &density);

private:
    /// Young's modulus
    double young_;

    /// Poisson's ratio
    double poisson_;

    /// Density
    double density_;
};