#pragma once
#include <math.h>
#include "General.h"

//MOTHER CLASS
class IntegQuadrature
{
public:
    ///Constructor - Defines a new IntegQuadrature
    IntegQuadrature();

    /// Destructor - Delete the IntegQuadrature
    virtual ~IntegQuadrature();

    /// Returns the number of one-dimensional integration points
    /// @return integration points number (1D)
    int get1DPointsNumber();

    /// Returns the number of two-dimensional integration points
    /// @return integration points number (2D)
    int get2DPointsNumber();

    /// Modifies the arguments with the integration points
    /// @param int point number
    /// @param double dimensionless coordinate
    /// @param double weight
    void get1DIntegrationPoint(const int &point, double &xsi1, double &weight);

    /// Modifies the arguments with the integration points
    /// @param int point number
    /// @param double dimensionless coordinate (direction 1)
    /// @param double dimensionless coordinate (direction 2)
    /// @param double weight
    void get2DIntegrationPoint(const int &point, double &xsi1, double &xsi2, double &weight);

    /// Modifies the arguments with the integration points
    /// @param int point number
    /// @param double dimensionless coordinate (direction 1)
    /// @param double dimensionless coordinate (direction 2)
    /// @param double dimensionless coordinate (direction 3)
    /// @param double weight
    void get3DIntegrationPoint(const int &point, double &xsi1, double &xsi2, double &xsi3, double &weight);

protected:
    /// Number of points in one-dimensional integration IntegQuadrature
    int nGauss1D_;

    /// Number of points in two-dimensional integration IntegQuadrature
    int nPoints2D_;

    /// Number of points in three-dimensional integration IntegQuadrature
    int nPoints3D_;

    /// One-dimensional IntegQuadrature ([0][i]: coordinate, [1][i]: weight)
    matDouble gaussPoints1D_;

    /// Two-dimensional IntegQuadrature ([0][i]: coordinate 1, [1][i]: coordinate 2, [2][i]: weight)
    matDouble integPoints2D_;

    /// Three-dimensional IntegQuadrature ([0][i]: coordinate 1, [1][i]: coordinate 2, [2][i]: coordinate 3, [3][i]: weight)
    matDouble integPoints3D_;
};

//GAUSS1D CLASS
class Gauss1D : public IntegQuadrature
{
public:
    /// Constructor - Defines a one-dimensional integration IntegQuadrature
    Gauss1D(const int &nPoints1D);

    /// Destructor - Delete the Gauss1D
    virtual ~Gauss1D();

private:
    /// Sets the one-dimensional integration IntegQuadrature to IntegQuadrature object
    void setIntegrationPoints();
};

//GAUSS2D CLASS
class Gauss2D : public Gauss1D
{
public:
    /// Constructor - Defines a unidimensional and bidimensional integration IntegQuadrature
    /// @param int points number of both integration IntegQuadratures
    Gauss2D(const int &nPoints);

    /// Constructor - Defines a unidimensional and bidimensional integration IntegQuadrature
    /// @param int points number of unidimensional integration IntegQuadrature
    /// @param int points number of bidimensional integration IntegQuadrature
    Gauss2D(const int &nPoints1D, const int &nPoints2D);

    /// Destructor - Delete Gauss2D
    virtual ~Gauss2D();

private:
    /// Sets the bidimensional integration IntegQuadrature to IntegQuadrature object. This function uses the one-dimensional square already calculated.
    void setIntegrationPoints();

    /// Sets the bidimensional integration IntegQuadrature to IntegQuadrature object. This function is used when nPoints1D_ * nPoints1D_ != nPoints2D_
    void setIntegrationPointsDIF();
};

// //GAUSS3D CLASS
// class Gauss3D : public Gauss2D
// {
// public:
//     /// Constructor - Defines a unidimensional, bidimensional and three-dimensional integration IntegQuadrature
//     /// @param int points number of integration IntegQuadratures
//     Gauss3D(const int &nPoints);

//     virtual ~Gauss3D();

//     void setIntegrationPoints();
// };

//HAMMER2D
class Hammer2D : public Gauss1D
{
public:
    /// Constructor - Defines a unidimensional and bidimensional integration IntegQuadrature
    /// @param int points number of unidimensional integration IntegQuadrature (GAUSS)
    /// @param int points number of bidimensional integration IntegQuadrature (HAMMER)
    Hammer2D(const int &nGauss, const int &nHammer);

    /// Destructor - Delete the Hammer2D
    virtual ~Hammer2D();

private:
    /// Sets the bidimensional integration IntegQuadrature of Hammer to IntegQuadrature object.
    void setIntegrationPoints();
};