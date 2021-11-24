#include "IntegrationQuadrature.h"

//MOTHER CLASS
IntegQuadrature::IntegQuadrature()
{
}

IntegQuadrature::~IntegQuadrature()
{
}

int IntegQuadrature::get1DPointsNumber()
{
    return nGauss1D_;
}

int IntegQuadrature::get2DPointsNumber()
{
    return nPoints2D_;
}

void IntegQuadrature::get1DIntegrationPoint(const int &point, double &xsi1, double &weight)
{
    xsi1 = gaussPoints1D_[0][point];
    weight = gaussPoints1D_[1][point];
}

void IntegQuadrature::get2DIntegrationPoint(const int &point, double &xsi1, double &xsi2, double &weight)
{
    xsi1 = integPoints2D_[0][point];
    xsi2 = integPoints2D_[1][point];
    weight = integPoints2D_[2][point];
}

void IntegQuadrature::get3DIntegrationPoint(const int &point, double &xsi1, double &xsi2, double &xsi3, double &weight)
{
    xsi1 = integPoints3D_[0][point];
    xsi2 = integPoints3D_[1][point];
    xsi3 = integPoints3D_[2][point];
    weight = integPoints3D_[3][point];
}

//GAUSS 1D
Gauss1D::Gauss1D(const int &nPoints) : IntegQuadrature()
{
    nGauss1D_ = nPoints;
    setIntegrationPoints();
}

Gauss1D::~Gauss1D()
{
}

void Gauss1D::setIntegrationPoints()
{
    setColumnToMatDoubleWith2Lines(gaussPoints1D_, nGauss1D_);

    const double pi = std::acos(-1.0);
    double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
    int mga = (nGauss1D_ + 1) / 2;
    xmga = 0.0;
    xlga = 1.0;

    for (int iga = 1; iga <= mga; iga++)
    {
        zga = std::cos(pi * (double(iga) - 0.25) / (double(nGauss1D_) + 0.5));
    g1:
        p1ga = 1.0;
        p2ga = 0.0;
        for (int jga = 1; jga <= nGauss1D_; jga++)
        {
            p3ga = p2ga;
            p2ga = p1ga;
            p1ga = ((2.0 * double(jga) - 1.0) * zga * p2ga - (double(jga) - 1.0) * p3ga) / (double(jga));
        }

        ppga = nGauss1D_ * (zga * p1ga - p2ga) / (zga * zga - 1.0);
        z1ga = zga;
        zga = z1ga - p1ga / ppga;

        if (fabs(zga - z1ga) > 1.0e-15)
            goto g1;

        gaussPoints1D_[0][iga - 1] = xmga - xlga * zga;
        gaussPoints1D_[0][nGauss1D_ - iga] = xmga + xlga * zga;
        gaussPoints1D_[1][iga - 1] = 2.0 * xlga / ((1.0 - zga * zga) * ppga * ppga);
        gaussPoints1D_[1][nGauss1D_ - iga] = gaussPoints1D_[1][iga - 1];
    }
}

//GAUSS 2D
Gauss2D::Gauss2D(const int &nPoints) : Gauss1D(nPoints)
{
    nPoints2D_ = nPoints * nPoints;
    setIntegrationPoints();
}

Gauss2D::Gauss2D(const int &nPoints1D, const int &nPoints2D) : Gauss1D(nPoints1D)
{
    nPoints2D_ = nPoints2D * nPoints2D;
    setIntegrationPointsDIF();
}

Gauss2D::~Gauss2D()
{
}

void Gauss2D::setIntegrationPoints()
{
    setColumnToMatDoubleWith3Lines(integPoints2D_, nPoints2D_);

    int aux = 0;

    for (int i = 0; i < nGauss1D_; i++)
    {
        for (int j = 0; j < nGauss1D_; j++)
        {
            integPoints2D_[0][aux + j] = gaussPoints1D_[0][i];

            integPoints2D_[1][aux + j] = gaussPoints1D_[0][j];

            integPoints2D_[2][aux + j] = gaussPoints1D_[1][j] * gaussPoints1D_[1][i];
        }
        aux += nGauss1D_;
    }
}

void Gauss2D::setIntegrationPointsDIF()
{
    int nPoints = sqrt(nPoints2D_);
    matDouble gauss1D;
    setColumnToMatDoubleWith2Lines(gauss1D, nPoints);

    const double pi = std::acos(-1.0);
    double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
    int mga = (nPoints + 1) / 2;
    xmga = 0.0;
    xlga = 1.0;

    for (int iga = 1; iga <= mga; iga++)
    {
        zga = std::cos(pi * (double(iga) - 0.25) / (double(nPoints) + 0.5));
    g1:
        p1ga = 1.0;
        p2ga = 0.0;
        for (int jga = 1; jga <= nPoints; jga++)
        {
            p3ga = p2ga;
            p2ga = p1ga;
            p1ga = ((2.0 * double(jga) - 1.0) * zga * p2ga - (double(jga) - 1.0) * p3ga) / (double(jga));
        }

        ppga = nPoints * (zga * p1ga - p2ga) / (zga * zga - 1.0);
        z1ga = zga;
        zga = z1ga - p1ga / ppga;

        if (fabs(zga - z1ga) > 1.0e-15)
            goto g1;

        gauss1D[0][iga - 1] = xmga - xlga * zga;
        gauss1D[0][nPoints - iga] = xmga + xlga * zga;
        gauss1D[1][iga - 1] = 2.0 * xlga / ((1.0 - zga * zga) * ppga * ppga);
        gauss1D[1][nPoints - iga] = gauss1D[1][iga - 1];
    }

    setColumnToMatDoubleWith3Lines(integPoints2D_, nPoints2D_);

    int aux = 0;

    for (int i = 0; i < nPoints; i++)
    {
        for (int j = 0; j < nPoints; j++)
        {
            integPoints2D_[0][aux + j] = gauss1D[0][i];

            integPoints2D_[1][aux + j] = gauss1D[0][j];

            integPoints2D_[2][aux + j] = gauss1D[1][j] * gauss1D[1][i];
        }
        aux += nPoints;
    }
}

// //GAUSS 3D
// Gauss3D::Gauss3D(const int &nPoints)
// {
// }

// Gauss3D::~Gauss3D();

// void Gauss3D::setIntegrationPoints();

//HAMMER 2D
Hammer2D::Hammer2D(const int &nGauss, const int &nHammer) : Gauss1D(nGauss)
{
    nPoints2D_ = nHammer;
    setIntegrationPoints();
}

Hammer2D::~Hammer2D()
{
}

void Hammer2D::setIntegrationPoints()
{
    setColumnToMatDoubleWith3Lines(integPoints2D_, nPoints2D_);

    if (nPoints2D_ == 7)
    {
        integPoints2D_[0][0] = 0.333333333333333;
        integPoints2D_[0][1] = 0.797426985353087;
        integPoints2D_[0][2] = 0.101286507323456;
        integPoints2D_[0][3] = 0.101286507323456;
        integPoints2D_[0][4] = 0.470142064105115;
        integPoints2D_[0][5] = 0.059715871789770;
        integPoints2D_[0][6] = 0.470142064105115;

        integPoints2D_[1][0] = 0.333333333333333;
        integPoints2D_[1][1] = 0.101286507323456;
        integPoints2D_[1][2] = 0.797426985353087;
        integPoints2D_[1][3] = 0.101286507323456;
        integPoints2D_[1][4] = 0.470142064105115;
        integPoints2D_[1][5] = 0.470142064105115;
        integPoints2D_[1][6] = 0.059715871789770;

        integPoints2D_[2][0] = 0.112500000000000;
        integPoints2D_[2][1] = 0.062969590272414;
        integPoints2D_[2][2] = 0.062969590272414;
        integPoints2D_[2][3] = 0.062969590272414;
        integPoints2D_[2][4] = 0.066197076394253;
        integPoints2D_[2][5] = 0.066197076394253;
        integPoints2D_[2][6] = 0.066197076394253;
    }
    else if (nPoints2D_ == 12)
    {
        integPoints2D_[0][0] = 0.501426509658179;
        integPoints2D_[0][1] = 0.249286745170910;
        integPoints2D_[0][2] = 0.249286745170910;
        integPoints2D_[0][3] = 0.873821971016996;
        integPoints2D_[0][4] = 0.063089014491502;
        integPoints2D_[0][5] = 0.063089014491502;
        integPoints2D_[0][6] = 0.053145049844816;
        integPoints2D_[0][7] = 0.310352451033785;
        integPoints2D_[0][8] = 0.636502499121399;
        integPoints2D_[0][9] = 0.310352451033785;
        integPoints2D_[0][10] = 0.636502499121399;
        integPoints2D_[0][11] = 0.053145049844816;

        integPoints2D_[1][0] = 0.249286745170910;
        integPoints2D_[1][1] = 0.249286745170910;
        integPoints2D_[1][2] = 0.501426509658179;
        integPoints2D_[1][3] = 0.063089014491502;
        integPoints2D_[1][4] = 0.063089014491502;
        integPoints2D_[1][5] = 0.873821971016996;
        integPoints2D_[1][6] = 0.310352451033785;
        integPoints2D_[1][7] = 0.636502499121399;
        integPoints2D_[1][8] = 0.053145049844816;
        integPoints2D_[1][9] = 0.053145049844816;
        integPoints2D_[1][10] = 0.310352451033785;
        integPoints2D_[1][11] = 0.636502499121399;

        integPoints2D_[2][0] = 0.058393137863190;
        integPoints2D_[2][1] = 0.058393137863190;
        integPoints2D_[2][2] = 0.058393137863190;
        integPoints2D_[2][3] = 0.025422453185104;
        integPoints2D_[2][4] = 0.025422453185104;
        integPoints2D_[2][5] = 0.025422453185104;
        integPoints2D_[2][6] = 0.041425537809187;
        integPoints2D_[2][7] = 0.041425537809187;
        integPoints2D_[2][8] = 0.041425537809187;
        integPoints2D_[2][9] = 0.041425537809187;
        integPoints2D_[2][10] = 0.041425537809187;
        integPoints2D_[2][11] = 0.041425537809187;
    }
    else if (nPoints2D_ == 14)
    {
        integPoints2D_[0][0] = 0.500000000000000;
        integPoints2D_[0][1] = 0.151929760985185;
        integPoints2D_[0][2] = 0.848070239014815;
        integPoints2D_[0][3] = 0.500000000000000;
        integPoints2D_[0][4] = 0.500000000000000;
        integPoints2D_[0][5] = 0.705213096157672;
        integPoints2D_[0][6] = 0.294786903842327;
        integPoints2D_[0][7] = 0.166666666666667;
        integPoints2D_[0][8] = 0.050643253661728;
        integPoints2D_[0][9] = 0.398713492676543;
        integPoints2D_[0][10] = 0.050643253661728;
        integPoints2D_[0][11] = 0.235071032052557;
        integPoints2D_[0][12] = 0.235071032052557;
        integPoints2D_[0][13] = 0.029857935894885;

        integPoints2D_[1][0] = 0.166666666666667;
        integPoints2D_[1][1] = 0.050643253661729;
        integPoints2D_[1][2] = 0.050643253661729;
        integPoints2D_[1][3] = 0.398713492676544;
        integPoints2D_[1][4] = 0.029857935894885;
        integPoints2D_[1][5] = 0.235071032052557;
        integPoints2D_[1][6] = 0.235071032052557;
        integPoints2D_[1][7] = 0.500000000000000;
        integPoints2D_[1][8] = 0.151929760985185;
        integPoints2D_[1][9] = 0.500000000000000;
        integPoints2D_[1][10] = 0.848070239014816;
        integPoints2D_[1][11] = 0.294786903842327;
        integPoints2D_[1][12] = 0.705213096157672;
        integPoints2D_[1][13] = 0.500000000000000;

        integPoints2D_[2][0] = 0.056250000000000;
        integPoints2D_[2][1] = 0.031484795136207;
        integPoints2D_[2][2] = 0.031484795136207;
        integPoints2D_[2][3] = 0.031484795136207;
        integPoints2D_[2][4] = 0.033098538197126;
        integPoints2D_[2][5] = 0.033098538197126;
        integPoints2D_[2][6] = 0.033098538197126;
        integPoints2D_[2][7] = 0.056250000000000;
        integPoints2D_[2][8] = 0.031484795136207;
        integPoints2D_[2][9] = 0.031484795136207;
        integPoints2D_[2][10] = 0.031484795136207;
        integPoints2D_[2][11] = 0.033098538197126;
        integPoints2D_[2][12] = 0.033098538197126;
        integPoints2D_[2][13] = 0.033098538197126;
    }
    else if (nPoints2D_ == 28)
    {
        integPoints2D_[0][0] = 0.166666666666667;
        integPoints2D_[0][1] = 0.050643253661728;
        integPoints2D_[0][2] = 0.398713492676543;
        integPoints2D_[0][3] = 0.050643253661728;
        integPoints2D_[0][4] = 0.235071032052557;
        integPoints2D_[0][5] = 0.235071032052557;
        integPoints2D_[0][6] = 0.029857935894885;
        integPoints2D_[0][7] = 0.333333333333333;
        integPoints2D_[0][8] = 0.101286507323457;
        integPoints2D_[0][9] = 0.449356746338272;
        integPoints2D_[0][10] = 0.449356746338272;
        integPoints2D_[0][11] = 0.264928967947442;
        integPoints2D_[0][12] = 0.470142064105115;
        integPoints2D_[0][13] = 0.264928967947442;
        integPoints2D_[0][14] = 0.666666666666667;
        integPoints2D_[0][15] = 0.550643253661728;
        integPoints2D_[0][16] = 0.898713492676543;
        integPoints2D_[0][17] = 0.550643253661728;
        integPoints2D_[0][18] = 0.735071032052557;
        integPoints2D_[0][19] = 0.735071032052558;
        integPoints2D_[0][20] = 0.529857935894885;
        integPoints2D_[0][21] = 0.166666666666667;
        integPoints2D_[0][22] = 0.050643253661728;
        integPoints2D_[0][23] = 0.398713492676543;
        integPoints2D_[0][24] = 0.050643253661728;
        integPoints2D_[0][25] = 0.235071032052557;
        integPoints2D_[0][26] = 0.235071032052557;
        integPoints2D_[0][27] = 0.029857935894885;

        integPoints2D_[1][0] = 0.166666666666667;
        integPoints2D_[1][1] = 0.050643253661729;
        integPoints2D_[1][2] = 0.050643253661729;
        integPoints2D_[1][3] = 0.398713492676544;
        integPoints2D_[1][4] = 0.029857935894885;
        integPoints2D_[1][5] = 0.235071032052557;
        integPoints2D_[1][6] = 0.235071032052557;
        integPoints2D_[1][7] = 0.333333333333333;
        integPoints2D_[1][8] = 0.449356746338272;
        integPoints2D_[1][9] = 0.101286507323457;
        integPoints2D_[1][10] = 0.449356746338272;
        integPoints2D_[1][11] = 0.264928967947442;
        integPoints2D_[1][12] = 0.264928967947442;
        integPoints2D_[1][13] = 0.470142064105115;
        integPoints2D_[1][14] = 0.166666666666667;
        integPoints2D_[1][15] = 0.050643253661729;
        integPoints2D_[1][16] = 0.050643253661729;
        integPoints2D_[1][17] = 0.398713492676544;
        integPoints2D_[1][18] = 0.029857935894885;
        integPoints2D_[1][19] = 0.235071032052557;
        integPoints2D_[1][20] = 0.235071032052557;
        integPoints2D_[1][21] = 0.666666666666667;
        integPoints2D_[1][22] = 0.550643253661729;
        integPoints2D_[1][23] = 0.550643253661729;
        integPoints2D_[1][24] = 0.898713492676544;
        integPoints2D_[1][25] = 0.529857935894885;
        integPoints2D_[1][26] = 0.735071032052558;
        integPoints2D_[1][27] = 0.735071032052557;

        integPoints2D_[2][0] = 0.028125000000000;
        integPoints2D_[2][1] = 0.015742397568103;
        integPoints2D_[2][2] = 0.015742397568103;
        integPoints2D_[2][3] = 0.015742397568103;
        integPoints2D_[2][4] = 0.016549269098563;
        integPoints2D_[2][5] = 0.016549269098563;
        integPoints2D_[2][6] = 0.016549269098563;
        integPoints2D_[2][7] = 0.028125000000000;
        integPoints2D_[2][8] = 0.015742397568103;
        integPoints2D_[2][9] = 0.015742397568103;
        integPoints2D_[2][10] = 0.015742397568103;
        integPoints2D_[2][11] = 0.016549269098563;
        integPoints2D_[2][12] = 0.016549269098563;
        integPoints2D_[2][13] = 0.016549269098563;
        integPoints2D_[2][14] = 0.028125000000000;
        integPoints2D_[2][15] = 0.015742397568103;
        integPoints2D_[2][16] = 0.015742397568103;
        integPoints2D_[2][17] = 0.015742397568103;
        integPoints2D_[2][18] = 0.016549269098563;
        integPoints2D_[2][19] = 0.016549269098563;
        integPoints2D_[2][20] = 0.016549269098563;
        integPoints2D_[2][21] = 0.028125000000000;
        integPoints2D_[2][22] = 0.015742397568103;
        integPoints2D_[2][23] = 0.015742397568103;
        integPoints2D_[2][24] = 0.015742397568103;
        integPoints2D_[2][25] = 0.016549269098563;
        integPoints2D_[2][26] = 0.016549269098563;
        integPoints2D_[2][27] = 0.016549269098563;
    }
    else if (nPoints2D_ == 112)
    {
        integPoints2D_[0][0] = 0.083333333333333;
        integPoints2D_[0][1] = 0.025321626830864;
        integPoints2D_[0][2] = 0.199356746338272;
        integPoints2D_[0][3] = 0.025321626830864;
        integPoints2D_[0][4] = 0.117535516026279;
        integPoints2D_[0][5] = 0.117535516026279;
        integPoints2D_[0][6] = 0.014928967947443;
        integPoints2D_[0][7] = 0.166666666666667;
        integPoints2D_[0][8] = 0.050643253661728;
        integPoints2D_[0][9] = 0.224678373169136;
        integPoints2D_[0][10] = 0.224678373169136;
        integPoints2D_[0][11] = 0.132464483973721;
        integPoints2D_[0][12] = 0.235071032052557;
        integPoints2D_[0][13] = 0.132464483973721;
        integPoints2D_[0][14] = 0.333333333333333;
        integPoints2D_[0][15] = 0.275321626830864;
        integPoints2D_[0][16] = 0.449356746338272;
        integPoints2D_[0][17] = 0.275321626830864;
        integPoints2D_[0][18] = 0.367535516026279;
        integPoints2D_[0][19] = 0.367535516026279;
        integPoints2D_[0][20] = 0.264928967947442;
        integPoints2D_[0][21] = 0.416666666666667;
        integPoints2D_[0][22] = 0.300643253661728;
        integPoints2D_[0][23] = 0.474678373169136;
        integPoints2D_[0][24] = 0.474678373169136;
        integPoints2D_[0][25] = 0.382464483973721;
        integPoints2D_[0][26] = 0.485071032052558;
        integPoints2D_[0][27] = 0.382464483973721;
        integPoints2D_[0][28] = 0.583333333333333;
        integPoints2D_[0][29] = 0.525321626830864;
        integPoints2D_[0][30] = 0.699356746338272;
        integPoints2D_[0][31] = 0.525321626830864;
        integPoints2D_[0][32] = 0.617535516026279;
        integPoints2D_[0][33] = 0.617535516026279;
        integPoints2D_[0][34] = 0.514928967947442;
        integPoints2D_[0][35] = 0.666666666666667;
        integPoints2D_[0][36] = 0.550643253661728;
        integPoints2D_[0][37] = 0.724678373169136;
        integPoints2D_[0][38] = 0.724678373169136;
        integPoints2D_[0][39] = 0.632464483973721;
        integPoints2D_[0][40] = 0.735071032052558;
        integPoints2D_[0][41] = 0.632464483973721;
        integPoints2D_[0][42] = 0.833333333333333;
        integPoints2D_[0][43] = 0.775321626830864;
        integPoints2D_[0][44] = 0.949356746338272;
        integPoints2D_[0][45] = 0.775321626830864;
        integPoints2D_[0][46] = 0.867535516026279;
        integPoints2D_[0][47] = 0.867535516026279;
        integPoints2D_[0][48] = 0.764928967947442;
        integPoints2D_[0][49] = 0.083333333333333;
        integPoints2D_[0][50] = 0.025321626830864;
        integPoints2D_[0][51] = 0.199356746338272;
        integPoints2D_[0][52] = 0.025321626830864;
        integPoints2D_[0][53] = 0.117535516026279;
        integPoints2D_[0][54] = 0.117535516026279;
        integPoints2D_[0][55] = 0.014928967947443;
        integPoints2D_[0][56] = 0.166666666666667;
        integPoints2D_[0][57] = 0.050643253661728;
        integPoints2D_[0][58] = 0.224678373169136;
        integPoints2D_[0][59] = 0.224678373169136;
        integPoints2D_[0][60] = 0.132464483973721;
        integPoints2D_[0][61] = 0.235071032052557;
        integPoints2D_[0][62] = 0.132464483973721;
        integPoints2D_[0][63] = 0.333333333333333;
        integPoints2D_[0][64] = 0.275321626830864;
        integPoints2D_[0][65] = 0.449356746338272;
        integPoints2D_[0][66] = 0.275321626830864;
        integPoints2D_[0][67] = 0.367535516026279;
        integPoints2D_[0][68] = 0.367535516026279;
        integPoints2D_[0][69] = 0.264928967947442;
        integPoints2D_[0][70] = 0.416666666666667;
        integPoints2D_[0][71] = 0.300643253661728;
        integPoints2D_[0][72] = 0.474678373169136;
        integPoints2D_[0][73] = 0.474678373169136;
        integPoints2D_[0][74] = 0.382464483973721;
        integPoints2D_[0][75] = 0.485071032052558;
        integPoints2D_[0][76] = 0.382464483973721;
        integPoints2D_[0][77] = 0.583333333333333;
        integPoints2D_[0][78] = 0.525321626830864;
        integPoints2D_[0][79] = 0.699356746338272;
        integPoints2D_[0][80] = 0.525321626830864;
        integPoints2D_[0][81] = 0.617535516026279;
        integPoints2D_[0][82] = 0.617535516026279;
        integPoints2D_[0][83] = 0.514928967947442;
        integPoints2D_[0][84] = 0.083333333333333;
        integPoints2D_[0][85] = 0.025321626830864;
        integPoints2D_[0][86] = 0.199356746338272;
        integPoints2D_[0][87] = 0.025321626830864;
        integPoints2D_[0][88] = 0.117535516026279;
        integPoints2D_[0][89] = 0.117535516026279;
        integPoints2D_[0][90] = 0.014928967947443;
        integPoints2D_[0][91] = 0.166666666666667;
        integPoints2D_[0][92] = 0.050643253661728;
        integPoints2D_[0][93] = 0.224678373169136;
        integPoints2D_[0][94] = 0.224678373169136;
        integPoints2D_[0][95] = 0.132464483973721;
        integPoints2D_[0][96] = 0.235071032052557;
        integPoints2D_[0][97] = 0.132464483973721;
        integPoints2D_[0][98] = 0.333333333333333;
        integPoints2D_[0][99] = 0.275321626830864;
        integPoints2D_[0][100] = 0.449356746338272;
        integPoints2D_[0][101] = 0.275321626830864;
        integPoints2D_[0][102] = 0.367535516026279;
        integPoints2D_[0][103] = 0.367535516026279;
        integPoints2D_[0][104] = 0.264928967947442;
        integPoints2D_[0][105] = 0.083333333333333;
        integPoints2D_[0][106] = 0.025321626830864;
        integPoints2D_[0][107] = 0.199356746338272;
        integPoints2D_[0][108] = 0.025321626830864;
        integPoints2D_[0][109] = 0.117535516026279;
        integPoints2D_[0][110] = 0.117535516026279;
        integPoints2D_[0][111] = 0.014928967947443;

        integPoints2D_[1][0] = 0.166666666666667;
        integPoints2D_[1][1] = 0.050643253661728;
        integPoints2D_[1][2] = 0.224678373169136;
        integPoints2D_[1][3] = 0.224678373169136;
        integPoints2D_[1][4] = 0.132464483973721;
        integPoints2D_[1][5] = 0.235071032052557;
        integPoints2D_[1][6] = 0.132464483973721;
        integPoints2D_[1][7] = 0.083333333333333;
        integPoints2D_[1][8] = 0.025321626830864;
        integPoints2D_[1][9] = 0.025321626830864;
        integPoints2D_[1][10] = 0.199356746338272;
        integPoints2D_[1][11] = 0.014928967947443;
        integPoints2D_[1][12] = 0.117535516026279;
        integPoints2D_[1][13] = 0.117535516026279;
        integPoints2D_[1][14] = 0.083333333333334;
        integPoints2D_[1][15] = 0.025321626830864;
        integPoints2D_[1][16] = 0.025321626830864;
        integPoints2D_[1][17] = 0.199356746338272;
        integPoints2D_[1][18] = 0.014928967947443;
        integPoints2D_[1][19] = 0.117535516026279;
        integPoints2D_[1][20] = 0.117535516026279;
        integPoints2D_[1][21] = 0.166666666666667;
        integPoints2D_[1][22] = 0.224678373169136;
        integPoints2D_[1][23] = 0.050643253661728;
        integPoints2D_[1][24] = 0.224678373169136;
        integPoints2D_[1][25] = 0.132464483973721;
        integPoints2D_[1][26] = 0.132464483973721;
        integPoints2D_[1][27] = 0.235071032052557;
        integPoints2D_[1][28] = 0.166666666666667;
        integPoints2D_[1][29] = 0.050643253661728;
        integPoints2D_[1][30] = 0.224678373169136;
        integPoints2D_[1][31] = 0.224678373169136;
        integPoints2D_[1][32] = 0.132464483973721;
        integPoints2D_[1][33] = 0.235071032052557;
        integPoints2D_[1][34] = 0.132464483973721;
        integPoints2D_[1][35] = 0.083333333333334;
        integPoints2D_[1][36] = 0.025321626830864;
        integPoints2D_[1][37] = 0.025321626830864;
        integPoints2D_[1][38] = 0.199356746338272;
        integPoints2D_[1][39] = 0.014928967947443;
        integPoints2D_[1][40] = 0.117535516026279;
        integPoints2D_[1][41] = 0.117535516026279;
        integPoints2D_[1][42] = 0.083333333333334;
        integPoints2D_[1][43] = 0.025321626830864;
        integPoints2D_[1][44] = 0.025321626830864;
        integPoints2D_[1][45] = 0.199356746338272;
        integPoints2D_[1][46] = 0.014928967947443;
        integPoints2D_[1][47] = 0.117535516026279;
        integPoints2D_[1][48] = 0.117535516026279;
        integPoints2D_[1][49] = 0.333333333333333;
        integPoints2D_[1][50] = 0.275321626830864;
        integPoints2D_[1][51] = 0.275321626830864;
        integPoints2D_[1][52] = 0.449356746338272;
        integPoints2D_[1][53] = 0.264928967947442;
        integPoints2D_[1][54] = 0.367535516026279;
        integPoints2D_[1][55] = 0.367535516026279;
        integPoints2D_[1][56] = 0.416666666666667;
        integPoints2D_[1][57] = 0.474678373169136;
        integPoints2D_[1][58] = 0.300643253661728;
        integPoints2D_[1][59] = 0.474678373169136;
        integPoints2D_[1][60] = 0.382464483973721;
        integPoints2D_[1][61] = 0.382464483973721;
        integPoints2D_[1][62] = 0.485071032052557;
        integPoints2D_[1][63] = 0.416666666666667;
        integPoints2D_[1][64] = 0.300643253661728;
        integPoints2D_[1][65] = 0.474678373169136;
        integPoints2D_[1][66] = 0.474678373169136;
        integPoints2D_[1][67] = 0.382464483973721;
        integPoints2D_[1][68] = 0.485071032052558;
        integPoints2D_[1][69] = 0.382464483973721;
        integPoints2D_[1][70] = 0.333333333333333;
        integPoints2D_[1][71] = 0.275321626830864;
        integPoints2D_[1][72] = 0.275321626830864;
        integPoints2D_[1][73] = 0.449356746338272;
        integPoints2D_[1][74] = 0.264928967947442;
        integPoints2D_[1][75] = 0.367535516026279;
        integPoints2D_[1][76] = 0.367535516026279;
        integPoints2D_[1][77] = 0.333333333333333;
        integPoints2D_[1][78] = 0.275321626830864;
        integPoints2D_[1][79] = 0.275321626830864;
        integPoints2D_[1][80] = 0.449356746338272;
        integPoints2D_[1][81] = 0.264928967947442;
        integPoints2D_[1][82] = 0.367535516026279;
        integPoints2D_[1][83] = 0.367535516026279;
        integPoints2D_[1][84] = 0.666666666666667;
        integPoints2D_[1][85] = 0.550643253661728;
        integPoints2D_[1][86] = 0.724678373169136;
        integPoints2D_[1][87] = 0.724678373169136;
        integPoints2D_[1][88] = 0.632464483973721;
        integPoints2D_[1][89] = 0.735071032052558;
        integPoints2D_[1][90] = 0.632464483973721;
        integPoints2D_[1][91] = 0.583333333333333;
        integPoints2D_[1][92] = 0.525321626830864;
        integPoints2D_[1][93] = 0.525321626830864;
        integPoints2D_[1][94] = 0.699356746338272;
        integPoints2D_[1][95] = 0.514928967947442;
        integPoints2D_[1][96] = 0.617535516026279;
        integPoints2D_[1][97] = 0.617535516026279;
        integPoints2D_[1][98] = 0.583333333333333;
        integPoints2D_[1][99] = 0.525321626830864;
        integPoints2D_[1][100] = 0.525321626830864;
        integPoints2D_[1][101] = 0.699356746338272;
        integPoints2D_[1][102] = 0.514928967947442;
        integPoints2D_[1][103] = 0.617535516026279;
        integPoints2D_[1][104] = 0.617535516026279;
        integPoints2D_[1][105] = 0.833333333333333;
        integPoints2D_[1][106] = 0.775321626830864;
        integPoints2D_[1][107] = 0.775321626830864;
        integPoints2D_[1][108] = 0.949356746338272;
        integPoints2D_[1][109] = 0.764928967947442;
        integPoints2D_[1][110] = 0.867535516026279;
        integPoints2D_[1][111] = 0.867535516026279;

        integPoints2D_[2][0] = 0.007031250000000;
        integPoints2D_[2][1] = 0.003935599392026;
        integPoints2D_[2][2] = 0.003935599392026;
        integPoints2D_[2][3] = 0.003935599392026;
        integPoints2D_[2][4] = 0.004137317274641;
        integPoints2D_[2][5] = 0.004137317274641;
        integPoints2D_[2][6] = 0.004137317274641;
        integPoints2D_[2][7] = 0.007031250000000;
        integPoints2D_[2][8] = 0.003935599392026;
        integPoints2D_[2][9] = 0.003935599392026;
        integPoints2D_[2][10] = 0.003935599392026;
        integPoints2D_[2][11] = 0.004137317274641;
        integPoints2D_[2][12] = 0.004137317274641;
        integPoints2D_[2][13] = 0.004137317274641;
        integPoints2D_[2][14] = 0.007031250000000;
        integPoints2D_[2][15] = 0.003935599392026;
        integPoints2D_[2][16] = 0.003935599392026;
        integPoints2D_[2][17] = 0.003935599392026;
        integPoints2D_[2][18] = 0.004137317274641;
        integPoints2D_[2][19] = 0.004137317274641;
        integPoints2D_[2][20] = 0.004137317274641;
        integPoints2D_[2][21] = 0.007031250000000;
        integPoints2D_[2][22] = 0.003935599392026;
        integPoints2D_[2][23] = 0.003935599392026;
        integPoints2D_[2][24] = 0.003935599392026;
        integPoints2D_[2][25] = 0.004137317274641;
        integPoints2D_[2][26] = 0.004137317274641;
        integPoints2D_[2][27] = 0.004137317274641;
        integPoints2D_[2][28] = 0.007031250000000;
        integPoints2D_[2][29] = 0.003935599392026;
        integPoints2D_[2][30] = 0.003935599392026;
        integPoints2D_[2][31] = 0.003935599392026;
        integPoints2D_[2][32] = 0.004137317274641;
        integPoints2D_[2][33] = 0.004137317274641;
        integPoints2D_[2][34] = 0.004137317274641;
        integPoints2D_[2][35] = 0.007031250000000;
        integPoints2D_[2][36] = 0.003935599392026;
        integPoints2D_[2][37] = 0.003935599392026;
        integPoints2D_[2][38] = 0.003935599392026;
        integPoints2D_[2][39] = 0.004137317274641;
        integPoints2D_[2][40] = 0.004137317274641;
        integPoints2D_[2][41] = 0.004137317274641;
        integPoints2D_[2][42] = 0.007031250000000;
        integPoints2D_[2][43] = 0.003935599392026;
        integPoints2D_[2][44] = 0.003935599392026;
        integPoints2D_[2][45] = 0.003935599392026;
        integPoints2D_[2][46] = 0.004137317274641;
        integPoints2D_[2][47] = 0.004137317274641;
        integPoints2D_[2][48] = 0.004137317274641;
        integPoints2D_[2][49] = 0.007031250000000;
        integPoints2D_[2][50] = 0.003935599392026;
        integPoints2D_[2][51] = 0.003935599392026;
        integPoints2D_[2][52] = 0.003935599392026;
        integPoints2D_[2][53] = 0.004137317274641;
        integPoints2D_[2][54] = 0.004137317274641;
        integPoints2D_[2][55] = 0.004137317274641;
        integPoints2D_[2][56] = 0.007031250000000;
        integPoints2D_[2][57] = 0.003935599392026;
        integPoints2D_[2][58] = 0.003935599392026;
        integPoints2D_[2][59] = 0.003935599392026;
        integPoints2D_[2][60] = 0.004137317274641;
        integPoints2D_[2][61] = 0.004137317274641;
        integPoints2D_[2][62] = 0.004137317274641;
        integPoints2D_[2][63] = 0.007031250000000;
        integPoints2D_[2][64] = 0.003935599392026;
        integPoints2D_[2][65] = 0.003935599392026;
        integPoints2D_[2][66] = 0.003935599392026;
        integPoints2D_[2][67] = 0.004137317274641;
        integPoints2D_[2][68] = 0.004137317274641;
        integPoints2D_[2][69] = 0.004137317274641;
        integPoints2D_[2][70] = 0.007031250000000;
        integPoints2D_[2][71] = 0.003935599392026;
        integPoints2D_[2][72] = 0.003935599392026;
        integPoints2D_[2][73] = 0.003935599392026;
        integPoints2D_[2][74] = 0.004137317274641;
        integPoints2D_[2][75] = 0.004137317274641;
        integPoints2D_[2][76] = 0.004137317274641;
        integPoints2D_[2][77] = 0.007031250000000;
        integPoints2D_[2][78] = 0.003935599392026;
        integPoints2D_[2][79] = 0.003935599392026;
        integPoints2D_[2][80] = 0.003935599392026;
        integPoints2D_[2][81] = 0.004137317274641;
        integPoints2D_[2][82] = 0.004137317274641;
        integPoints2D_[2][83] = 0.004137317274641;
        integPoints2D_[2][84] = 0.007031250000000;
        integPoints2D_[2][85] = 0.003935599392026;
        integPoints2D_[2][86] = 0.003935599392026;
        integPoints2D_[2][87] = 0.003935599392026;
        integPoints2D_[2][88] = 0.004137317274641;
        integPoints2D_[2][89] = 0.004137317274641;
        integPoints2D_[2][90] = 0.004137317274641;
        integPoints2D_[2][91] = 0.007031250000000;
        integPoints2D_[2][92] = 0.003935599392026;
        integPoints2D_[2][93] = 0.003935599392026;
        integPoints2D_[2][94] = 0.003935599392026;
        integPoints2D_[2][95] = 0.004137317274641;
        integPoints2D_[2][96] = 0.004137317274641;
        integPoints2D_[2][97] = 0.004137317274641;
        integPoints2D_[2][98] = 0.007031250000000;
        integPoints2D_[2][99] = 0.003935599392026;
        integPoints2D_[2][100] = 0.003935599392026;
        integPoints2D_[2][101] = 0.003935599392026;
        integPoints2D_[2][102] = 0.004137317274641;
        integPoints2D_[2][103] = 0.004137317274641;
        integPoints2D_[2][104] = 0.004137317274641;
        integPoints2D_[2][105] = 0.007031250000000;
        integPoints2D_[2][106] = 0.003935599392026;
        integPoints2D_[2][107] = 0.003935599392026;
        integPoints2D_[2][108] = 0.003935599392026;
        integPoints2D_[2][109] = 0.004137317274641;
        integPoints2D_[2][110] = 0.004137317274641;
        integPoints2D_[2][111] = 0.004137317274641;
    }
    else if (nPoints2D_ == 224)
    {
        integPoints2D_[0][0] = 0.125000000000000;
        integPoints2D_[0][1] = 0.037982440246296;
        integPoints2D_[0][2] = 0.212017559753704;
        integPoints2D_[0][3] = 0.125000000000000;
        integPoints2D_[0][4] = 0.125000000000000;
        integPoints2D_[0][5] = 0.176303274039418;
        integPoints2D_[0][6] = 0.073696725960582;
        integPoints2D_[0][7] = 0.375000000000000;
        integPoints2D_[0][8] = 0.287982440246296;
        integPoints2D_[0][9] = 0.462017559753704;
        integPoints2D_[0][10] = 0.375000000000000;
        integPoints2D_[0][11] = 0.375000000000000;
        integPoints2D_[0][12] = 0.426303274039418;
        integPoints2D_[0][13] = 0.323696725960582;
        integPoints2D_[0][14] = 0.625000000000000;
        integPoints2D_[0][15] = 0.537982440246296;
        integPoints2D_[0][16] = 0.712017559753704;
        integPoints2D_[0][17] = 0.625000000000000;
        integPoints2D_[0][18] = 0.625000000000000;
        integPoints2D_[0][19] = 0.676303274039418;
        integPoints2D_[0][20] = 0.573696725960582;
        integPoints2D_[0][21] = 0.875000000000000;
        integPoints2D_[0][22] = 0.787982440246296;
        integPoints2D_[0][23] = 0.962017559753704;
        integPoints2D_[0][24] = 0.875000000000000;
        integPoints2D_[0][25] = 0.875000000000000;
        integPoints2D_[0][26] = 0.926303274039418;
        integPoints2D_[0][27] = 0.823696725960582;
        integPoints2D_[0][28] = 0.041666666666667;
        integPoints2D_[0][29] = 0.012660813415432;
        integPoints2D_[0][30] = 0.099678373169136;
        integPoints2D_[0][31] = 0.012660813415432;
        integPoints2D_[0][32] = 0.058767758013139;
        integPoints2D_[0][33] = 0.058767758013139;
        integPoints2D_[0][34] = 0.007464483973721;
        integPoints2D_[0][35] = 0.208333333333333;
        integPoints2D_[0][36] = 0.237339186584568;
        integPoints2D_[0][37] = 0.237339186584568;
        integPoints2D_[0][38] = 0.150321626830864;
        integPoints2D_[0][39] = 0.242535516026279;
        integPoints2D_[0][40] = 0.191232241986861;
        integPoints2D_[0][41] = 0.191232241986861;
        integPoints2D_[0][42] = 0.291666666666667;
        integPoints2D_[0][43] = 0.262660813415432;
        integPoints2D_[0][44] = 0.349678373169136;
        integPoints2D_[0][45] = 0.262660813415432;
        integPoints2D_[0][46] = 0.308767758013139;
        integPoints2D_[0][47] = 0.308767758013139;
        integPoints2D_[0][48] = 0.257464483973721;
        integPoints2D_[0][49] = 0.458333333333333;
        integPoints2D_[0][50] = 0.487339186584568;
        integPoints2D_[0][51] = 0.487339186584568;
        integPoints2D_[0][52] = 0.400321626830864;
        integPoints2D_[0][53] = 0.492535516026279;
        integPoints2D_[0][54] = 0.441232241986861;
        integPoints2D_[0][55] = 0.441232241986861;
        integPoints2D_[0][56] = 0.541666666666667;
        integPoints2D_[0][57] = 0.512660813415432;
        integPoints2D_[0][58] = 0.599678373169136;
        integPoints2D_[0][59] = 0.512660813415432;
        integPoints2D_[0][60] = 0.558767758013139;
        integPoints2D_[0][61] = 0.558767758013139;
        integPoints2D_[0][62] = 0.507464483973721;
        integPoints2D_[0][63] = 0.708333333333333;
        integPoints2D_[0][64] = 0.737339186584568;
        integPoints2D_[0][65] = 0.737339186584568;
        integPoints2D_[0][66] = 0.650321626830864;
        integPoints2D_[0][67] = 0.742535516026279;
        integPoints2D_[0][68] = 0.691232241986861;
        integPoints2D_[0][69] = 0.691232241986861;
        integPoints2D_[0][70] = 0.791666666666667;
        integPoints2D_[0][71] = 0.762660813415432;
        integPoints2D_[0][72] = 0.849678373169136;
        integPoints2D_[0][73] = 0.762660813415432;
        integPoints2D_[0][74] = 0.808767758013139;
        integPoints2D_[0][75] = 0.808767758013139;
        integPoints2D_[0][76] = 0.757464483973721;
        integPoints2D_[0][77] = 0.125000000000000;
        integPoints2D_[0][78] = 0.037982440246296;
        integPoints2D_[0][79] = 0.125000000000000;
        integPoints2D_[0][80] = 0.212017559753704;
        integPoints2D_[0][81] = 0.073696725960582;
        integPoints2D_[0][82] = 0.176303274039418;
        integPoints2D_[0][83] = 0.125000000000000;
        integPoints2D_[0][84] = 0.375000000000000;
        integPoints2D_[0][85] = 0.287982440246296;
        integPoints2D_[0][86] = 0.375000000000000;
        integPoints2D_[0][87] = 0.462017559753704;
        integPoints2D_[0][88] = 0.323696725960582;
        integPoints2D_[0][89] = 0.426303274039418;
        integPoints2D_[0][90] = 0.375000000000000;
        integPoints2D_[0][91] = 0.541666666666667;
        integPoints2D_[0][92] = 0.338625693908024;
        integPoints2D_[0][93] = 0.599678373169136;
        integPoints2D_[0][94] = 0.686695932922840;
        integPoints2D_[0][95] = 0.456161209934303;
        integPoints2D_[0][96] = 0.661374306091976;
        integPoints2D_[0][97] = 0.507464483973721;
        integPoints2D_[0][98] = 0.125000000000000;
        integPoints2D_[0][99] = 0.037982440246296;
        integPoints2D_[0][100] = 0.212017559753704;
        integPoints2D_[0][101] = 0.125000000000000;
        integPoints2D_[0][102] = 0.125000000000000;
        integPoints2D_[0][103] = 0.176303274039418;
        integPoints2D_[0][104] = 0.073696725960582;
        integPoints2D_[0][105] = 0.375000000000000;
        integPoints2D_[0][106] = 0.287982440246296;
        integPoints2D_[0][107] = 0.462017559753704;
        integPoints2D_[0][108] = 0.375000000000000;
        integPoints2D_[0][109] = 0.375000000000000;
        integPoints2D_[0][110] = 0.426303274039418;
        integPoints2D_[0][111] = 0.323696725960582;
        integPoints2D_[0][112] = 0.625000000000000;
        integPoints2D_[0][113] = 0.537982440246296;
        integPoints2D_[0][114] = 0.712017559753704;
        integPoints2D_[0][115] = 0.625000000000000;
        integPoints2D_[0][116] = 0.625000000000000;
        integPoints2D_[0][117] = 0.676303274039418;
        integPoints2D_[0][118] = 0.573696725960582;
        integPoints2D_[0][119] = 0.041666666666667;
        integPoints2D_[0][120] = 0.012660813415432;
        integPoints2D_[0][121] = 0.099678373169136;
        integPoints2D_[0][122] = 0.012660813415432;
        integPoints2D_[0][123] = 0.058767758013139;
        integPoints2D_[0][124] = 0.058767758013139;
        integPoints2D_[0][125] = 0.007464483973721;
        integPoints2D_[0][126] = 0.208333333333333;
        integPoints2D_[0][127] = 0.150321626830864;
        integPoints2D_[0][128] = 0.237339186584568;
        integPoints2D_[0][129] = 0.237339186584568;
        integPoints2D_[0][130] = 0.191232241986861;
        integPoints2D_[0][131] = 0.242535516026279;
        integPoints2D_[0][132] = 0.191232241986861;
        integPoints2D_[0][133] = 0.291666666666667;
        integPoints2D_[0][134] = 0.262660813415432;
        integPoints2D_[0][135] = 0.349678373169136;
        integPoints2D_[0][136] = 0.262660813415432;
        integPoints2D_[0][137] = 0.308767758013139;
        integPoints2D_[0][138] = 0.308767758013139;
        integPoints2D_[0][139] = 0.257464483973721;
        integPoints2D_[0][140] = 0.458333333333333;
        integPoints2D_[0][141] = 0.400321626830864;
        integPoints2D_[0][142] = 0.487339186584568;
        integPoints2D_[0][143] = 0.487339186584568;
        integPoints2D_[0][144] = 0.441232241986861;
        integPoints2D_[0][145] = 0.492535516026279;
        integPoints2D_[0][146] = 0.441232241986861;
        integPoints2D_[0][147] = 0.541666666666667;
        integPoints2D_[0][148] = 0.512660813415432;
        integPoints2D_[0][149] = 0.599678373169136;
        integPoints2D_[0][150] = 0.512660813415432;
        integPoints2D_[0][151] = 0.558767758013139;
        integPoints2D_[0][152] = 0.558767758013139;
        integPoints2D_[0][153] = 0.507464483973721;
        integPoints2D_[0][154] = 0.125000000000000;
        integPoints2D_[0][155] = 0.037982440246296;
        integPoints2D_[0][156] = 0.125000000000000;
        integPoints2D_[0][157] = 0.212017559753704;
        integPoints2D_[0][158] = 0.073696725960582;
        integPoints2D_[0][159] = 0.176303274039418;
        integPoints2D_[0][160] = 0.125000000000000;
        integPoints2D_[0][161] = 0.375000000000000;
        integPoints2D_[0][162] = 0.287982440246296;
        integPoints2D_[0][163] = 0.375000000000000;
        integPoints2D_[0][164] = 0.462017559753704;
        integPoints2D_[0][165] = 0.323696725960582;
        integPoints2D_[0][166] = 0.426303274039418;
        integPoints2D_[0][167] = 0.375000000000000;
        integPoints2D_[0][168] = 0.125000000000000;
        integPoints2D_[0][169] = 0.037982440246296;
        integPoints2D_[0][170] = 0.212017559753704;
        integPoints2D_[0][171] = 0.125000000000000;
        integPoints2D_[0][172] = 0.125000000000000;
        integPoints2D_[0][173] = 0.176303274039418;
        integPoints2D_[0][174] = 0.073696725960582;
        integPoints2D_[0][175] = 0.375000000000000;
        integPoints2D_[0][176] = 0.287982440246296;
        integPoints2D_[0][177] = 0.462017559753704;
        integPoints2D_[0][178] = 0.375000000000000;
        integPoints2D_[0][179] = 0.375000000000000;
        integPoints2D_[0][180] = 0.426303274039418;
        integPoints2D_[0][181] = 0.323696725960582;
        integPoints2D_[0][182] = 0.041666666666667;
        integPoints2D_[0][183] = 0.012660813415432;
        integPoints2D_[0][184] = 0.099678373169136;
        integPoints2D_[0][185] = 0.012660813415432;
        integPoints2D_[0][186] = 0.058767758013139;
        integPoints2D_[0][187] = 0.058767758013139;
        integPoints2D_[0][188] = 0.007464483973721;
        integPoints2D_[0][189] = 0.208333333333333;
        integPoints2D_[0][190] = 0.150321626830864;
        integPoints2D_[0][191] = 0.237339186584568;
        integPoints2D_[0][192] = 0.237339186584568;
        integPoints2D_[0][193] = 0.191232241986861;
        integPoints2D_[0][194] = 0.242535516026279;
        integPoints2D_[0][195] = 0.191232241986861;
        integPoints2D_[0][196] = 0.291666666666667;
        integPoints2D_[0][197] = 0.262660813415432;
        integPoints2D_[0][198] = 0.349678373169136;
        integPoints2D_[0][199] = 0.262660813415432;
        integPoints2D_[0][200] = 0.308767758013139;
        integPoints2D_[0][201] = 0.308767758013139;
        integPoints2D_[0][202] = 0.257464483973721;
        integPoints2D_[0][203] = 0.125000000000000;
        integPoints2D_[0][204] = 0.037982440246296;
        integPoints2D_[0][205] = 0.125000000000000;
        integPoints2D_[0][206] = 0.212017559753704;
        integPoints2D_[0][207] = 0.073696725960582;
        integPoints2D_[0][208] = 0.176303274039418;
        integPoints2D_[0][209] = 0.125000000000000;
        integPoints2D_[0][210] = 0.125000000000000;
        integPoints2D_[0][211] = 0.037982440246296;
        integPoints2D_[0][212] = 0.212017559753704;
        integPoints2D_[0][213] = 0.125000000000000;
        integPoints2D_[0][214] = 0.125000000000000;
        integPoints2D_[0][215] = 0.176303274039418;
        integPoints2D_[0][216] = 0.073696725960582;
        integPoints2D_[0][217] = 0.041666666666667;
        integPoints2D_[0][218] = 0.012660813415432;
        integPoints2D_[0][219] = 0.099678373169136;
        integPoints2D_[0][220] = 0.012660813415432;
        integPoints2D_[0][221] = 0.058767758013139;
        integPoints2D_[0][222] = 0.058767758013139;
        integPoints2D_[0][223] = 0.007464483973721;

        integPoints2D_[1][0] = 0.041666666666667;
        integPoints2D_[1][1] = 0.012660813415432;
        integPoints2D_[1][2] = 0.012660813415432;
        integPoints2D_[1][3] = 0.099678373169136;
        integPoints2D_[1][4] = 0.007464483973721;
        integPoints2D_[1][5] = 0.058767758013139;
        integPoints2D_[1][6] = 0.058767758013139;
        integPoints2D_[1][7] = 0.041666666666667;
        integPoints2D_[1][8] = 0.012660813415432;
        integPoints2D_[1][9] = 0.012660813415432;
        integPoints2D_[1][10] = 0.099678373169136;
        integPoints2D_[1][11] = 0.007464483973721;
        integPoints2D_[1][12] = 0.058767758013139;
        integPoints2D_[1][13] = 0.058767758013139;
        integPoints2D_[1][14] = 0.041666666666667;
        integPoints2D_[1][15] = 0.012660813415432;
        integPoints2D_[1][16] = 0.012660813415432;
        integPoints2D_[1][17] = 0.099678373169136;
        integPoints2D_[1][18] = 0.007464483973721;
        integPoints2D_[1][19] = 0.058767758013139;
        integPoints2D_[1][20] = 0.058767758013139;
        integPoints2D_[1][21] = 0.041666666666667;
        integPoints2D_[1][22] = 0.012660813415432;
        integPoints2D_[1][23] = 0.012660813415432;
        integPoints2D_[1][24] = 0.099678373169136;
        integPoints2D_[1][25] = 0.007464483973721;
        integPoints2D_[1][26] = 0.058767758013139;
        integPoints2D_[1][27] = 0.058767758013139;
        integPoints2D_[1][28] = 0.125000000000000;
        integPoints2D_[1][29] = 0.037982440246296;
        integPoints2D_[1][30] = 0.125000000000000;
        integPoints2D_[1][31] = 0.212017559753704;
        integPoints2D_[1][32] = 0.073696725960582;
        integPoints2D_[1][33] = 0.176303274039418;
        integPoints2D_[1][34] = 0.125000000000000;
        integPoints2D_[1][35] = 0.125000000000000;
        integPoints2D_[1][36] = 0.037982440246296;
        integPoints2D_[1][37] = 0.212017559753704;
        integPoints2D_[1][38] = 0.125000000000000;
        integPoints2D_[1][39] = 0.125000000000000;
        integPoints2D_[1][40] = 0.176303274039418;
        integPoints2D_[1][41] = 0.073696725960582;
        integPoints2D_[1][42] = 0.125000000000000;
        integPoints2D_[1][43] = 0.037982440246296;
        integPoints2D_[1][44] = 0.125000000000000;
        integPoints2D_[1][45] = 0.212017559753704;
        integPoints2D_[1][46] = 0.073696725960582;
        integPoints2D_[1][47] = 0.176303274039418;
        integPoints2D_[1][48] = 0.125000000000000;
        integPoints2D_[1][49] = 0.125000000000000;
        integPoints2D_[1][50] = 0.037982440246296;
        integPoints2D_[1][51] = 0.212017559753704;
        integPoints2D_[1][52] = 0.125000000000000;
        integPoints2D_[1][53] = 0.125000000000000;
        integPoints2D_[1][54] = 0.176303274039418;
        integPoints2D_[1][55] = 0.073696725960582;
        integPoints2D_[1][56] = 0.125000000000000;
        integPoints2D_[1][57] = 0.037982440246296;
        integPoints2D_[1][58] = 0.125000000000000;
        integPoints2D_[1][59] = 0.212017559753704;
        integPoints2D_[1][60] = 0.073696725960582;
        integPoints2D_[1][61] = 0.176303274039418;
        integPoints2D_[1][62] = 0.125000000000000;
        integPoints2D_[1][63] = 0.125000000000000;
        integPoints2D_[1][64] = 0.037982440246296;
        integPoints2D_[1][65] = 0.212017559753704;
        integPoints2D_[1][66] = 0.125000000000000;
        integPoints2D_[1][67] = 0.125000000000000;
        integPoints2D_[1][68] = 0.176303274039418;
        integPoints2D_[1][69] = 0.073696725960582;
        integPoints2D_[1][70] = 0.125000000000000;
        integPoints2D_[1][71] = 0.037982440246296;
        integPoints2D_[1][72] = 0.125000000000000;
        integPoints2D_[1][73] = 0.212017559753704;
        integPoints2D_[1][74] = 0.073696725960582;
        integPoints2D_[1][75] = 0.176303274039418;
        integPoints2D_[1][76] = 0.125000000000000;
        integPoints2D_[1][77] = 0.208333333333333;
        integPoints2D_[1][78] = 0.237339186584568;
        integPoints2D_[1][79] = 0.150321626830864;
        integPoints2D_[1][80] = 0.237339186584568;
        integPoints2D_[1][81] = 0.191232241986861;
        integPoints2D_[1][82] = 0.191232241986861;
        integPoints2D_[1][83] = 0.242535516026279;
        integPoints2D_[1][84] = 0.208333333333333;
        integPoints2D_[1][85] = 0.237339186584568;
        integPoints2D_[1][86] = 0.150321626830864;
        integPoints2D_[1][87] = 0.237339186584568;
        integPoints2D_[1][88] = 0.191232241986861;
        integPoints2D_[1][89] = 0.191232241986861;
        integPoints2D_[1][90] = 0.242535516026279;
        integPoints2D_[1][91] = 0.291666666666667;
        integPoints2D_[1][92] = 0.436695932922840;
        integPoints2D_[1][93] = 0.175643253661728;
        integPoints2D_[1][94] = 0.262660813415432;
        integPoints2D_[1][95] = 0.308767758013139;
        integPoints2D_[1][96] = 0.206161209934303;
        integPoints2D_[1][97] = 0.360071032052558;
        integPoints2D_[1][98] = 0.291666666666667;
        integPoints2D_[1][99] = 0.262660813415432;
        integPoints2D_[1][100] = 0.262660813415432;
        integPoints2D_[1][101] = 0.349678373169136;
        integPoints2D_[1][102] = 0.257464483973721;
        integPoints2D_[1][103] = 0.308767758013139;
        integPoints2D_[1][104] = 0.308767758013139;
        integPoints2D_[1][105] = 0.291666666666667;
        integPoints2D_[1][106] = 0.262660813415432;
        integPoints2D_[1][107] = 0.262660813415432;
        integPoints2D_[1][108] = 0.349678373169136;
        integPoints2D_[1][109] = 0.257464483973721;
        integPoints2D_[1][110] = 0.308767758013139;
        integPoints2D_[1][111] = 0.308767758013139;
        integPoints2D_[1][112] = 0.291666666666667;
        integPoints2D_[1][113] = 0.262660813415432;
        integPoints2D_[1][114] = 0.262660813415432;
        integPoints2D_[1][115] = 0.349678373169136;
        integPoints2D_[1][116] = 0.257464483973721;
        integPoints2D_[1][117] = 0.308767758013139;
        integPoints2D_[1][118] = 0.308767758013139;
        integPoints2D_[1][119] = 0.375000000000000;
        integPoints2D_[1][120] = 0.287982440246296;
        integPoints2D_[1][121] = 0.375000000000000;
        integPoints2D_[1][122] = 0.462017559753704;
        integPoints2D_[1][123] = 0.323696725960582;
        integPoints2D_[1][124] = 0.426303274039418;
        integPoints2D_[1][125] = 0.375000000000000;
        integPoints2D_[1][126] = 0.375000000000000;
        integPoints2D_[1][127] = 0.375000000000000;
        integPoints2D_[1][128] = 0.287982440246296;
        integPoints2D_[1][129] = 0.462017559753704;
        integPoints2D_[1][130] = 0.323696725960582;
        integPoints2D_[1][131] = 0.375000000000000;
        integPoints2D_[1][132] = 0.426303274039418;
        integPoints2D_[1][133] = 0.375000000000000;
        integPoints2D_[1][134] = 0.287982440246296;
        integPoints2D_[1][135] = 0.375000000000000;
        integPoints2D_[1][136] = 0.462017559753704;
        integPoints2D_[1][137] = 0.323696725960582;
        integPoints2D_[1][138] = 0.426303274039418;
        integPoints2D_[1][139] = 0.375000000000000;
        integPoints2D_[1][140] = 0.375000000000000;
        integPoints2D_[1][141] = 0.375000000000000;
        integPoints2D_[1][142] = 0.287982440246296;
        integPoints2D_[1][143] = 0.462017559753704;
        integPoints2D_[1][144] = 0.323696725960582;
        integPoints2D_[1][145] = 0.375000000000000;
        integPoints2D_[1][146] = 0.426303274039418;
        integPoints2D_[1][147] = 0.375000000000000;
        integPoints2D_[1][148] = 0.287982440246296;
        integPoints2D_[1][149] = 0.375000000000000;
        integPoints2D_[1][150] = 0.462017559753704;
        integPoints2D_[1][151] = 0.323696725960582;
        integPoints2D_[1][152] = 0.426303274039418;
        integPoints2D_[1][153] = 0.375000000000000;
        integPoints2D_[1][154] = 0.458333333333333;
        integPoints2D_[1][155] = 0.487339186584568;
        integPoints2D_[1][156] = 0.400321626830864;
        integPoints2D_[1][157] = 0.487339186584568;
        integPoints2D_[1][158] = 0.441232241986861;
        integPoints2D_[1][159] = 0.441232241986861;
        integPoints2D_[1][160] = 0.492535516026279;
        integPoints2D_[1][161] = 0.458333333333333;
        integPoints2D_[1][162] = 0.487339186584568;
        integPoints2D_[1][163] = 0.400321626830864;
        integPoints2D_[1][164] = 0.487339186584568;
        integPoints2D_[1][165] = 0.441232241986861;
        integPoints2D_[1][166] = 0.441232241986861;
        integPoints2D_[1][167] = 0.492535516026279;
        integPoints2D_[1][168] = 0.541666666666667;
        integPoints2D_[1][169] = 0.512660813415432;
        integPoints2D_[1][170] = 0.512660813415432;
        integPoints2D_[1][171] = 0.599678373169136;
        integPoints2D_[1][172] = 0.507464483973721;
        integPoints2D_[1][173] = 0.558767758013139;
        integPoints2D_[1][174] = 0.558767758013139;
        integPoints2D_[1][175] = 0.541666666666667;
        integPoints2D_[1][176] = 0.512660813415432;
        integPoints2D_[1][177] = 0.512660813415432;
        integPoints2D_[1][178] = 0.599678373169136;
        integPoints2D_[1][179] = 0.507464483973721;
        integPoints2D_[1][180] = 0.558767758013139;
        integPoints2D_[1][181] = 0.558767758013139;
        integPoints2D_[1][182] = 0.625000000000000;
        integPoints2D_[1][183] = 0.537982440246296;
        integPoints2D_[1][184] = 0.625000000000000;
        integPoints2D_[1][185] = 0.712017559753704;
        integPoints2D_[1][186] = 0.573696725960582;
        integPoints2D_[1][187] = 0.676303274039418;
        integPoints2D_[1][188] = 0.625000000000000;
        integPoints2D_[1][189] = 0.625000000000000;
        integPoints2D_[1][190] = 0.625000000000000;
        integPoints2D_[1][191] = 0.537982440246296;
        integPoints2D_[1][192] = 0.712017559753704;
        integPoints2D_[1][193] = 0.573696725960582;
        integPoints2D_[1][194] = 0.625000000000000;
        integPoints2D_[1][195] = 0.676303274039418;
        integPoints2D_[1][196] = 0.625000000000000;
        integPoints2D_[1][197] = 0.537982440246296;
        integPoints2D_[1][198] = 0.625000000000000;
        integPoints2D_[1][199] = 0.712017559753704;
        integPoints2D_[1][200] = 0.573696725960582;
        integPoints2D_[1][201] = 0.676303274039418;
        integPoints2D_[1][202] = 0.625000000000000;
        integPoints2D_[1][203] = 0.708333333333333;
        integPoints2D_[1][204] = 0.737339186584568;
        integPoints2D_[1][205] = 0.650321626830864;
        integPoints2D_[1][206] = 0.737339186584568;
        integPoints2D_[1][207] = 0.691232241986861;
        integPoints2D_[1][208] = 0.691232241986861;
        integPoints2D_[1][209] = 0.742535516026279;
        integPoints2D_[1][210] = 0.791666666666667;
        integPoints2D_[1][211] = 0.762660813415432;
        integPoints2D_[1][212] = 0.762660813415432;
        integPoints2D_[1][213] = 0.849678373169136;
        integPoints2D_[1][214] = 0.757464483973721;
        integPoints2D_[1][215] = 0.808767758013139;
        integPoints2D_[1][216] = 0.808767758013139;
        integPoints2D_[1][217] = 0.875000000000000;
        integPoints2D_[1][218] = 0.787982440246296;
        integPoints2D_[1][219] = 0.875000000000000;
        integPoints2D_[1][220] = 0.962017559753704;
        integPoints2D_[1][221] = 0.823696725960582;
        integPoints2D_[1][222] = 0.926303274039418;
        integPoints2D_[1][223] = 0.875000000000000;

        integPoints2D_[2][0] = 0.003515625000000;
        integPoints2D_[2][1] = 0.001967799696013;
        integPoints2D_[2][2] = 0.001967799696013;
        integPoints2D_[2][3] = 0.001967799696013;
        integPoints2D_[2][4] = 0.002068658637320;
        integPoints2D_[2][5] = 0.002068658637320;
        integPoints2D_[2][6] = 0.002068658637320;
        integPoints2D_[2][7] = 0.003515625000000;
        integPoints2D_[2][8] = 0.001967799696013;
        integPoints2D_[2][9] = 0.001967799696013;
        integPoints2D_[2][10] = 0.001967799696013;
        integPoints2D_[2][11] = 0.002068658637320;
        integPoints2D_[2][12] = 0.002068658637320;
        integPoints2D_[2][13] = 0.002068658637320;
        integPoints2D_[2][14] = 0.003515625000000;
        integPoints2D_[2][15] = 0.001967799696013;
        integPoints2D_[2][16] = 0.001967799696013;
        integPoints2D_[2][17] = 0.001967799696013;
        integPoints2D_[2][18] = 0.002068658637320;
        integPoints2D_[2][19] = 0.002068658637320;
        integPoints2D_[2][20] = 0.002068658637320;
        integPoints2D_[2][21] = 0.003515625000000;
        integPoints2D_[2][22] = 0.001967799696013;
        integPoints2D_[2][23] = 0.001967799696013;
        integPoints2D_[2][24] = 0.001967799696013;
        integPoints2D_[2][25] = 0.002068658637320;
        integPoints2D_[2][26] = 0.002068658637320;
        integPoints2D_[2][27] = 0.002068658637320;
        integPoints2D_[2][28] = 0.003515625000000;
        integPoints2D_[2][29] = 0.001967799696013;
        integPoints2D_[2][30] = 0.001967799696013;
        integPoints2D_[2][31] = 0.001967799696013;
        integPoints2D_[2][32] = 0.002068658637320;
        integPoints2D_[2][33] = 0.002068658637320;
        integPoints2D_[2][34] = 0.002068658637320;
        integPoints2D_[2][35] = 0.003515625000000;
        integPoints2D_[2][36] = 0.001967799696013;
        integPoints2D_[2][37] = 0.001967799696013;
        integPoints2D_[2][38] = 0.001967799696013;
        integPoints2D_[2][39] = 0.002068658637320;
        integPoints2D_[2][40] = 0.002068658637320;
        integPoints2D_[2][41] = 0.002068658637320;
        integPoints2D_[2][42] = 0.003515625000000;
        integPoints2D_[2][43] = 0.001967799696013;
        integPoints2D_[2][44] = 0.001967799696013;
        integPoints2D_[2][45] = 0.001967799696013;
        integPoints2D_[2][46] = 0.002068658637320;
        integPoints2D_[2][47] = 0.002068658637320;
        integPoints2D_[2][48] = 0.002068658637320;
        integPoints2D_[2][49] = 0.003515625000000;
        integPoints2D_[2][50] = 0.001967799696013;
        integPoints2D_[2][51] = 0.001967799696013;
        integPoints2D_[2][52] = 0.001967799696013;
        integPoints2D_[2][53] = 0.002068658637320;
        integPoints2D_[2][54] = 0.002068658637320;
        integPoints2D_[2][55] = 0.002068658637320;
        integPoints2D_[2][56] = 0.003515625000000;
        integPoints2D_[2][57] = 0.001967799696013;
        integPoints2D_[2][58] = 0.001967799696013;
        integPoints2D_[2][59] = 0.001967799696013;
        integPoints2D_[2][60] = 0.002068658637320;
        integPoints2D_[2][61] = 0.002068658637320;
        integPoints2D_[2][62] = 0.002068658637320;
        integPoints2D_[2][63] = 0.003515625000000;
        integPoints2D_[2][64] = 0.001967799696013;
        integPoints2D_[2][65] = 0.001967799696013;
        integPoints2D_[2][66] = 0.001967799696013;
        integPoints2D_[2][67] = 0.002068658637320;
        integPoints2D_[2][68] = 0.002068658637320;
        integPoints2D_[2][69] = 0.002068658637320;
        integPoints2D_[2][70] = 0.003515625000000;
        integPoints2D_[2][71] = 0.001967799696013;
        integPoints2D_[2][72] = 0.001967799696013;
        integPoints2D_[2][73] = 0.001967799696013;
        integPoints2D_[2][74] = 0.002068658637320;
        integPoints2D_[2][75] = 0.002068658637320;
        integPoints2D_[2][76] = 0.002068658637320;
        integPoints2D_[2][77] = 0.003515625000000;
        integPoints2D_[2][78] = 0.001967799696013;
        integPoints2D_[2][79] = 0.001967799696013;
        integPoints2D_[2][80] = 0.001967799696013;
        integPoints2D_[2][81] = 0.002068658637320;
        integPoints2D_[2][82] = 0.002068658637320;
        integPoints2D_[2][83] = 0.002068658637320;
        integPoints2D_[2][84] = 0.003515625000000;
        integPoints2D_[2][85] = 0.001967799696013;
        integPoints2D_[2][86] = 0.001967799696013;
        integPoints2D_[2][87] = 0.001967799696013;
        integPoints2D_[2][88] = 0.002068658637320;
        integPoints2D_[2][89] = 0.002068658637320;
        integPoints2D_[2][90] = 0.002068658637320;
        integPoints2D_[2][91] = 0.003515625000000;
        integPoints2D_[2][92] = 0.001967799696013;
        integPoints2D_[2][93] = 0.001967799696013;
        integPoints2D_[2][94] = 0.001967799696013;
        integPoints2D_[2][95] = 0.002068658637320;
        integPoints2D_[2][96] = 0.002068658637320;
        integPoints2D_[2][97] = 0.002068658637320;
        integPoints2D_[2][98] = 0.003515625000000;
        integPoints2D_[2][99] = 0.001967799696013;
        integPoints2D_[2][100] = 0.001967799696013;
        integPoints2D_[2][101] = 0.001967799696013;
        integPoints2D_[2][102] = 0.002068658637320;
        integPoints2D_[2][103] = 0.002068658637320;
        integPoints2D_[2][104] = 0.002068658637320;
        integPoints2D_[2][105] = 0.003515625000000;
        integPoints2D_[2][106] = 0.001967799696013;
        integPoints2D_[2][107] = 0.001967799696013;
        integPoints2D_[2][108] = 0.001967799696013;
        integPoints2D_[2][109] = 0.002068658637320;
        integPoints2D_[2][110] = 0.002068658637320;
        integPoints2D_[2][111] = 0.002068658637320;
        integPoints2D_[2][112] = 0.003515625000000;
        integPoints2D_[2][113] = 0.001967799696013;
        integPoints2D_[2][114] = 0.001967799696013;
        integPoints2D_[2][115] = 0.001967799696013;
        integPoints2D_[2][116] = 0.002068658637320;
        integPoints2D_[2][117] = 0.002068658637320;
        integPoints2D_[2][118] = 0.002068658637320;
        integPoints2D_[2][119] = 0.003515625000000;
        integPoints2D_[2][120] = 0.001967799696013;
        integPoints2D_[2][121] = 0.001967799696013;
        integPoints2D_[2][122] = 0.001967799696013;
        integPoints2D_[2][123] = 0.002068658637320;
        integPoints2D_[2][124] = 0.002068658637320;
        integPoints2D_[2][125] = 0.002068658637320;
        integPoints2D_[2][126] = 0.003515625000000;
        integPoints2D_[2][127] = 0.001967799696013;
        integPoints2D_[2][128] = 0.001967799696013;
        integPoints2D_[2][129] = 0.001967799696013;
        integPoints2D_[2][130] = 0.002068658637320;
        integPoints2D_[2][131] = 0.002068658637320;
        integPoints2D_[2][132] = 0.002068658637320;
        integPoints2D_[2][133] = 0.003515625000000;
        integPoints2D_[2][134] = 0.001967799696013;
        integPoints2D_[2][135] = 0.001967799696013;
        integPoints2D_[2][136] = 0.001967799696013;
        integPoints2D_[2][137] = 0.002068658637320;
        integPoints2D_[2][138] = 0.002068658637320;
        integPoints2D_[2][139] = 0.002068658637320;
        integPoints2D_[2][140] = 0.003515625000000;
        integPoints2D_[2][141] = 0.001967799696013;
        integPoints2D_[2][142] = 0.001967799696013;
        integPoints2D_[2][143] = 0.001967799696013;
        integPoints2D_[2][144] = 0.002068658637320;
        integPoints2D_[2][145] = 0.002068658637320;
        integPoints2D_[2][146] = 0.002068658637320;
        integPoints2D_[2][147] = 0.003515625000000;
        integPoints2D_[2][148] = 0.001967799696013;
        integPoints2D_[2][149] = 0.001967799696013;
        integPoints2D_[2][150] = 0.001967799696013;
        integPoints2D_[2][151] = 0.002068658637320;
        integPoints2D_[2][152] = 0.002068658637320;
        integPoints2D_[2][153] = 0.002068658637320;
        integPoints2D_[2][154] = 0.003515625000000;
        integPoints2D_[2][155] = 0.001967799696013;
        integPoints2D_[2][156] = 0.001967799696013;
        integPoints2D_[2][157] = 0.001967799696013;
        integPoints2D_[2][158] = 0.002068658637320;
        integPoints2D_[2][159] = 0.002068658637320;
        integPoints2D_[2][160] = 0.002068658637320;
        integPoints2D_[2][161] = 0.003515625000000;
        integPoints2D_[2][162] = 0.001967799696013;
        integPoints2D_[2][163] = 0.001967799696013;
        integPoints2D_[2][164] = 0.001967799696013;
        integPoints2D_[2][165] = 0.002068658637320;
        integPoints2D_[2][166] = 0.002068658637320;
        integPoints2D_[2][167] = 0.002068658637320;
        integPoints2D_[2][168] = 0.003515625000000;
        integPoints2D_[2][169] = 0.001967799696013;
        integPoints2D_[2][170] = 0.001967799696013;
        integPoints2D_[2][171] = 0.001967799696013;
        integPoints2D_[2][172] = 0.002068658637320;
        integPoints2D_[2][173] = 0.002068658637320;
        integPoints2D_[2][174] = 0.002068658637320;
        integPoints2D_[2][175] = 0.003515625000000;
        integPoints2D_[2][176] = 0.001967799696013;
        integPoints2D_[2][177] = 0.001967799696013;
        integPoints2D_[2][178] = 0.001967799696013;
        integPoints2D_[2][179] = 0.002068658637320;
        integPoints2D_[2][180] = 0.002068658637320;
        integPoints2D_[2][181] = 0.002068658637320;
        integPoints2D_[2][182] = 0.003515625000000;
        integPoints2D_[2][183] = 0.001967799696013;
        integPoints2D_[2][184] = 0.001967799696013;
        integPoints2D_[2][185] = 0.001967799696013;
        integPoints2D_[2][186] = 0.002068658637320;
        integPoints2D_[2][187] = 0.002068658637320;
        integPoints2D_[2][188] = 0.002068658637320;
        integPoints2D_[2][189] = 0.003515625000000;
        integPoints2D_[2][190] = 0.001967799696013;
        integPoints2D_[2][191] = 0.001967799696013;
        integPoints2D_[2][192] = 0.001967799696013;
        integPoints2D_[2][193] = 0.002068658637320;
        integPoints2D_[2][194] = 0.002068658637320;
        integPoints2D_[2][195] = 0.002068658637320;
        integPoints2D_[2][196] = 0.003515625000000;
        integPoints2D_[2][197] = 0.001967799696013;
        integPoints2D_[2][198] = 0.001967799696013;
        integPoints2D_[2][199] = 0.001967799696013;
        integPoints2D_[2][200] = 0.002068658637320;
        integPoints2D_[2][201] = 0.002068658637320;
        integPoints2D_[2][202] = 0.002068658637320;
        integPoints2D_[2][203] = 0.003515625000000;
        integPoints2D_[2][204] = 0.001967799696013;
        integPoints2D_[2][205] = 0.001967799696013;
        integPoints2D_[2][206] = 0.001967799696013;
        integPoints2D_[2][207] = 0.002068658637320;
        integPoints2D_[2][208] = 0.002068658637320;
        integPoints2D_[2][209] = 0.002068658637320;
        integPoints2D_[2][210] = 0.003515625000000;
        integPoints2D_[2][211] = 0.001967799696013;
        integPoints2D_[2][212] = 0.001967799696013;
        integPoints2D_[2][213] = 0.001967799696013;
        integPoints2D_[2][214] = 0.002068658637320;
        integPoints2D_[2][215] = 0.002068658637320;
        integPoints2D_[2][216] = 0.002068658637320;
        integPoints2D_[2][217] = 0.003515625000000;
        integPoints2D_[2][218] = 0.001967799696013;
        integPoints2D_[2][219] = 0.001967799696013;
        integPoints2D_[2][220] = 0.001967799696013;
        integPoints2D_[2][221] = 0.002068658637320;
        integPoints2D_[2][222] = 0.002068658637320;
        integPoints2D_[2][223] = 0.002068658637320;
    }
    else
    {
        std::cout << "SELECT A VALID NUMBER FOR THE AMOUNT OF HAMMER POINTS. THE SELECTED NUMBER IS" << nPoints2D_ << "." << std::endl;
    }
}