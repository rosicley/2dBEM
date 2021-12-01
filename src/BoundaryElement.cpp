#include "BoundaryElement.h"

BoundaryElement::BoundaryElement(const int &index, const std::vector<Node *> &geoConnection, const std::vector<CollocationPoint *> &collocConnection, IntegQuadrature *quadrature)
{
    index_ = index;
    geoConnection_ = geoConnection;
    collocConnection_ = collocConnection;
    quadrature_ = quadrature;
    order_ = geoConnection_.size() - 1;

    xsis_.resize(order_ + 1);
    xsis_[0] = -1.0;
    xsis_[1] = 1.0;
    double aux = 2.0 / static_cast<double>(order_);
    for (int i = 0; i < order_ - 1; i++)
    {
        xsis_[2 + i] = xsis_[0] + static_cast<double>(i + 1) * aux;
    }
}

BoundaryElement::BoundaryElement(const int &order)
{
    order_ = order;
    xsis_.resize(order_ + 1);
    xsis_[0] = -1.0;
    xsis_[1] = 1.0;
    double aux = 2.0 / static_cast<double>(order_);
    for (int i = 0; i < order_ - 1; i++)
    {
        xsis_[2 + i] = xsis_[0] + static_cast<double>(i + 1) * aux;
    }
}

BoundaryElement::~BoundaryElement()
{
}

std::vector<Node *> BoundaryElement::getConnection()
{
    return geoConnection_;
}

Node *BoundaryElement::getNode(const int &index)
{
    return geoConnection_[index];
}

int BoundaryElement::getIndex()
{
    return index_;
}

std::vector<CollocationPoint *> BoundaryElement::getCollocationConnection()
{
    return collocConnection_;
}

void BoundaryElement::interpolateGeometricalCoordinate(const double &xsi, vecDouble &coord)
{
    vecDouble phi(geoConnection_.size());
    getShapeFunction(xsi, phi);
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;

    int cont = 0;
    for (Node *n : geoConnection_)
    {
        coord[0] += n->getCoordinate(0) * phi[cont];
        coord[1] += n->getCoordinate(1) * phi[cont];
        coord[2] += n->getCoordinate(2) * phi[cont];
        cont++;
    }
}

void BoundaryElement::coordOfSourcePoints(const double &xsi, const double &offset, vecDouble &coord)
{
    int nNodes = geoConnection_.size();
    vecDouble phi(nNodes), dphi_dxsi(nNodes);
    getShapeFunction(xsi, phi);
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;

    getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

    vecDouble tangent(2), norm(2), coordNode(3);
    for (int in = 0; in < nNodes; in++)
    {
        geoConnection_[in]->getCoordinate(coordNode);
        coord[0] += phi[in] * coordNode[0];
        coord[1] += phi[in] * coordNode[1];
        tangent[0] += dphi_dxsi[in] * coordNode[0];
        tangent[1] += dphi_dxsi[in] * coordNode[1];
    }
    double jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
    norm[0] = tangent[1] / jac;
    norm[1] = -tangent[0] / jac;

    coord[0] += offset * norm[0];
    coord[1] += offset * norm[1];
}

void BoundaryElement::getShapeFunction(const double &xsi, vecDouble &phi)
{
    setValueVecDouble(phi, 1.0);
    for (int i = 0; i <= order_; i++)
    {
        for (int j = 0; j <= order_; j++)
        {
            if (i != j)
            {
                phi[i] *= (xsi - xsis_[j]) / (xsis_[i] - xsis_[j]);
            }
        }
    }
}

void BoundaryElement::getShapeFunctionAndDerivate(const double &xsi, vecDouble &phi, vecDouble &dphi_dxsi)
{
    setValueVecDouble(phi, 1.0);
    setValueVecDouble(dphi_dxsi, 0.0);
    for (int i = 0; i <= order_; i++)
    {
        for (int j = 0; j <= order_; j++)
        {
            if (i != j)
            {
                double aux = 1.0;
                for (int k = 0; k <= order_; k++)
                {
                    if (i != k and j != k)
                    {
                        aux *= (xsi - xsis_[k]) / (xsis_[i] - xsis_[k]);
                    }
                }
                phi[i] *= (xsi - xsis_[j]) / (xsis_[i] - xsis_[j]);
                dphi_dxsi[i] += (1.0 / (xsis_[i] - xsis_[j]) * aux);
            }
        }
    }
}

void BoundaryElement::crossSectionProperties(double &area, double &perim, double &Ix, double &Iy, double &Mex, double &Mey)
{
    const int nNodes = geoConnection_.size();
    double xsi, weight, jac, JACW;
    vecDouble phi(nNodes), dphi_dxsi(nNodes), initialCoord(3);

    for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
    {
        quadrature_->get1DIntegrationPoint(ip, xsi, weight);
        getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

        vecDouble tangent(2, 0.0), norm(2, 0.0), coord(2, 0.0);
        for (int in = 0; in < nNodes; in++)
        {
            geoConnection_[in]->getCoordinate(initialCoord);
            coord[0] += phi[in] * initialCoord[0];
            coord[1] += phi[in] * initialCoord[1];
            tangent[0] += dphi_dxsi[in] * initialCoord[0];
            tangent[1] += dphi_dxsi[in] * initialCoord[1];
        }
        jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
        norm[0] = tangent[1] / jac;
        norm[1] = -tangent[0] / jac;

        JACW = jac * weight;

        area += (coord[0] * 0.5 * norm[0] + coord[1] * 0.5 * norm[1]) * JACW;
        perim += JACW;
        Iy += (pow(coord[0], 3) / 3.0 * norm[0]) * JACW;
        Ix += (pow(coord[1], 3) / 3.0 * norm[1]) * JACW;
        Mex += (pow(coord[0], 2) / 2.0 * norm[0]) * JACW;
        Mey += (pow(coord[1], 2) / 2.0 * norm[1]) * JACW;
    }
}

void BoundaryElement::potentialContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matGlocal, matDouble &matHlocal)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matGlocal.size() != nSource)
    {
        setSizeMatDouble(matGlocal, nSource, nNodes);
        setSizeMatDouble(matHlocal, nSource, nNodes);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, Rx, Ry, R, Q, U, Usingular, VPC, jacSP, pi = 3.141592653589793;
    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3), phiSP(nNodes), dphiSP_dxsi(nNodes);

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            Rx = coordIP[0] - coordSourceP[0];
            Ry = coordIP[1] - coordSourceP[1];
            R = sqrt(Rx * Rx + Ry * Ry);

            Q = -1.0 / (2.0 * pi * R) * (Rx / R * normal[0] + Ry / R * normal[1]) * JACW;

            if (!sourcePoints[isp]->checkElement(index_)) ///PONTO FONTE NÃO PERTENCE AO ELEMENTO
            {
                U = -1.0 / (2.0 * pi) * log(R) * JACW;
                Usingular = 0.0;
                VPC = 0.0;
            }
            else
            {
                int indexSP = sourcePoints[isp]->getIndex();
                double xsiSP;
                for (int in = 0; in < nNodes; in++)
                {
                    if (indexSP == collocConnection_[in]->getIndex())
                    {
                        xsiSP = xsis_[in];
                        break;
                    }
                }
                double eps = xsi - xsiSP;
                double Rast = jac * eps;
                getShapeFunctionAndDerivate(xsiSP, phiSP, dphiSP_dxsi);
                vecDouble tangentSP(2, 0.0);
                for (int in = 0; in < nNodes; in++)
                {
                    geoConnection_[in]->getCoordinate(coordNode);
                    tangentSP[0] += dphiSP_dxsi[in] * coordNode[0];
                    tangentSP[1] += dphiSP_dxsi[in] * coordNode[1];
                }
                jacSP = sqrt(tangentSP[0] * tangentSP[0] + tangentSP[1] * tangentSP[1]);
                U = -1.0 / (2.0 * pi) * log(R) * JACW;
                Usingular = -1.0 / (2.0 * pi) * log(fabs(Rast)) * jacSP * weight;

                if (ip == 0)
                {
                    if (fabs(xsiSP - 1.0) <= 1.0e-12 or fabs(xsiSP + 1.0) <= 1.0e-12)
                    {
                        VPC = -1.0 / (2.0 * pi) * jacSP * (2.0 * log(2.0 * jacSP) - 2.0);
                    }
                    else
                    {
                        VPC = -1.0 / (2.0 * pi) * jacSP * ((1.0 + xsiSP) * log(jacSP * (1.0 + xsiSP)) + (1.0 - xsiSP) * log(jacSP * (1.0 - xsiSP)) - (1.0 + xsiSP) - (1.0 - xsiSP));
                    }
                    for (int in = 0; in < nNodes; in++)
                    {
                        matGlocal[isp][in] += VPC * phiSP[in];
                    }
                }
            }
            for (int in = 0; in < nNodes; in++)
            {
                matHlocal[isp][in] += Q * phi[in];
                matGlocal[isp][in] += U * phi[in] + (-Usingular) * phiSP[in];
            }
        }
    }
}

void BoundaryElement::calculateInternalFlux(const std::vector<SourcePoint *> &sourcePoints, matDouble &matGlocal, matDouble &matHlocal)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matGlocal.size() != 2 * nSource)
    {
        setSizeMatDouble(matGlocal, 2 * nSource, nNodes);
        setSizeMatDouble(matHlocal, 2 * nSource, nNodes);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, Rx, Ry, R, S1, S2, D1, D2, pi = 3.141592653589793;
    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            Rx = coordIP[0] - coordSourceP[0];
            Ry = coordIP[1] - coordSourceP[1];
            R = sqrt(Rx * Rx + Ry * Ry);

            double dR_dx = Rx / R;
            double dR_dy = Ry / R;
            double DRDN = dR_dx * normal[0] + dR_dy * normal[1];

            S1 = 1.0 / (2.0 * pi * R * R) * (normal[0] - 2.0 * dR_dx * DRDN) * JACW;
            S2 = 1.0 / (2.0 * pi * R * R) * (normal[1] - 2.0 * dR_dy * DRDN) * JACW;

            D1 = 1.0 / (2.0 * pi * R) * dR_dx * JACW;
            D2 = 1.0 / (2.0 * pi * R) * dR_dy * JACW;

            for (int in = 0; in < nNodes; in++)
            {
                matHlocal[2 * isp][in] += S1 * phi[in];
                matHlocal[2 * isp + 1][in] += S2 * phi[in];

                matGlocal[2 * isp][in] += D1 * phi[in];
                matGlocal[2 * isp + 1][in] += D2 * phi[in];
            }
        }
    }
}

void BoundaryElement::elasticityContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matGlocal, matDouble &matHlocal, Material *material)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matGlocal.size() == 0)
    {
        setSizeMatDouble(matGlocal, 2 * nSource, 2 * nNodes);
        setSizeMatDouble(matHlocal, 2 * nSource, 2 * nNodes);
    }
    else if (matGlocal.size() != 2 * nSource or matGlocal[0].size() != 2 * nNodes)
    {
        setSizeMatDouble(matGlocal, 2 * nSource, 2 * nNodes);
        setSizeMatDouble(matHlocal, 2 * nSource, 2 * nNodes);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, normRadius, pi = 3.141592653589793, poisson, young, density, auxU, auxU2, auxP, G, dRadius_dNormal;
    double radius[2], dRadius_dXi[2];
    double Pfund[2][2], Ufund[2][2];
    const double identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));
    auxU = 1.0 / (8.0 * pi * G * (1.0 - poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            // vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            double tangent[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, coordIP[2] = {0.0, 0.0};
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            auxU2 = -(3.0 - 4.0 * poisson) * log(normRadius);
            auxP = -1.0 / (4.0 * pi * (1.0 - poisson) * normRadius);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;
            dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];

            if (!sourcePoints[isp]->checkElement(index_)) ///PONTO FONTE NÃO PERTENCE AO ELEMENTO
            {
                for (int l = 0; l < 2; l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Ufund[l][k] = auxU * (auxU2 * identity[l][k] - (7.0 - 8.0 * poisson) * 0.5 * identity[l][k] + dRadius_dXi[l] * dRadius_dXi[k]) * JACW;
                        Pfund[l][k] = auxP * (dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[l][k] + 2.0 * dRadius_dXi[l] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (normal[l] * dRadius_dXi[k] - normal[k] * dRadius_dXi[l])) * JACW;
                    }
                }
                for (int l = 0; l < 2; l++)
                {
                    for (int in = 0; in < nNodes; in++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            matHlocal[2 * isp + l][2 * in + k] += Pfund[l][k] * phi[in];
                            matGlocal[2 * isp + l][2 * in + k] += Ufund[l][k] * phi[in];
                        }
                    }
                }
            }
            else
            {
                double Usingular[2][2], Psingular[2][2], dRadiusAst_dXi[2], VPC_U[2][2] = {{0.0, 0.0}, {0.0, 0.0}}, VPC_P[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
                double xsiSP;
                for (int in = 0, indexSP = sourcePoints[isp]->getIndex(); in < nNodes; in++)
                {
                    if (indexSP == collocConnection_[in]->getIndex())
                    {
                        xsiSP = xsis_[in];
                        break;
                    }
                }

                vecDouble phiSP(nNodes), dphiSP_dxsi(nNodes);
                getShapeFunctionAndDerivate(xsiSP, phiSP, dphiSP_dxsi);

                double tangentSP[2] = {0.0, 0.0};
                for (int in = 0; in < nNodes; in++)
                {
                    geoConnection_[in]->getCoordinate(coordNode);
                    tangentSP[0] += dphiSP_dxsi[in] * coordNode[0];
                    tangentSP[1] += dphiSP_dxsi[in] * coordNode[1];
                }

                double eps = xsi - xsiSP;
                double radiusAst = jac * eps;

                double jacSP = sqrt(tangentSP[0] * tangentSP[0] + tangentSP[1] * tangentSP[1]);
                double normalSP[2] = {tangentSP[1] / jacSP, -tangentSP[0] / jacSP};
                double jacSPW = jacSP * weight;

                dRadiusAst_dXi[0] = tangentSP[0] / jacSP;
                dRadiusAst_dXi[1] = tangentSP[1] / jacSP;

                for (int l = 0; l < 2; l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Ufund[l][k] = auxU * (auxU2 * identity[l][k] - (7.0 - 8.0 * poisson) * 0.5 * identity[l][k] + dRadius_dXi[l] * dRadius_dXi[k]) * JACW;
                        Pfund[l][k] = auxP * (dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[l][k] + 2.0 * dRadius_dXi[l] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (normal[l] * dRadius_dXi[k] - normal[k] * dRadius_dXi[l])) * JACW;

                        Usingular[l][k] = auxU * (-(3.0 - 4.0 * poisson) * log(fabs(radiusAst)) * identity[l][k]) * jacSPW;
                        Psingular[l][k] = -1.0 / (4.0 * pi * (1.0 - poisson) * radiusAst) * ((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) * jacSPW;
                    }
                }

                if (ip == 0)
                {
                    if (fabs(xsiSP - 1.0) <= 1.0e-12 or fabs(xsiSP + 1.0) <= 1.0e-12)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                VPC_U[l][k] = -auxU * (3.0 - 4.0 * poisson) * identity[l][k] * jacSP * (2.0 * log(2.0 * jacSP) - 2.0);

                                VPC_P[l][k] = ((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) / (4.0 * pi * (1.0 - poisson)) * log(2.0) * xsiSP;
                            }
                        }
                    }
                    else
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                VPC_U[l][k] = -auxU * (3.0 - 4.0 * poisson) * identity[l][k] * jacSP * ((1.0 + xsiSP) * log(jacSP * (1.0 + xsiSP)) + (1.0 - xsiSP) * log(jacSP * (1.0 - xsiSP)) - 2.0);

                                VPC_P[l][k] = -((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) / (4.0 * pi * (1.0 - poisson)) * (log(1.0 - xsiSP) - log(1.0 + xsiSP));
                            }
                        }
                    }
                }

                for (int l = 0; l < 2; l++)
                {
                    for (int in = 0; in < nNodes; in++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            matHlocal[2 * isp + l][2 * in + k] += Pfund[l][k] * phi[in] + (-Psingular[l][k] + VPC_P[l][k]) * phiSP[in];
                            matGlocal[2 * isp + l][2 * in + k] += Ufund[l][k] * phi[in] + (-Usingular[l][k] + VPC_U[l][k]) * phiSP[in];
                        }
                    }
                }
            }
        }
    }
}

void BoundaryElement::calculateInternalStress(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matG.size() == 0)
    {
        setSizeMatDouble(matG, 4 * nSource, 2 * nNodes);
        setSizeMatDouble(matH, 4 * nSource, 2 * nNodes);
    }
    else if (matG.size() != 4 * nSource or matG[0].size() != 2 * nNodes)
    {
        setSizeMatDouble(matG, 4 * nSource, 2 * nNodes);
        setSizeMatDouble(matH, 4 * nSource, 2 * nNodes);
    }
    setValueMatDouble(matG, 0.0);
    setValueMatDouble(matH, 0.0);

    double xsi, weight, jac, JACW, normRadius, pi = 3.141592653589793, poisson, young, density, auxU, auxU2, auxP, G, dRadius_dNormal;
    double radius[2], dRadius_dXi[2];
    double D[2][2][2], S[2][2][2], DAUX[4][2], SAUX[4][2];
    const double identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));
    auxU = 1.0 / (8.0 * pi * G * (1.0 - poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            double tangent[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, coordIP[2] = {0.0, 0.0};
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            auxU2 = -(3.0 - 4.0 * poisson) * log(normRadius);
            auxP = -1.0 / (4.0 * pi * (1.0 - poisson) * normRadius);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;
            dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];

            for (int k = 0; k < 2; k++)
            {
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        D[k][i][j] = -auxP * ((1.0 - 2.0 * poisson) * (identity[k][i] * dRadius_dXi[j] + identity[k][j] * dRadius_dXi[i] - identity[i][j] * dRadius_dXi[k]) + 2.0 * dRadius_dXi[i] * dRadius_dXi[j] * dRadius_dXi[k]) * JACW;

                        S[k][i][j] = G / (2.0 * pi * (1.0 - poisson) * normRadius * normRadius) * (2.0 * dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[i][j] * dRadius_dXi[k] + poisson * (identity[i][k] * dRadius_dXi[j] + identity[j][k] * dRadius_dXi[i]) - 4.0 * dRadius_dXi[i] * dRadius_dXi[j] * dRadius_dXi[k]) + 2.0 * poisson * (normal[i] * dRadius_dXi[j] * dRadius_dXi[k] + normal[j] * dRadius_dXi[i] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (2.0 * normal[k] * dRadius_dXi[i] * dRadius_dXi[j] + normal[j] * identity[i][k] + normal[i] * identity[j][k]) - (1.0 - 4.0 * poisson) * normal[k] * identity[i][j]) * JACW;
                    }
                }
            }
            DAUX[0][0] = D[0][0][0];
            DAUX[0][1] = D[1][0][0];
            DAUX[1][0] = D[0][0][1];
            DAUX[1][1] = D[1][0][1];
            DAUX[2][0] = D[0][1][0];
            DAUX[2][1] = D[1][1][0];
            DAUX[3][0] = D[0][1][1];
            DAUX[3][1] = D[1][1][1];
            SAUX[0][0] = S[0][0][0];
            SAUX[0][1] = S[1][0][0];
            SAUX[1][0] = S[0][0][1];
            SAUX[1][1] = S[1][0][1];
            SAUX[2][0] = S[0][1][0];
            SAUX[2][1] = S[1][1][0];
            SAUX[3][0] = S[0][1][1];
            SAUX[3][1] = S[1][1][1];
            for (int l = 0; l < 4; l++)
            {
                for (int in = 0; in < nNodes; in++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        matH[4 * isp + l][2 * in + k] += SAUX[l][k] * phi[in];
                        matG[4 * isp + l][2 * in + k] += DAUX[l][k] * phi[in];
                    }
                }
            }
        }
    }
}

void BoundaryElement::elasticityBodyForceContribution(const std::vector<SourcePoint *> &sourcePoints, vecDouble &vecB, Material *material, const vecDouble &force)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (vecB.size() != 2 * nSource)
    {
        vecB.resize(2 * nSource);
    }
    setValueVecDouble(vecB, 0.0);

    double xsi, weight, jac, JACW, normRadius, pi = 3.141592653589793, poisson, young, density, G;
    double radius[2], dRadius_dXi[2];

    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            double tangent[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, coordIP[2] = {0.0, 0.0};
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;

            double dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];
            double dRadius_dB = dRadius_dXi[0] * force[0] + dRadius_dXi[1] * force[1];

            vecB[2 * isp] += (-normRadius / (8.0 * pi * G) * (2.0 * log(normRadius) + 1.0) *
                              (dRadius_dNormal * force[0] - 1.0 / (2.0 * (1.0 - poisson)) * dRadius_dB * normal[0])) *
                             JACW;

            vecB[2 * isp + 1] += (-normRadius / (8.0 * pi * G) * (2.0 * log(normRadius) + 1.0) *
                                  (dRadius_dNormal * force[1] - 1.0 / (2.0 * (1.0 - poisson)) * dRadius_dB * normal[1])) *
                                 JACW;
        }
    }
}

void BoundaryElement::calculateInternalStressBodyForce(const std::vector<SourcePoint *> &sourcePoints, vecDouble &vecB, Material *material, const vecDouble &force)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (vecB.size() != 4 * nSource)
    {
        vecB.resize(4 * nSource);
    }
    setValueVecDouble(vecB, 0.0);

    double xsi, weight, jac, JACW, normRadius, pi = 3.141592653589793, poisson, young, density, G;
    double radius[2], dRadius_dXi[2], D[2][2];
    const double identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);

            double tangent[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, coordIP[2] = {0.0, 0.0};
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;

            double dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];
            double dRadius_dB = dRadius_dXi[0] * force[0] + dRadius_dXi[1] * force[1];
            double bm_nm = force[0] * normal[0] + force[1] * normal[1];

            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    double aux0 = 2.0 * dRadius_dNormal * (force[i] * dRadius_dXi[j] + force[j] * dRadius_dXi[i]);
                    double aux1 = poisson * identity[i][j] * (2.0 * dRadius_dNormal * dRadius_dB + (1.0 + 2.0 * log(normRadius)) * bm_nm);
                    double aux2 = -dRadius_dB * (normal[i] * dRadius_dXi[j] + normal[j] * dRadius_dXi[i]);
                    double aux3 = (1.0 - 2.0 * poisson) * 0.5 * (1.0 + 2.0 * log(normRadius)) * (force[i] * normal[j] + force[j] * normal[i]);
                    D[i][j] = 1.0 / (8.0 * pi) * (aux0 + 1.0 / (1.0 - poisson) * (aux1 + aux2 + aux3)) * JACW;
                }
            }

            vecB[4 * isp] += D[0][0];
            vecB[4 * isp + 1] += D[0][1];
            vecB[4 * isp + 2] += D[1][0];
            vecB[4 * isp + 3] += D[1][1];
        }
    }
}

/////////////////DISCONT

DiscontBoundaryElement::DiscontBoundaryElement(const int &index, const std::vector<Node *> &geoConnection, const std::vector<CollocationPoint *> &collocConnection, IntegQuadrature *quadrature, const std::string &discont, const double &collocParam) : BoundaryElement(index, geoConnection, collocConnection, quadrature)
{
    discont_ = discont;
    collocParam_ = collocParam;
    collocXsis_ = xsis_;
    double aux = 2.0 / static_cast<double>(order_) * collocParam_;
    if (discont_ == "both")
    {
        collocXsis_[0] += aux;
        collocXsis_[1] -= aux;
    }
    else if (discont_ == "left")
    {
        collocXsis_[0] += aux;
    }
    else if (discont_ == "right")
    {
        collocXsis_[1] -= aux;
    }
    else
    {
        std::cout << "ERROR..\n";
    }
}

DiscontBoundaryElement::DiscontBoundaryElement(const int &order, const std::string &discont, const double &param) : BoundaryElement(order)
{
    discont_ = discont;
    collocParam_ = param;
    collocXsis_ = xsis_;
    double aux = 2.0 / static_cast<double>(order_) * collocParam_;
    if (discont_ == "both")
    {
        collocXsis_[0] += aux;
        collocXsis_[1] -= aux;
    }
    else if (discont_ == "left")
    {
        collocXsis_[0] += aux;
    }
    else if (discont_ == "right")
    {
        collocXsis_[1] -= aux;
    }
    else
    {
        std::cout << "ERROR..\n";
    }
}

DiscontBoundaryElement::~DiscontBoundaryElement()
{
}

void DiscontBoundaryElement::getCollocationShapeFunction(const double &xsi, vecDouble &phi)
{
    setValueVecDouble(phi, 1.0);
    for (int i = 0; i <= order_; i++)
    {
        for (int j = 0; j <= order_; j++)
        {
            if (i != j)
            {
                phi[i] *= (xsi - collocXsis_[j]) / (collocXsis_[i] - collocXsis_[j]);
            }
        }
    }
}

void DiscontBoundaryElement::potentialContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matGlocal, matDouble &matHlocal)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matGlocal.size() != nSource)
    {
        setSizeMatDouble(matGlocal, nSource, nNodes);
        setSizeMatDouble(matHlocal, nSource, nNodes);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, Rx, Ry, R, Q, U, Usingular, VPC, jacSP, pi = 3.141592653589793;
    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3), phiColloc(nNodes), phiSP(nNodes), dphiSP_dxsi(nNodes), phiCollocSP(nNodes);

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);
            getCollocationShapeFunction(xsi, phiColloc);

            vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            Rx = coordIP[0] - coordSourceP[0];
            Ry = coordIP[1] - coordSourceP[1];
            R = sqrt(Rx * Rx + Ry * Ry);

            Q = -1.0 / (2.0 * pi * R) * (Rx / R * normal[0] + Ry / R * normal[1]) * JACW;

            if (!sourcePoints[isp]->checkElement(index_)) ///PONTO FONTE NÃO PERTENCE AO ELEMENTO
            {
                U = -1.0 / (2.0 * pi) * log(R) * JACW;
                Usingular = 0.0;
                VPC = 0.0;
            }
            else
            {
                int indexSP = sourcePoints[isp]->getIndex();
                double xsiSP;
                for (int in = 0; in < nNodes; in++)
                {
                    if (indexSP == collocConnection_[in]->getIndex())
                    {
                        xsiSP = collocXsis_[in];
                        break;
                    }
                }
                double eps = xsi - xsiSP;
                double Rast = jac * eps;
                getShapeFunctionAndDerivate(xsiSP, phiSP, dphiSP_dxsi);
                getCollocationShapeFunction(xsiSP, phiCollocSP);
                vecDouble tangentSP(2, 0.0);
                for (int in = 0; in < nNodes; in++)
                {
                    geoConnection_[in]->getCoordinate(coordNode);
                    tangentSP[0] += dphiSP_dxsi[in] * coordNode[0];
                    tangentSP[1] += dphiSP_dxsi[in] * coordNode[1];
                }
                jacSP = sqrt(tangentSP[0] * tangentSP[0] + tangentSP[1] * tangentSP[1]);
                U = -1.0 / (2.0 * pi) * log(R) * JACW;
                Usingular = -1.0 / (2.0 * pi) * log(fabs(Rast)) * jacSP * weight;

                if (ip == 0)
                {
                    if (fabs(xsiSP - 1.0) <= 1.0e-12 or fabs(xsiSP + 1.0) <= 1.0e-12)
                    {
                        VPC = -1.0 / (2.0 * pi) * jacSP * (2.0 * log(2.0 * jacSP) - 2.0);
                    }
                    else
                    {
                        VPC = -1.0 / (2.0 * pi) * jacSP * ((1.0 + xsiSP) * log(jacSP * (1.0 + xsiSP)) + (1.0 - xsiSP) * log(jacSP * (1.0 - xsiSP)) - (1.0 + xsiSP) - (1.0 - xsiSP));
                    }
                    for (int in = 0; in < nNodes; in++)
                    {
                        matGlocal[isp][in] += VPC * phiCollocSP[in];
                    }
                }
            }
            for (int in = 0; in < nNodes; in++)
            {
                matHlocal[isp][in] += Q * phiColloc[in];
                matGlocal[isp][in] += U * phiColloc[in] + (-Usingular) * phiCollocSP[in];
            }
        }
    }
}

void DiscontBoundaryElement::calculateInternalFlux(const std::vector<SourcePoint *> &internalPoints, matDouble &matGlocal, matDouble &matHlocal)
{
    const int nNodes = geoConnection_.size();
    const int nInternalP = internalPoints.size();

    if (matGlocal.size() != 2 * nInternalP)
    {
        setSizeMatDouble(matGlocal, 2 * nInternalP, nNodes, 0.0);
        setSizeMatDouble(matHlocal, 2 * nInternalP, nNodes, 0.0);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, Rx, Ry, R, S1, S2, D1, D2, pi = 3.141592653589793;
    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3), phiColloc(nNodes);

    for (int isp = 0; isp < nInternalP; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);
            getCollocationShapeFunction(xsi, phiColloc);

            vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            internalPoints[isp]->getCoordinate(coordSourceP);

            Rx = coordIP[0] - coordSourceP[0];
            Ry = coordIP[1] - coordSourceP[1];
            R = sqrt(Rx * Rx + Ry * Ry);

            double dR_dx = Rx / R;
            double dR_dy = Ry / R;
            double DRDN = dR_dx * normal[0] + dR_dy * normal[1];

            S1 = 1.0 / (2.0 * pi * R * R) * (normal[0] - 2.0 * dR_dx * DRDN) * JACW;
            S2 = 1.0 / (2.0 * pi * R * R) * (normal[1] - 2.0 * dR_dy * DRDN) * JACW;

            D1 = 1.0 / (2.0 * pi * R) * dR_dx * JACW;
            D2 = 1.0 / (2.0 * pi * R) * dR_dy * JACW;

            for (int in = 0; in < nNodes; in++)
            {
                matHlocal[2 * isp][in] += S1 * phiColloc[in];
                matHlocal[2 * isp + 1][in] += S2 * phiColloc[in];

                matGlocal[2 * isp][in] += D1 * phiColloc[in];
                matGlocal[2 * isp + 1][in] += D2 * phiColloc[in];
            }
        }
    }
}

void DiscontBoundaryElement::elasticityContribution(const std::vector<SourcePoint *> &sourcePoints, matDouble &matGlocal, matDouble &matHlocal, Material *material)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matGlocal.size() == 0)
    {
        setSizeMatDouble(matGlocal, 2 * nSource, 2 * nNodes);
        setSizeMatDouble(matHlocal, 2 * nSource, 2 * nNodes);
    }
    else if (matGlocal.size() != 2 * nSource or matGlocal[0].size() != 2 * nNodes)
    {
        setSizeMatDouble(matGlocal, 2 * nSource, 2 * nNodes);
        setSizeMatDouble(matHlocal, 2 * nSource, 2 * nNodes);
    }
    setValueMatDouble(matGlocal, 0.0);
    setValueMatDouble(matHlocal, 0.0);

    double xsi, weight, jac, JACW, normRadius, jacSP, pi = 3.141592653589793, poisson, young, density, auxU, auxU2, auxP, G, dRadius_dNormal;
    double radius[2], dRadius_dXi[2];
    double Pfund[2][2], Ufund[2][2], identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    vecDouble phi(nNodes), phiColloc(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));
    auxU = 1.0 / (8.0 * pi * G * (1.0 - poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);
            getCollocationShapeFunction(xsi, phiColloc);

            vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            auxU2 = -(3.0 - 4.0 * poisson) * log(normRadius);
            auxP = -1.0 / (4.0 * pi * (1.0 - poisson) * normRadius);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;
            dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];

            if (!sourcePoints[isp]->checkElement(index_)) ///PONTO FONTE NÃO PERTENCE AO ELEMENTO
            {
                for (int l = 0; l < 2; l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Ufund[l][k] = auxU * (auxU2 * identity[l][k] - (7.0 - 8.0 * poisson) * 0.5 * identity[l][k] + dRadius_dXi[l] * dRadius_dXi[k]) * JACW;
                        Pfund[l][k] = auxP * (dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[l][k] + 2.0 * dRadius_dXi[l] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (normal[l] * dRadius_dXi[k] - normal[k] * dRadius_dXi[l])) * JACW;
                    }
                }
                for (int l = 0; l < 2; l++)
                {
                    for (int in = 0; in < nNodes; in++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            matHlocal[2 * isp + l][2 * in + k] += Pfund[l][k] * phiColloc[in];
                            matGlocal[2 * isp + l][2 * in + k] += Ufund[l][k] * phiColloc[in];
                        }
                    }
                }
            }
            else
            {
                double Usingular[2][2], Psingular[2][2], dRadiusAst_dXi[2], VPC_U[2][2] = {{0.0, 0.0}, {0.0, 0.0}}, VPC_P[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
                double xsiSP;
                for (int in = 0, indexSP = sourcePoints[isp]->getIndex(); in < nNodes; in++)
                {
                    if (indexSP == collocConnection_[in]->getIndex())
                    {
                        xsiSP = collocXsis_[in];
                        break;
                    }
                }

                vecDouble phiSP(nNodes), dphiSP_dxsi(nNodes), phiCollocSP(nNodes);
                getShapeFunctionAndDerivate(xsiSP, phiSP, dphiSP_dxsi);
                getCollocationShapeFunction(xsiSP, phiCollocSP);

                double tangentSP[2] = {0.0, 0.0};
                for (int in = 0; in < nNodes; in++)
                {
                    geoConnection_[in]->getCoordinate(coordNode);
                    tangentSP[0] += dphiSP_dxsi[in] * coordNode[0];
                    tangentSP[1] += dphiSP_dxsi[in] * coordNode[1];
                }

                double eps = xsi - xsiSP;
                double radiusAst = jac * eps;

                double jacSP = sqrt(tangentSP[0] * tangentSP[0] + tangentSP[1] * tangentSP[1]);
                double normalSP[2] = {tangentSP[1] / jacSP, -tangentSP[0] / jacSP};
                double jacSPW = jacSP * weight;

                dRadiusAst_dXi[0] = tangentSP[0] / jacSP;
                dRadiusAst_dXi[1] = tangentSP[1] / jacSP;

                for (int l = 0; l < 2; l++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        Ufund[l][k] = auxU * (auxU2 * identity[l][k] - (7.0 - 8.0 * poisson) * 0.5 * identity[l][k] + dRadius_dXi[l] * dRadius_dXi[k]) * JACW;
                        Pfund[l][k] = auxP * (dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[l][k] + 2.0 * dRadius_dXi[l] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (normal[l] * dRadius_dXi[k] - normal[k] * dRadius_dXi[l])) * JACW;

                        Usingular[l][k] = auxU * (-(3.0 - 4.0 * poisson) * log(fabs(radiusAst)) * identity[l][k]) * jacSPW;
                        Psingular[l][k] = -1.0 / (4.0 * pi * (1.0 - poisson) * radiusAst) * ((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) * jacSPW;
                    }
                }

                if (ip == 0)
                {
                    if (fabs(xsiSP - 1.0) <= 1.0e-12 or fabs(xsiSP + 1.0) <= 1.0e-12)
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                VPC_U[l][k] = -auxU * (3.0 - 4.0 * poisson) * identity[l][k] * jacSP * (2.0 * log(2.0 * jacSP) - 2.0);

                                VPC_P[l][k] = ((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) / (4.0 * pi * (1.0 - poisson)) * log(2.0) * xsiSP;
                            }
                        }
                    }
                    else
                    {
                        for (int l = 0; l < 2; l++)
                        {
                            for (int k = 0; k < 2; k++)
                            {
                                VPC_U[l][k] = -auxU * (3.0 - 4.0 * poisson) * identity[l][k] * jacSP * ((1.0 + xsiSP) * log(jacSP * (1.0 + xsiSP)) + (1.0 - xsiSP) * log(jacSP * (1.0 - xsiSP)) - 2.0);

                                VPC_P[l][k] = -((1.0 - 2.0 * poisson) * (normalSP[l] * dRadiusAst_dXi[k] - normalSP[k] * dRadiusAst_dXi[l])) / (4.0 * pi * (1.0 - poisson)) * (log(1.0 - xsiSP) - log(1.0 + xsiSP));
                            }
                        }
                    }
                }

                for (int l = 0; l < 2; l++)
                {
                    for (int in = 0; in < nNodes; in++)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            matHlocal[2 * isp + l][2 * in + k] += Pfund[l][k] * phiColloc[in] + (-Psingular[l][k] + VPC_P[l][k]) * phiCollocSP[in];
                            matGlocal[2 * isp + l][2 * in + k] += Ufund[l][k] * phiColloc[in] + (-Usingular[l][k] + VPC_U[l][k]) * phiCollocSP[in];
                        }
                    }
                }
            }
        }
    }
}

void DiscontBoundaryElement::calculateInternalStress(const std::vector<SourcePoint *> &sourcePoints, matDouble &matG, matDouble &matH, Material *material)
{
    const int nNodes = geoConnection_.size();
    const int nSource = sourcePoints.size();

    if (matG.size() == 0)
    {
        setSizeMatDouble(matG, 4 * nSource, 2 * nNodes);
        setSizeMatDouble(matH, 4 * nSource, 2 * nNodes);
    }
    else if (matG.size() != 4 * nSource or matG[0].size() != 2 * nNodes)
    {
        setSizeMatDouble(matG, 4 * nSource, 2 * nNodes);
        setSizeMatDouble(matH, 4 * nSource, 2 * nNodes);
    }
    setValueMatDouble(matG, 0.0);
    setValueMatDouble(matH, 0.0);

    double xsi, weight, jac, JACW, normRadius, pi = 3.141592653589793, poisson, young, density, auxU, auxU2, auxP, G, dRadius_dNormal;
    double radius[2], dRadius_dXi[2];
    double D[2][2][2], S[2][2][2], DAUX[4][2], SAUX[4][2];
    const double identity[2][2] = {{1.0, 0.0}, {0.0, 1.0}};

    vecDouble phi(nNodes), dphi_dxsi(nNodes), coordNode(3), coordSourceP(3), phiColloc(nNodes);

    material->getProperties(young, poisson, density);
    G = young / (2.0 * (1.0 + poisson));
    auxU = 1.0 / (8.0 * pi * G * (1.0 - poisson));

    for (int isp = 0; isp < nSource; isp++)
    {
        for (int ip = 0, nip = quadrature_->get1DPointsNumber(); ip < nip; ip++)
        {
            quadrature_->get1DIntegrationPoint(ip, xsi, weight);
            getShapeFunctionAndDerivate(xsi, phi, dphi_dxsi);
            getCollocationShapeFunction(xsi, phiColloc);

            // vecDouble tangent(2, 0.0), normal(2, 0.0), coordIP(2, 0.0);
            double tangent[2] = {0.0, 0.0}, normal[2] = {0.0, 0.0}, coordIP[2] = {0.0, 0.0};
            for (int in = 0; in < nNodes; in++)
            {
                geoConnection_[in]->getCoordinate(coordNode);
                coordIP[0] += phi[in] * coordNode[0];
                coordIP[1] += phi[in] * coordNode[1];
                tangent[0] += dphi_dxsi[in] * coordNode[0];
                tangent[1] += dphi_dxsi[in] * coordNode[1];
            }
            jac = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1]);
            normal[0] = tangent[1] / jac;
            normal[1] = -tangent[0] / jac;
            JACW = jac * weight;

            sourcePoints[isp]->getCoordinate(coordSourceP);

            radius[0] = coordIP[0] - coordSourceP[0];
            radius[1] = coordIP[1] - coordSourceP[1];
            normRadius = sqrt(radius[0] * radius[0] + radius[1] * radius[1]);

            auxU2 = -(3.0 - 4.0 * poisson) * log(normRadius);
            auxP = -1.0 / (4.0 * pi * (1.0 - poisson) * normRadius);

            dRadius_dXi[0] = radius[0] / normRadius;
            dRadius_dXi[1] = radius[1] / normRadius;
            dRadius_dNormal = dRadius_dXi[0] * normal[0] + dRadius_dXi[1] * normal[1];

            for (int k = 0; k < 2; k++)
            {
                for (int i = 0; i < 2; i++)
                {
                    for (int j = 0; j < 2; j++)
                    {
                        D[k][i][j] = -auxP * ((1.0 - 2.0 * poisson) * (identity[k][i] * dRadius_dXi[j] + identity[k][j] * dRadius_dXi[i] - identity[i][j] * dRadius_dXi[k]) + 2.0 * dRadius_dXi[i] * dRadius_dXi[j] * dRadius_dXi[k]) * JACW;

                        S[k][i][j] = G / (2.0 * pi * (1.0 - poisson) * normRadius * normRadius) * (2.0 * dRadius_dNormal * ((1.0 - 2.0 * poisson) * identity[i][j] * dRadius_dXi[k] + poisson * (identity[i][k] * dRadius_dXi[j] + identity[j][k] * dRadius_dXi[i]) - 4.0 * dRadius_dXi[i] * dRadius_dXi[j] * dRadius_dXi[k]) + 2.0 * poisson * (normal[i] * dRadius_dXi[j] * dRadius_dXi[k] + normal[j] * dRadius_dXi[i] * dRadius_dXi[k]) + (1.0 - 2.0 * poisson) * (2.0 * normal[k] * dRadius_dXi[i] * dRadius_dXi[j] + normal[j] * identity[i][k] + normal[i] * identity[j][k]) - (1.0 - 4.0 * poisson) * normal[k] * identity[i][j]) * JACW;
                    }
                }
            }
            DAUX[0][0] = D[0][0][0];
            DAUX[0][1] = D[1][0][0];
            DAUX[1][0] = D[0][0][1];
            DAUX[1][1] = D[1][0][1];
            DAUX[2][0] = D[0][1][0];
            DAUX[2][1] = D[1][1][0];
            DAUX[3][0] = D[0][1][1];
            DAUX[3][1] = D[1][1][1];
            SAUX[0][0] = S[0][0][0];
            SAUX[0][1] = S[1][0][0];
            SAUX[1][0] = S[0][0][1];
            SAUX[1][1] = S[1][0][1];
            SAUX[2][0] = S[0][1][0];
            SAUX[2][1] = S[1][1][0];
            SAUX[3][0] = S[0][1][1];
            SAUX[3][1] = S[1][1][1];
            for (int l = 0; l < 4; l++)
            {
                for (int in = 0; in < nNodes; in++)
                {
                    for (int k = 0; k < 2; k++)
                    {
                        matH[4 * isp + l][2 * in + k] += SAUX[l][k] * phiColloc[in];
                        matG[4 * isp + l][2 * in + k] += DAUX[l][k] * phiColloc[in];
                    }
                }
            }
        }
    }
}
