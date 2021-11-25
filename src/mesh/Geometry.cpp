#include "Geometry.h"

Geometry::Geometry(const std::string &name)
{
    name_ = name;
}

Geometry::~Geometry() {}

std::string Geometry::getName()
{
    return name_;
}

Point *Geometry::addPoint(std::vector<double> coordinates, const double &lcar, const bool &discretization)
{
    std::string name = "p" + std::to_string(points_.size());
    Point *p = new Point(name, coordinates, lcar, discretization);
    points_[name] = p;
    return p;
}

Line *Geometry::addLine(std::vector<Point *> points, const bool &discretization)
{
    std::string name = "l" + std::to_string(lines_.size());
    Line *l = new Line(name, points, discretization);
    lines_[name] = l;
    return l;
}

Line *Geometry::addCircle(std::vector<Point *> points, const bool &discretization)
{
    std::string name = "l" + std::to_string(lines_.size());
    Line *l = new Circle(name, points, discretization);
    lines_[name] = l;
    return l;
}

Line *Geometry::addSpline(std::vector<Point *> points, const bool &discretization)
{
    std::string name = "l" + std::to_string(lines_.size());
    Line *l = new Spline(name, points, discretization);
    lines_[name] = l;
    return l;
}

LineLoop *Geometry::addLineLoop(std::vector<Line *> lines, const bool &verify)
{
    std::string name = "ll" + std::to_string(lineLoops_.size());
    LineLoop *ll = new LineLoop(name, lines);
    if (verify)
    {
        ll->verification();
    }
    lineLoops_[name] = ll;
    return ll;
}

// Surface *Geometry::addSurface(std::vector<LineLoop *> lineLoop, const int &indexMaterial, const double &thickness)
// {
//     std::string name = "s" + std::to_string(surfaces_.size());
//     Surface *s = new Surface(name, lineLoop, indexMaterial, thickness);
//     surfaces_[name] = s;
//     return s;
// }

Surface *Geometry::addSurface(std::vector<Line *> lines, const int &indexMaterial, const double &thickness)
{
    LineLoop *ll = addLineLoop(lines);
    std::string name = "s" + std::to_string(surfaces_.size());
    Surface *s = new Surface(name, {ll}, indexMaterial, thickness);
    surfaces_[name] = s;
    return s;
}

SurfaceLoop *Geometry::addSurfaceLoop(std::vector<Surface *> surfaces, const bool &verify)
{
    std::string name = "sl" + std::to_string(surfaceLoops_.size());
    SurfaceLoop *sl = new SurfaceLoop(name, surfaces);
    if (verify)
    {
        sl->verification();
    }
    surfaceLoops_[name] = sl;
    return sl;
}

void Geometry::addNeumannCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2, const std::vector<double> &directionX3)
{
    vecDouble force(3);
    force[0] = (directionX1.size() > 0) ? directionX1[0] : 0.0;
    force[1] = (directionX2.size() > 0) ? directionX2[0] : 0.0;
    force[2] = (directionX3.size() > 0) ? directionX3[0] : 0.0;
    neumannConditions_[line->getName()] = force;
}

void Geometry::addNeumannCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2, const std::vector<double> &directionX3)
{
    vecDouble force(3);
    force[0] = (directionX1.size() > 0) ? directionX1[0] : 0.0;
    force[1] = (directionX2.size() > 0) ? directionX2[0] : 0.0;
    force[2] = (directionX3.size() > 0) ? directionX3[0] : 0.0;
    neumannConditions_[point->getName()] = force;
}

void Geometry::addDirichletCondition(Line *line, const std::vector<double> &directionX1, const std::vector<double> &directionX2, const std::vector<double> &directionX3)
{
    vecDouble desloc(3);
    desloc[0] = (directionX1.size() > 0) ? directionX1[0] : 1.0e-240;
    desloc[1] = (directionX2.size() > 0) ? directionX2[0] : 1.0e-240;
    desloc[2] = (directionX3.size() > 0) ? directionX3[0] : 1.0e-240;
    diricheletConditions_[line->getName()] = desloc;
}

void Geometry::addDirichletCondition(Point *point, const std::vector<double> &directionX1, const std::vector<double> &directionX2, const std::vector<double> &directionX3)
{
    vecDouble desloc(3);
    desloc[0] = (directionX1.size() > 0) ? directionX1[0] : 1.0e-240;
    desloc[1] = (directionX2.size() > 0) ? directionX2[0] : 1.0e-240;
    desloc[2] = (directionX3.size() > 0) ? directionX3[0] : 1.0e-240;
    diricheletConditions_[point->getName()] = desloc;
}

std::unordered_map<std::string, vecDouble> Geometry::getDirichletCondition()
{
    return diricheletConditions_;
}

std::unordered_map<std::string, vecDouble> Geometry::getNeumannCondition()
{
    return neumannConditions_;
}

Surface *Geometry::getSurface(const std::string &name)
{
    return surfaces_[name];
}

// void Geometry::addCrackPoint(const std::string &crackName, const vecDouble &coordinates)
// {
//     if (domain_ == "GLOBAL")
//     {
//         int index = points_.size();
//         std::stringstream name;
//         name << "p" << index;
//         Point *p = new Point(name.str(), coordinates, crackes_[crackName]->getLastPoint()->getlcar(), true);
//         points_[name.str()] = p;

//         crackes_[crackName]->addLastPoint(p);
//     }
//     else //LOCAL
//     {
//         std::stringstream name;
//         name << "p" << points_.size();
//         Point *p = new Point(name.str(), coordinates, crackes_[crackName]->getLastPoint()->getlcar(), true);
//         points_[name.str()] = p;

//         crackes_[crackName]->addLastPoint(p);
//     }
// }

// void Geometry::addCrack(std::vector<Point *> points, PlaneSurface *surface, const std::string &openBoundary, const double &jradius, const double &lcarJ)
// {
//     int index = crackes_.size();
//     std::stringstream name;
//     name << "c" << index;

//     Crack *c = new Crack(name.str(), points, surface, openBoundary);
//     crackes_[name.str()] = c;

//     Jradius_ = jradius;
//     lcarJ_ = lcarJ;
// }

// void Geometry::addCrackOnGlobal(std::vector<Point *> points, const std::string &openBoundary, const double &jradius, const double &lcarJ, const double &lengthOffset, const double &lcarOfLocalBoundary)
// {
//     int index = crackes_.size();
//     std::stringstream name;
//     name << "c" << index;

//     PlaneSurface *surface;

//     Crack *c = new Crack(name.str(), points, surface, openBoundary);
//     crackes_[name.str()] = c;

//     Jradius_ = jradius;
//     lcarJ_ = lcarJ;
//     lengthOffset_ = lengthOffset;
//     lcarOfLocalBoundary_ = lcarOfLocalBoundary;
// }

// // std::unordered_map<std::string, Crack *> Geometry::getCrackes()
// // {
// //     return crackes_;
// // }

// double Geometry::getRadiusJintegral()
// {
//     return Jradius_;
// }

void Geometry::transfiniteLine(std::vector<Line *> lines, const int &divisions, const double &progression)
{
    transfiniteLines_.push_back(lines);
    divisionLines_.push_back(divisions);
    progressionLines_.push_back(progression);
}

// // void Geometry::transfiniteSurface(std::vector<PlaneSurface *> planeSurfaces, std::string orientation, std::vector<Point *> points)
// // {
// //     std::stringstream text;
// //     text << "Transfinite Surface {";
// //     for (size_t i = 0; i < planeSurfaces.size(); i++)
// //     {
// //         text << planeSurfaces[i]->getName();
// //         if (i != (planeSurfaces.size() - 1))
// //             text << ", ";
// //     }
// //     text << "} ";
// //     if (points.size() != 0)
// //     {
// //         text << "= {";
// //         for (size_t i = 0; i < points.size(); i++)
// //         {
// //             text << points[i]->getName();
// //             if (i != (points.size() - 1))
// //                 text << ", ";
// //         }
// //         text << "} " << orientation << ";\n//\n";
// //     }
// //     else
// //     {
// //         text << orientation << ";\n//\n";
// //     }
// //     gmshCode_ += text.str();
// // }

std::string Geometry::createGmshCode()
{
    std::string gmshCode;
    for (int in = 0, nPoints = points_.size(); in < nPoints; in++)
    {
        std::string name = "p" + std::to_string(in);
        gmshCode += points_[name]->getGmshCode();
    }

    for (int il = 0, nLines = lines_.size(); il < nLines; il++)
    {
        std::string name = "l" + std::to_string(il);
        gmshCode += lines_[name]->getGmshCode();
    }
    for (int ill = 0, nLineLoops = lineLoops_.size(); ill < nLineLoops; ill++)
    {
        std::string name = "ll" + std::to_string(ill);
        gmshCode += lineLoops_[name]->getGmshCode();
    }
    // for (int is = 0, nSurface = surfaces_.size(); is < nSurface; is++)
    // {
    //     std::string name = "s" + std::to_string(is);
    //     gmshCode += surfaces_[name]->getGmshCode();
    // }
    for (int itl = 0, nTL = transfiniteLines_.size(); itl < nTL; itl++)
    {
        gmshCode += "Transfinite Curve {";
        for (int il = 0, nL = transfiniteLines_[itl].size(); il < nL; il++)
        {
            gmshCode += transfiniteLines_[itl][il]->getName();
            if (nL - il != 1)
            {
                gmshCode += ", ";
            }
        }
        gmshCode += "} = " + std::to_string(divisionLines_[itl] + 1) + " Using Progression " + std::to_string(progressionLines_[itl]) + ";\n//\n";
    }
    return gmshCode;
}

std::unordered_map<std::string, Point *> Geometry::getPoints()
{
    return points_;
}

Point *Geometry::getPoint(const std::string &name)
{
    return points_[name];
}

std::unordered_map<std::string, Surface *> Geometry::getSurfaces()
{
    return surfaces_;
}

int Geometry::getIndexMaterial(const std::string &name)
{
    return surfaces_[name]->getIndexMaterial();
}

std::unordered_map<std::string, std::vector<std::string>> Geometry::getRegions()
{
    std::unordered_map<std::string, std::vector<std::string>> regions;

    for (auto &s : surfaces_)
    {
        std::vector<std::string> linesName;
        for (auto &ll : s.second->getLineLoops())
        {
            for (auto &lines : ll->getLines())
            {
                linesName.push_back(lines->getName());
            }
        }
        regions[s.first] = linesName;
    }
    return regions;
}

bool Geometry::verifyDuplicatedPoint(const std::string &namePoint)
{
    bool duplicated = false;
    vecDouble coordPoint(3), coordAux(3);
    points_[namePoint]->getCoordinates(coordPoint);
    for (auto &point : points_)
    {
        if (point.first != namePoint)
        {
            point.second->getCoordinates(coordAux);
            if (fabs(coordPoint[0] - coordAux[0]) <= 1.0e-12 and fabs(coordPoint[1] - coordAux[1]) <= 1.0e-12)
                duplicated = true;
        }
    }
    return duplicated;
}

// // int Geometry::getNumberOfPoints()
// // {
// //     return points_.size();
// // }

// // int Geometry::getNumberOfCrackes()
// // {
// //     return crackes_.size();
// // }

// void Geometry::createGeometryFromCrack()
// {
//     cout << "IMPLEMENTAR FUNÇÃO createGeometryFromCrack()\n";
// }

// // void Geometry::createGeometryFromCrack()
// // {
// //     if (remeshNumber_ == 0) //Primeira malha criada
// //     {
// //         for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes_.begin(); crack != crackes_.end(); crack++)
// //         {
// //             std::vector<Point *> crackPoints = crack->second->getPoints();
// //             std::string open = crack->second->getOpenBoundary();

// //             if (crackPoints.size() > 2)
// //             {
// //                 std::cout << "NÃO FOI IMPLEMENTADO FISSURAS INICIAIS COM MAIS DE DOIS PONTOS." << std::endl;
// //             }
// //             double initialAngle = crack->second->getLastAngle();
// //             const double pi = 3.14159265359;

// //             bounded_vector<double, 2> coord1, coord2, coordInitial, coordLast;

// //             coordInitial = crackPoints[0]->getCoordinates();
// //             coordLast = crackPoints[1]->getCoordinates();

// //             if (open == "first")
// //             {
// //                 coord1 = coordInitial;
// //                 coord2 = coordInitial;

// //                 coord1(0) += cos(initialAngle - 0.5 * pi) * lengthOffset_;
// //                 coord1(1) += sin(initialAngle - 0.5 * pi) * lengthOffset_;

// //                 coord2(0) += cos(initialAngle + 0.5 * pi) * lengthOffset_;
// //                 coord2(1) += sin(initialAngle + 0.5 * pi) * lengthOffset_;
// //             }
// //             else
// //             {
// //                 coord1 = coordInitial;
// //                 coord2 = coordInitial;

// //                 double length = lengthOffset_ / sin(0.25 * pi);

// //                 coord1(0) += cos(initialAngle - 0.75 * pi) * length;
// //                 coord1(1) += sin(initialAngle - 0.75 * pi) * length;

// //                 coord2(0) += cos(initialAngle + 0.75 * pi) * length;
// //                 coord2(1) += sin(initialAngle + 0.75 * pi) * length;
// //             }

// //             Point *paux1_1 = addPoint({coord1(0), coord1(1)}, lcarOfLocalBoundary_, false);
// //             Point *paux1_2 = addPoint({coord2(0), coord2(1)}, lcarOfLocalBoundary_, false);

// //             crack->second->addLastPointsOfLocalGeometry(paux1_1, paux1_2);

// //             if (open == "second")
// //             {
// //                 coord1 = coordLast;
// //                 coord2 = coordLast;

// //                 coord1(0) += cos(initialAngle - 0.5 * pi) * lengthOffset_;
// //                 coord1(1) += sin(initialAngle - 0.5 * pi) * lengthOffset_;

// //                 coord2(0) += cos(initialAngle + 0.5 * pi) * lengthOffset_;
// //                 coord2(1) += sin(initialAngle + 0.5 * pi) * lengthOffset_;
// //             }
// //             else
// //             {
// //                 coord1 = coordLast;
// //                 coord2 = coordLast;

// //                 double length = lengthOffset_ / sin(0.25 * pi);

// //                 coord1(0) += cos(initialAngle - 0.25 * pi) * length;
// //                 coord1(1) += sin(initialAngle - 0.25 * pi) * length;

// //                 coord2(0) += cos(initialAngle + 0.25 * pi) * length;
// //                 coord2(1) += sin(initialAngle + 0.25 * pi) * length;
// //             }
// //             Point *paux1 = addPoint({coord1(0), coord1(1)}, lcarOfLocalBoundary_, false);
// //             Point *paux2 = addPoint({coord2(0), coord2(1)}, lcarOfLocalBoundary_, false);

// //             crack->second->addLastPointsOfLocalGeometry(paux1, paux2);

// //             std::pair<std::vector<Point *>, std::vector<Point *>> boundaryPoints = crack->second->getPointsOfLocalGeometry();

// //             std::vector<Point *> conec;

// //             int aux = boundaryPoints.first.size();
// //             if (open == "first")
// //             {
// //                 bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[aux - 1]->getCoordinates() + boundaryPoints.second[aux - 1]->getCoordinates());
// //                 Point *paux = addPoint({auxCoord(0), auxCoord(1)}, lcarOfLocalBoundary_, false);
// //                 crack->second->setAuxPoint(paux);
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.first[i]);
// //                 }
// //                 conec.push_back(paux);
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                 }
// //             }
// //             else if (open == "second")
// //             {
// //                 bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[0]->getCoordinates() + boundaryPoints.second[0]->getCoordinates());
// //                 Point *paux = addPoint({auxCoord(0), auxCoord(1)}, lcarOfLocalBoundary_, false);
// //                 crack->second->setAuxPoint(paux);

// //                 conec.push_back(crackPoints[1]);
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                 }
// //                 conec.push_back(paux);
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.first[i]);
// //                 }
// //                 conec.push_back(crackPoints[1]);
// //             }
// //             else
// //             {
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.first[i]);
// //                 }
// //                 for (int i = 0; i < aux; i++)
// //                 {
// //                     conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                 }
// //                 conec.push_back(boundaryPoints.first[0]);
// //             }

// //             Line *l0 = addBSpline(conec, true);
// //             LineLoop *ll;
// //             if (open == "first")
// //             {
// //                 Line *l1 = addLine({boundaryPoints.second[0], crackPoints[0]});
// //                 Line *l2 = addLine({crackPoints[0], boundaryPoints.first[0]});
// //                 ll = addLineLoop({l0, l1, l2});
// //             }
// //             else if (open == "second")
// //             {
// //                 Line *l1 = addLine({boundaryPoints.first[aux - 1], crackPoints[aux - 1]});
// //                 Line *l2 = addLine({crackPoints[aux - 1], boundaryPoints.second[aux - 1]});
// //                 ll = addLineLoop({l0, l1, l2});
// //             }
// //             else
// //             {
// //                 ll = addLineLoop({l0});
// //             }

// //             PlaneSurface *s = addPlaneSurface({ll}, indexMaterial_, thickness_);
// //             crack->second->setSurface(s);
// //         }
// //         remeshNumber_ += 1;
// //     }
// //     else //AQUI PRECISAMOS ATUALIZAR A POSIÇÃO DO ÚLTIMO BOUNDARYPOINTS E ADICIONAR MAIS UM
// //     {
// //         for (std::unordered_map<std::string, Crack *>::const_iterator crack = crackes_.begin(); crack != crackes_.end(); crack++)
// //         {
// //             std::vector<Point *> crackPoints = crack->second->getPoints();
// //             std::string open = crack->second->getOpenBoundary();

// //             std::pair<std::vector<Point *>, std::vector<Point *>> boundaryPoints = crack->second->getPointsOfLocalGeometry();

// //             if (crackPoints.size() > boundaryPoints.first.size() and crackPoints.size() > boundaryPoints.second.size()) //FOI ADICIONADO OUTRO PONTO NA FISSURA, LOGO PRECISAMOS RECONSTRUIR O CONTORNO
// //             {
// //                 const double pi = 3.14159265359;
// //                 bounded_vector<double, 2> coord0, coord1, coord2;
// //                 int aux = crackPoints.size();
// //                 coord0 = crackPoints[aux - 3]->getCoordinates();
// //                 coord1 = crackPoints[aux - 2]->getCoordinates();
// //                 coord2 = crackPoints[aux - 1]->getCoordinates();

// //                 double teta1, teta0;

// //                 teta0 = atan2(coord1(1) - coord0(1), coord1(0) - coord0(0));
// //                 teta1 = atan2(coord2(1) - coord1(1), coord2(0) - coord1(0));

// //                 bounded_matrix<double, 2, 2> mat;
// //                 mat(0, 0) = cos(teta0);
// //                 mat(1, 0) = sin(teta0);
// //                 mat(0, 1) = -cos(teta1);
// //                 mat(1, 1) = -sin(teta1);

// //                 //MODIFICANDO ÚLTIMO DO CONTORNO 1
// //                 bounded_vector<double, 2> boundCoord0, auxCoord1;
// //                 boundCoord0 = boundaryPoints.first[aux - 3]->getCoordinates();
// //                 auxCoord1(0) = coord1(0) + cos(teta1 - 0.5 * pi) * lengthOffset_;
// //                 auxCoord1(1) = coord1(1) + sin(teta1 - 0.5 * pi) * lengthOffset_;

// //                 bounded_vector<double, 2> auxDelta = auxCoord1 - boundCoord0;
// //                 bounded_vector<double, 2> values = prod(inverseMatrix(mat), auxDelta);

// //                 if (values(0) < 0.0)
// //                 {
// //                     values(0) = 0.0;
// //                 }

// //                 bounded_vector<double, 2> newCoord;
// //                 newCoord(0) = boundCoord0(0) + cos(teta0) * values(0);
// //                 newCoord(1) = boundCoord0(1) + sin(teta0) * values(0);
// //                 boundaryPoints.first[aux - 2]->setCoordinates(newCoord);

// //                 //MODIFICANDO ÚLTIMO DO CONTORNO 2
// //                 boundCoord0 = boundaryPoints.second[aux - 3]->getCoordinates();
// //                 auxCoord1(0) = coord1(0) + cos(teta1 + 0.5 * pi) * lengthOffset_;
// //                 auxCoord1(1) = coord1(1) + sin(teta1 + 0.5 * pi) * lengthOffset_;

// //                 auxDelta = auxCoord1 - boundCoord0;
// //                 values = prod(inverseMatrix(mat), auxDelta);

// //                 if (values(0) < 0.0)
// //                 {
// //                     values(0) = 0.0;
// //                 }

// //                 newCoord(0) = boundCoord0(0) + cos(teta0) * values(0);
// //                 newCoord(1) = boundCoord0(1) + sin(teta0) * values(0);
// //                 boundaryPoints.second[aux - 2]->setCoordinates(newCoord);

// //                 //CRIANDO DOIS ÚLTIMOS PONTOS DO CONTORNO
// //                 bounded_vector<double, 2> coord2_1, coord2_2;
// //                 coord2_1 = coord2;
// //                 coord2_2 = coord2;

// //                 double length = lengthOffset_ / sin(0.25 * pi);

// //                 coord2_1(0) += cos(teta1 - 0.25 * pi) * length;
// //                 coord2_1(1) += sin(teta1 - 0.25 * pi) * length;

// //                 coord2_2(0) += cos(teta1 + 0.25 * pi) * length;
// //                 coord2_2(1) += sin(teta1 + 0.25 * pi) * length;

// //                 Point *paux1 = addPoint({coord2_1(0), coord2_1(1)}, lcarOfLocalBoundary_, false);
// //                 Point *paux2 = addPoint({coord2_2(0), coord2_2(1)}, lcarOfLocalBoundary_, false);
// //                 crack->second->addLastPointsOfLocalGeometry(paux1, paux2);

// //                 //EDITANDO BSPLINE DO CONTORNO LOCAL

// //                 std::vector<Point *> conec;

// //                 boundaryPoints = crack->second->getPointsOfLocalGeometry();

// //                 if (open == "first") //IMPLEMENTADO SOMENTE PARA ESSE TIPO DE PROPAGAÇÃO!
// //                 {
// //                     bounded_vector<double, 2> auxCoord = 0.5 * (boundaryPoints.first[aux - 1]->getCoordinates() + boundaryPoints.second[aux - 1]->getCoordinates());
// //                     crack->second->getAuxPoint()->setCoordinates(auxCoord);
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.first[i]);
// //                     }
// //                     conec.push_back(crack->second->getAuxPoint());
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                     }
// //                 }
// //                 else if (open == "second")
// //                 {
// //                     conec.push_back(crackPoints[1]);
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                     }
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.first[i]);
// //                     }
// //                     conec.push_back(crackPoints[1]);
// //                 }
// //                 else
// //                 {
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.first[i]);
// //                     }
// //                     for (int i = 0; i < aux; i++)
// //                     {
// //                         conec.push_back(boundaryPoints.second[aux - 1 - i]);
// //                     }
// //                     conec.push_back(boundaryPoints.first[0]);
// //                 }

// //                 crack->second->getPlaneSurface()->getLineLoop(0)->getLine(0)->setPoints(conec);
// //             }
// //             else
// //             {
// //                 break;
// //             }
// //         }
// //         remeshNumber_ += 1;
// //     }
// // }