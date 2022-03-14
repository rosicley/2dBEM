#pragma once
#include "../src/Problem.h"

Geometry *geo = new Geometry("geo1"); //LOCAL OR GLOBAL
double l = 1.0, h = 1.0;

///CRIANDO SUB-REGIÃO 1
///PONTOS
Point *p0 = geo->addPoint({0.0, 0.0});
Point *p1 = geo->addPoint({l * 2.0, 0.0});
Point *p2 = geo->addPoint({l * 2.0, 0.0});
Point *p3 = geo->addPoint({l * 2.0, h});
Point *p4 = geo->addPoint({l * 2.0, h});
Point *p5 = geo->addPoint({0.0, h});
Point *p6 = geo->addPoint({0.0, h});
Point *p7 = geo->addPoint({0.0, 0.0});
///LINHAS
Line *l0 = geo->addLine({p0, p1});
Line *l1 = geo->addLine({p2, p3});
Line *l2 = geo->addLine({p4, p5});
Line *l3 = geo->addLine({p6, p7});
///SUPERFÍCIE
Surface *s0 = geo->addSurface({l0, l1, l2, l3}, 1);


///CRIANDO SUB-REGIÃO 2
///PONTOS
Point *p8 = geo->addPoint({l * 2.0, 0.0});
Point *p9 = geo->addPoint({l * 4.0, 0.0});
Point *p10 = geo->addPoint({l * 4.0, 0.0});
Point *p11 = geo->addPoint({l * 4.0, h});
Point *p12 = geo->addPoint({l * 4.0, h});
Point *p13 = geo->addPoint({l * 2.0, h});
Point *p14 = geo->addPoint({l * 2.0, h});
Point *p15 = geo->addPoint({l * 2.0, 0.0});
///LINHAS
Line *l4 = geo->addLine({p8, p9});
Line *l5 = geo->addLine({p10, p11});
Line *l6 = geo->addLine({p12, p13});
Line *l7 = geo->addLine({p14, p15});
///SUPERFÍCIE
Surface *s1 = geo->addSurface({l4, l5, l6, l7}, 0);


///ESCOLHENDO A QUANTIDADE DE LINHAS QUE O GMSH DIVIDIRÁ CADA LINHA
geo->transfiniteLine({l0, l2, l4, l6}, 2);
geo->transfiniteLine({l1, l3, l5, l7}, 1);


///APLICANDO CONDIÇÕES DE CONTORNO
geo->addDirichletCondition(l3, {0.0}, {0.0});
// geo->addNeummanCondition(l5, {1.0}, {});


///CRIANDO O OBJETO DO PROBLEMA
Problem *problem = new Problem(40, 0.25, 0.0);


///ADICIONANDO MATERIAIS
problem->addMaterial(2.0, 0.0);
problem->addMaterial(5.0, 0.0);


///GERANDO A MALHA GEOMÉTRICA E DE COLOCAÇÃO
problem->generateMesh(geo, 3);


///ADICIONANDO PONTOS INTERNOS
problem->addInternalPoints(s0, {{0.25, 0.5}, {0.5, 0.5}, {0.75, 0.5}, {1.0, 0.5}, {1.25, 0.5}, {1.5, 0.5}, {1.75, 0.5}});
problem->addInternalPoints(s1, {{2.25, 0.5}, {2.5, 0.5}, {2.75, 0.5}, {3.0, 0.5}, {3.25, 0.5}, {3.5, 0.5}, {3.75, 0.5}});


///ADICIONANDO FORÇA DE DOMÍNIO À SUPERFÍCIE
problem->addBodyForceToSurface(s1, {1.0, 0.0});


///ACOPLANDO AS LINHAS
problem->coupleLines({{l1, l7}});


///SOLUCIONANDO O PROBLEMA
problem->solveElasticityProblem("EPD");