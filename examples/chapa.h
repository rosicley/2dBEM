Geometry *geo = new Geometry("geo1"); //LOCAL OR GLOBAL

double l = 1.0, h = 1.0;
Point *p0 = geo->addPoint({0.0, 0.0});
Point *p1 = geo->addPoint({l, 0.0});
Point *p2 = geo->addPoint({l, 0.0});
Point *p3 = geo->addPoint({l, h});
Point *p4 = geo->addPoint({l, h});
Point *p5 = geo->addPoint({0.0, h});
Point *p6 = geo->addPoint({0.0, h});
Point *p7 = geo->addPoint({0.0, 0.0});
Line *l0 = geo->addLine({p0, p1});
Line *l1 = geo->addLine({p2, p3});
Line *l2 = geo->addLine({p4, p5});
Line *l3 = geo->addLine({p6, p7});
Surface *s0 = geo->addSurface({l0, l1, l2, l3});

Point *p8 = geo->addPoint({l, 0.0});
Point *p9 = geo->addPoint({l * 2.0, 0.0});
Point *p10 = geo->addPoint({l * 2.0, 0.0});
Point *p11 = geo->addPoint({l * 2.0, h});
Point *p12 = geo->addPoint({l * 2.0, h});
Point *p13 = geo->addPoint({l, h});
Point *p14 = geo->addPoint({l, h});
Point *p15 = geo->addPoint({l, 0.0});
Line *l4 = geo->addLine({p8, p9});
Line *l5 = geo->addLine({p10, p11});
Line *l6 = geo->addLine({p12, p13});
Line *l7 = geo->addLine({p14, p15});
Surface *s1 = geo->addSurface({l4, l5, l6, l7});

geo->transfiniteLine({l0, l1, l3, l2}, 5);
geo->transfiniteLine({l4, l5, l6, l7}, 5);

geo->addNeumannCondition(l2, {}, {1.0});
geo->addNeumannCondition(l6, {}, {1.0});
geo->addDirichletCondition(l3, {0.0}, {});
geo->addDirichletCondition(l0, {}, {0.0});
geo->addDirichletCondition(l4, {}, {0.0});

// geo->addNeumannCondition(l5, {}, {0.010});
// geo->addDirichletCondition(l3, {0.0}, {0.0});

Problem *problem = new Problem(40, 0.25, 0.0);

problem->addMaterial(1.0, 0.0);

// problem->addInternalPoints({{2.0, 3.0}, {3.0, 3.0}, {4.0, 3.0}, {1.0, 3.0}, {5.0, 3.0}});

problem->generateMesh(geo, 3, "AUTO", false, true); // ordem <= 10

problem->addInternalPoints(s0, {{0.25, 0.5}, {0.5, 0.5}, {0.75, 0.5}});
problem->addInternalPoints(s1, {{1.5, 0.25}, {1.5, 0.75}});

// problem->addBodyForceToSurface(s0, {1.0, 0.0});
// problem->addBodyForceToSurface(s1, {1.0, 0.0});

problem->coupleLines({{l1, l7}});

problem->solveElasticityProblem("EPD");