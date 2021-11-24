Geometry *geo = new Geometry("geo1"); //LOCAL OR GLOBAL

double l = 6.0;
Point *p0 = geo->addPoint({0.0, 0.0});
Point *p1 = geo->addPoint({l, 0.0});
Point *p1c = geo->addPoint({l, 0.0});
Point *p2 = geo->addPoint({l, l});
Point *p2c = geo->addPoint({l, l});
Point *p3 = geo->addPoint({0.0, l});
Point *p3c = geo->addPoint({0.0, l});
Point *p0c = geo->addPoint({0.0, 0.0});

Line *l0 = geo->addLine({p0, p1});
Line *l1 = geo->addLine({p1c, p2});
Line *l2 = geo->addLine({p2c, p3});
Line *l3 = geo->addLine({p3c, p0c});

geo->transfiniteLine({l0, l1, l2, l3}, 5);

geo->addNeumannCondition(l0, {0.0});
geo->addNeumannCondition(l1, {-50.0});
geo->addNeumannCondition(l2, {0.0});

geo->addDirichletCondition(l3, {300.0});

Problem *problem = new Problem(12, 0.25, 1.0);

// problem->addMaterial(100.0);

problem->addInternalPoints({{2.0, 3.0}, {3.0, 3.0}, {4.0, 3.0}, {1.0, 3.0}, {5.0, 3.0}});

problem->generateMesh(geo, 5, "AUTO", true); // ordem <= 10

problem->solvePotentialProblem();

problem->teste();