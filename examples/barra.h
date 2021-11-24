Geometry *geo = new Geometry("geo1"); //LOCAL OR GLOBAL

double l = 10.0, h = 1.0;
Point *p0 = geo->addPoint({0.0, 0.0});
Point *p1 = geo->addPoint({l, 0.0});
Point *p1c = geo->addPoint({l, 0.0});
Point *p2 = geo->addPoint({l, h});
Point *p2c = geo->addPoint({l, h});
Point *p3 = geo->addPoint({0.0, h});
Point *p3c = geo->addPoint({0.0, h});
Point *p0c = geo->addPoint({0.0, 0.0});
Line *l0 = geo->addLine({p0, p1});
Line *l1 = geo->addLine({p1c, p2});
Line *l2 = geo->addLine({p2c, p3});
Line *l3 = geo->addLine({p3c, p0c});

geo->transfiniteLine({l0, l2}, 50);
geo->transfiniteLine({l1, l3}, 5);

geo->addNeumannCondition(l1, {}, {-100.00});
geo->addDirichletCondition(l3, {0.0}, {0.0});
// geo->addDirichletCondition(l0, {}, {0.0});
// geo->addDirichletCondition(l2, {}, {0.0});
// geo->addDirichletCondition(l1, {},{0.0});

Problem *problem = new Problem(40, 0.25, 0.0);

problem->addMaterial(1000.0e03, 0.35);

// problem->addInternalPoints({{2.0, 3.0}, {3.0, 3.0}, {4.0, 3.0}, {1.0, 3.0}, {5.0, 3.0}});

problem->generateMesh(geo, 3, "AUTO", false); // ordem <= 10

problem->solveElasticityProblem("EPD");

problem->teste();