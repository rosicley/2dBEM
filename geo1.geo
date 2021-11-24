p0 = newp; Point(p0) = {0.000000, 0.000000, 0.0, 1.000000}; Physical Point('p0') = {p0};
//
p1 = newp; Point(p1) = {10.000000, 0.000000, 0.0, 1.000000}; Physical Point('p1') = {p1};
//
p2 = newp; Point(p2) = {10.000000, 0.000000, 0.0, 1.000000}; Physical Point('p2') = {p2};
//
p3 = newp; Point(p3) = {10.000000, 1.000000, 0.0, 1.000000}; Physical Point('p3') = {p3};
//
p4 = newp; Point(p4) = {10.000000, 1.000000, 0.0, 1.000000}; Physical Point('p4') = {p4};
//
p5 = newp; Point(p5) = {0.000000, 1.000000, 0.0, 1.000000}; Physical Point('p5') = {p5};
//
p6 = newp; Point(p6) = {0.000000, 1.000000, 0.0, 1.000000}; Physical Point('p6') = {p6};
//
p7 = newp; Point(p7) = {0.000000, 0.000000, 0.0, 1.000000}; Physical Point('p7') = {p7};
//
l0 = newl; Line(l0) = {p0, p1}; Physical Line('l0') = {l0};
//
l1 = newl; Line(l1) = {p2, p3}; Physical Line('l1') = {l1};
//
l2 = newl; Line(l2) = {p4, p5}; Physical Line('l2') = {l2};
//
l3 = newl; Line(l3) = {p6, p7}; Physical Line('l3') = {l3};
//
Transfinite Curve {l0, l2} = 51 Using Progression 1;
//
Transfinite Curve {l1, l3} = 6 Using Progression 1;
//
Mesh.ElementOrder = 3;
//
Mesh.Algorithm = 2;
//
Mesh 1;
//
Mesh.MshFileVersion = 2.2;
Save "geo1.msh";
