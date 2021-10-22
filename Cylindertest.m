%draw a cylinder
r = 1;
[X,Y,Z] = cylinder(r);
h = 10;
Z = Z * h;
surf(X,Y,Z);
axis equal