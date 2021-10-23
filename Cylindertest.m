%draw a cylinder
r = 1;
h = 10;

cyl_maker(r,h)

function cyl_maker(r,h)
[X,Y,Z] = cylinder(r);
Z = Z * h;
surf(X,Y,Z);
axis equal
end