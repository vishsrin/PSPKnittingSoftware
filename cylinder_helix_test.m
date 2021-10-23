draw_cyl_with_helix(3,10,45)

function draw_cyl_with_helix(r_cyl,l,alpha)
alpha = deg2rad(alpha); % wind angle, radians
t_size = 1000;

%phi0 = 0; % starting phi
%zeta = 0.5 * h / r;
%theta0 = atan(1 / zeta); % starting theta
%theta_span = theta0:0.01:(pi / 2); % theta span vectors

t = linspace(0,l,t_size);
r = zeros(size(t));
theta = zeros(size(t));
z = zeros(size(t));

for i = 1:length(t)
    r(i) = r_cyl;
    theta(i) = (t(i) * pi) / (2 * r_cyl * sin(alpha));
    x(i) = t(i);
end

y = r .* cos(theta);
z = r .* sin(theta);

plot3(x,y,z);
axis equal
hold on
[y_cyl,z_cyl,x_cyl] = cylinder(r_cyl);
x_cyl = x_cyl * l;
mesh(x_cyl,y_cyl,z_cyl);
end