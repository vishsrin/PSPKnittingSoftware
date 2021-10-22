r_cyl = 0.5; % radius of cylinder, in
h = 10; % height of cylinder, in
alpha = deg2rad(20); % wind angle, radians
t_size = 1000;

%phi0 = 0; % starting phi
%zeta = 0.5 * h / r;
%theta0 = atan(1 / zeta); % starting theta
%theta_span = theta0:0.01:(pi / 2); % theta span vectors

t = linspace(0,h,t_size);
r = zeros(size(t));
theta = zeros(size(t));
z = zeros(size(t));

for i = 1:length(t)
    r(i) = r_cyl;
    theta(i) = (t(i) * pi) / (2 * r_cyl * sin(alpha));
    z(i) = t(i);
end

x = r .* cos(theta);
y = r .* sin(theta);

plot3(x,y,z);
axis equal
hold on
[x_cyl,y_cyl,z_cyl] = cylinder(r_cyl);
z_cyl = z_cyl * h;
mesh(x_cyl,y_cyl,z_cyl);