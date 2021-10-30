draw_cyl_with_helix(3, 10, 45, 0, 0.1)
%calc_ota(3, 10, 1, 45, 0.1)

function draw_cyl_with_helix(r_cyl,l_cyl,wind_angle, r_axis, fiber_width)
alpha = deg2rad(wind_angle); % wind angle, radians
t_size = 1000;

%phi0 = 0; % starting phi
%zeta = 0.5 * h / r;
%theta0 = atan(1 / zeta); % starting theta
%theta_span = theta0:0.01:(pi / 2); % theta span vectors

t = linspace(0,l_cyl,t_size);
r = zeros(size(t));
theta = zeros(size(t));
x = zeros(size(t));

for i = 1:t_size
    r(i) = r_cyl;
    theta(i) = (t(i) * pi) / (2 * r_cyl * cot(alpha));
    x(i) = t(i);
end

r2 = zeros(2); %stores the rs needed to plot OTA line
theta2 = zeros(2); %stores thetas needed for OTA line
x2 = zeros(2); %stores xs for OTA

%first point for OTA is last point from helix curve
r2(1) = r_cyl;
theta2(1) = theta(t_size);
x2(1) = x(t_size);

%second point for OTA is the first point but the theta moves forward one
%OTA
r2(2) = r_cyl;
theta2(2) = theta2(1) + calc_ota(r_cyl, l_cyl, wind_angle, r_axis, fiber_width);
x2(2) = x2(1);


y = r .* cos(theta);
z = r .* sin(theta);

y2 = r2 .* cos(theta2);
z2 = r2 .* sin(theta2);

r3_axis = zeros(32);
theta_axis = linspace(0, 2 * pi, 32);
x_axis = zeros(32);

r3_axis = r3_axis + r_axis;
x_axis = x_axis + l_cyl;
y_axis = r3_axis .* cos(theta_axis);
z_axis = r3_axis .* sin(theta_axis);

plot3(x,y,z);
axis equal
hold on
plot3(x2, y2, z2);
plot3(x_axis, y_axis, z_axis);
[y_cyl,z_cyl,x_cyl] = cylinder(r_cyl);
x_cyl = x_cyl * l_cyl;
mesh(x_cyl,y_cyl,z_cyl);
end

function ota = calc_ota(r_cyl, l_cyl, wind_angle, r_axis, fiber_width);
%turn angle into radians
wind_angle = deg2rad(wind_angle);
%find the apparent width of the fiber at the corner of the cylinder-
%accounts for the wind angle
width_relative = sec(wind_angle) * fiber_width;
%calculate the number of fibers (F) that can fit around the circumference
%of the cylinder
num_fibers = 2 * pi * r_cyl / width_relative;
%make this an integer (ceil rather than floor bc overlap better than gap)
num_fibers = ceil(num_fibers)
%guess an initial OTA using just the maximum possible OTA (known from
%axis/cyl radius)
ota = 2 * acos(r_axis / r_cyl)
%calculate spiral angle (taken from helix parametric equation)
spiral_angle = l_cyl * pi / (2 * r_cyl * cot(wind_angle))
%calculate fiber steps (k) undertaken every 2OTA + 2spiral cycle.
%total angle spun every cycle
total_angle = mod(2 * (spiral_angle + ota), 2 * pi);
%K = 2(OTA + Spiral) * F/2pi
fiber_step = total_angle * num_fibers / (2 * pi);
%turn this into an integer
fiber_step = floor(fiber_step);
%find step size for which relative prime rules apply
fiber_step = find_fiber_step(fiber_step, num_fibers);
%solve for ota with final calculated step size
ota = (fiber_step * pi / num_fibers) - spiral_angle;
return
end

function fiber_step = find_fiber_step(fiber_step, num_fibers)
    GCD = gcd(fiber_step, num_fibers);
    if GCD == 0
        return
    end
    
    fiber_step_inc = fiber_step;
    fiber_step_dec = fiber_step;
    while true
        fiber_step_inc = fiber_step_inc + 1;
        GCD_inc = gcd(fiber_step_inc, num_fibers);
        if GCD_inc == 1
            fiber_step = fiber_step_inc;
            return
        end
        
        fiber_step_dec = fiber_step_dec - 1;
        GCD_dec = gcd(fiber_step_dec, num_fibers);
        if GCD_dec == 1
            fiber_step = fiber_step_dec;
            return
        end
    end
end