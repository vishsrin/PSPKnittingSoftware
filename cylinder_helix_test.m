draw_wind_pass(3, 20, 45, 1, 01)

function theta_step = draw_wind_pass(r_cyl, l_cyl, wind_angle, r_axis, fiber_width)
wind_angle = deg2rad(wind_angle); % wind angle, radians
t_size = 1000;

%phi0 = 0; % starting phi
%zeta = 0.5 * h / r;
%theta0 = atan(1 / zeta); % starting theta
%theta_span = theta0:0.01:(pi / 2); % theta span vectors

t = linspace(0,l_cyl,t_size); 
r = zeros(size(t));
theta_up = zeros(size(t));
theta_down = zeros(size(t));
x = zeros(size(t));
r_cyl;
l_cyl;
wind_angle;
r_axis;
fiber_width;
wind_angle = calc_wind_angle(r_cyl, l_cyl, wind_angle, fiber_width);
ota = pi - (2 * wind_angle);


theta_step = 0;
num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width)

for p = 1:num_fibers
    p;
    theta_step;
    for i = 1:t_size % drawing #fibers times
        r(i) = r_cyl;
        theta_up(i) = (t(i) * pi) / (2 * r_cyl * cot(wind_angle)) + theta_step; %theta_up values
        x(i) = t(i);
    end
    
    r_coord_ota = zeros(2); %stores the rs needed to plot OTA line
    theta_coord_ota = zeros(2); %stores thetas needed for OTA line
    x_coord_ota = zeros(2); %stores xs for OTA
    
    %first point for OTA is last point from helix curve
    r_coord_ota(1) = r_cyl;
    theta_coord_ota(1) = theta_up(t_size);
    x_coord_ota(1) = x(t_size);
    
    %second point for OTA is the first point but the theta moves forward one
    %OTA
    r_coord_ota(2) = r_cyl;
    theta_coord_ota(2) = theta_up(t_size) + ota;
    x_coord_ota(2) = x_coord_ota(1);
    
    % since theta_down pass moves the other way, array must be indexed
    % backwards
    for i = 1:t_size
        n = t_size - i + 1;
        theta_down(i) = (t(n) * pi) / (2 * r_cyl * cot(wind_angle)) + theta_coord_ota(2);
    end
    
    y_coord_up = r .* cos(theta_up);
    z_coord_up = r .* sin(theta_up);
    
    y_coord_down = r .* cos(theta_down);
    z_coord_down = r .* sin(theta_down);
    
    y_coord_ota = r_coord_ota .* cos(theta_coord_ota);
    z_coord_ota = r_coord_ota .* sin(theta_coord_ota);
    
    plot3(x,y_coord_up,z_coord_up); % cylinder up path
    axis equal
    hold on
    plot3(x, y_coord_down, z_coord_down); % cylinder down path
    plot3(x_coord_ota, y_coord_ota, z_coord_ota); % ota up path
    
    % re-assigning OTA values for down-OTA
    r_coord_ota(1) = r_cyl;
    theta_coord_ota(1) = theta_down(1);
    x_coord_ota(1) = x(1);
    
    r_coord_ota(2) = r_cyl;
    theta_coord_ota(2) = theta_down(1) + ota;
    x_coord_ota(2) = x_coord_ota(1);
    
    y_coord_ota = r_coord_ota .* cos(theta_coord_ota);
    z_coord_ota = r_coord_ota .* sin(theta_coord_ota);
    
    plot3(x_coord_ota, y_coord_ota, z_coord_ota); % OTA down path
    theta_step = theta_down(1) + ota;
end

[y_cyl,z_cyl,x_cyl] = cylinder(r_cyl); % creating cylinder object
x_cyl = x_cyl * l_cyl; % scaling
mesh(x_cyl,y_cyl,z_cyl);

return
end

function num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width)
width_relative = sec(wind_angle) * fiber_width;
num_fibers = 2 * pi * r_cyl / width_relative;
num_fibers = ceil(num_fibers);
end

function new_theta = calc_wind_angle(r_cyl, l_cyl, wind_angle, fiber_width)
theta_inc = 0.01;
theta_plus = wind_angle;
theta_minus = theta_plus;

num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width);
fiber_step = calc_fiber_step(r_cyl, l_cyl, wind_angle, fiber_width);

while (calc_step_prime(fiber_step, num_fibers))
    num_fibers = calc_num_fibers(r_cyl, theta_plus, fiber_width);
    fiber_step = calc_fiber_step(r_cyl, l_cyl, theta_plus, fiber_width);
    theta_plus = theta_plus + theta_inc;
end

num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width);
fiber_step = calc_fiber_step(r_cyl, l_cyl, wind_angle, fiber_width);

while (calc_step_prime(fiber_step, num_fibers))
    num_fibers = calc_num_fibers(r_cyl, theta_minus, fiber_width);
    fiber_step = calc_fiber_step(r_cyl, l_cyl, theta_minus, fiber_width);
    theta_minus = theta_minus - theta_inc;
end

if (theta_plus - wind_angle) > (wind_angle - theta_minus)
    new_theta = theta_minus;
else
    new_theta = theta_plus;
end

return
end

function fiber_step = calc_fiber_step(r_cyl, l_cyl, theta, fiber_width)

num_fibers = calc_num_fibers(r_cyl, theta, fiber_width);
%guess an initial OTA using just the maximum possible OTA (known from
%axis/cyl radius)
ota = pi - (2 * theta);
%calculate spiral angle (taken from helix parametric equation)
spiral_angle = l_cyl * pi / (2 * r_cyl * cot(theta));
%calculate fiber steps (k) undertaken every 2OTA + 2spiral cycle.
%total angle spun every cycle
total_angle = mod(2 * (spiral_angle + ota), 2 * pi);
%K = 2(OTA + Spiral) * F/2pi
fiber_step = total_angle * num_fibers / (2 * pi);
%turn this into an integer
fiber_step = floor(fiber_step);

return
end

function is_prime = calc_step_prime(fiber_step, num_fibers)
    GCD = gcd(fiber_step, num_fibers);
    if GCD == 0
        return
    end
    
    if GCD == 1
        is_prime = 1;
    else
        is_prime = 0;
    end
return
end