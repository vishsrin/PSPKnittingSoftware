draw_full_path(3, 20, 33, 1);

function draw_full_path(r_cyl, l_cyl, wind_angle, fiber_width)
wind_angle = deg2rad(wind_angle);

[y_cyl,z_cyl,x_cyl] = cylinder(r_cyl); % creating cylinder object
x_cyl = x_cyl * l_cyl; % scaling
mesh(x_cyl,y_cyl,z_cyl);
axis equal
hold on


num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width);

path_coordinates = draw_one_path(r_cyl, l_cyl, wind_angle, fiber_width, 0, num_fibers);

step_angle = path_coordinates(2, length(path_coordinates(2, :)));

for i = 1 : 1
path_coordinates = horzcat(path_coordinates, draw_one_path(r_cyl, l_cyl, wind_angle, fiber_width, i * step_angle, num_fibers));

graph_path(path_coordinates);

end
end

function path_coordinates = draw_one_path(r_cyl, l_cyl, wind_angle, fiber_width, starting_offset_angle, num_fibers)
t_size = 1000;
t = linspace(0, l_cyl, t_size);
r = zeros (size(t));
theta_up = zeros(size(t));
theta_down = zeros(size(t));
x_up = t;
x_down = t;

wind_angle = calc_wind_angle(r_cyl, l_cyl, wind_angle, fiber_width)
wind_angle_down = wind_angle;
ota = pi - (2 * wind_angle);
arc_width = (2 * pi / num_fibers)

%populate wind-up vectors
for i = 1 : t_size
    r(i) = r_cyl;
    theta_up(i) = starting_offset_angle + (t(i) * pi) / (2 * r_cyl * cot(wind_angle));
    x_up(i) = t(i);
end

r_coord_ota = zeros(2);
theta_coord_ota = zeros(2);
x_coord_ota = zeros(2);

%initialize path-start marker
r_coord_marker = zeros(2);
theta_coord_marker = zeros(2);
x_coord_marker = zeros(2);

%populate OTA-top vectors
r_coord_ota = [r_cyl r_cyl];
theta_coord_ota = [theta_up(t_size) theta_up(t_size) + ota];
x_coord_ota = [x_up(t_size) x_up(t_size)];

%create path-start marker
r_coord_marker = [r_cyl * (2 / 3) r_cyl];
theta_coord_marker = [theta_up(1) theta_up(1)];
x_coord_marker = [x_up(1) x_up(1)];

%create fiber_width visualizer
r_coord_width = [r_cyl r_cyl];
theta_coord_width = [theta_up(1) - arc_width theta_up(1) + arc_width];
x_coord_width = [x_up(1) x_up(1)];

%populate wind-down vectors
for i = 1 : t_size
    theta_down(i) = (t(i) * pi) / (2 * r_cyl * cot(wind_angle_down)) + theta_coord_ota(2);
    x_down(i) = l_cyl - x_up(i);
end

r_coord_ota_2(1) = r_cyl;
theta_coord_ota_2(1) = theta_down(t_size);

r_coord_ota_2(2) = r_cyl;
theta_coord_ota_2(2) = theta_coord_ota_2(1) + ota;

x_coord_ota_2 = [x_down(t_size) x_down(t_size)];

path_coordinates_r = [r_coord_marker, r_coord_width, r, r_coord_ota, r, r_coord_ota_2];
path_coordinates_theta = [theta_coord_marker, theta_coord_width, theta_up, theta_coord_ota, theta_down, theta_coord_ota_2];
path_coordinates_x = [x_coord_marker, x_coord_width, x_up, x_coord_ota, x_down, x_coord_ota_2];

% first column r's, second column thetas, third column x's
path_coordinates = [path_coordinates_r;path_coordinates_theta;path_coordinates_x];

end

function cartesian_coordinates = graph_path(path_coordinates)
rs = path_coordinates(1, :);
thetas = path_coordinates(2, :);
xs = path_coordinates(3, :);

max(rs)
max(thetas)
max(xs)


for q = 1 : numel(xs)
ys(q) = rs(q) .* cos(thetas(q));
zs(q) = rs(q) .* sin(thetas(q));
end

plot3(xs,ys,zs);
hold on

return

end

%use the initial desired wind angle to determine a wind angle that
%satisfies 'relative primity'
function new_theta = calc_wind_angle(r_cyl, l_cyl, wind_angle, fiber_width)
theta_inc = 0.001;
theta_plus = wind_angle;
theta_minus = theta_plus;

%find the number of fiber widths that can fit around circumference of
%mandrel
fprintf("before: ");
num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width)
%calculate the total fiber step between wind passes at the starting wind
%angle
fiber_step = calc_fiber_step(r_cyl, l_cyl, wind_angle, fiber_width)

temp_num_fibers = num_fibers;
temp_fiber_step = num_fibers;

if (~calc_step_prime(fiber_step, num_fibers))
    fprintf("not rel_prime initially...");
else
    fprintf("rel_prime initially...");
end

%increment the wind angle until relative primity is achieved
while (~calc_step_prime(fiber_step, num_fibers))
    num_fibers = calc_num_fibers(r_cyl, theta_plus, fiber_width);
    fiber_step = calc_fiber_step(r_cyl, l_cyl, theta_plus, fiber_width);
    theta_plus = theta_plus + theta_inc
end

%this ensures the theta-minus half actually gets called
num_fibers = temp_num_fibers;
fiber_step = temp_fiber_step;

%decrement the wind angle until relative primity is achieved
while (~calc_step_prime(fiber_step, num_fibers))
    num_fibers = calc_num_fibers(r_cyl, theta_minus, fiber_width);
    fiber_step = calc_fiber_step(r_cyl, l_cyl, theta_minus, fiber_width);
    theta_minus = theta_minus - theta_inc
end


%pick the relatively prime wind angle closest to the original and return
if (theta_plus - wind_angle) > (wind_angle - theta_minus)
    new_theta = theta_minus;
    fprintf("\ni have chosen theta_minus\n");
else
    new_theta = theta_plus;
    fprintf("\ni have chosen theta_plus\n");
end

fprintf("\nnum_fibers: %d, fiber_step: %d, new_theta: %.2d", num_fibers, fiber_step, (new_theta * 180 / pi));

return
end

%calculate the number of fibers that can fit around the circumference of
%the mandrel at the specified wind angle
function num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width)
width_relative = sec(wind_angle) * fiber_width;
num_fibers = 2 * pi * r_cyl / width_relative;
num_fibers = ceil(num_fibers);
end

function fiber_step = calc_fiber_step(r_cyl, l_cyl, wind_angle, fiber_width)

num_fibers = calc_num_fibers(r_cyl, wind_angle, fiber_width);
%guess an initial OTA using just the maximum possible OTA (known from
%axis/cyl radius)
ota = pi - (2 * wind_angle);
%calculate spiral angle (taken from helix parametric equation)
spiral_angle = l_cyl * pi / (2 * r_cyl * cot(wind_angle));
%calculate fiber steps (k) undertaken every 2OTA + 2spiral cycle.
%total angle spun every cycle
total_angle = mod(2 * (spiral_angle + ota), 2 * pi);
%K = 2(OTA + Spiral) * F/2pi
fiber_step = total_angle * num_fibers / (2 * pi);
%turn this into an integer
fiber_step = floor(fiber_step);

return
end

%determine whether two numbers (our total fiber step and our number of
%fibers) are relatively prime or not
function is_prime = calc_step_prime(fiber_step, num_fibers)
    GCD = gcd(fiber_step, num_fibers);
    if GCD == 1
        is_prime = 1;
    else
        is_prime = 0;
    end
return
end

