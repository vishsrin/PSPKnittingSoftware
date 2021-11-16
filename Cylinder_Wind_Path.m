

function draw_full_path(r_cyl, l_cyl, wind_angle, r_axis, fiber_width)
wind_angle = deg2rad(wind_angle);


end


function path_coordinates = draw_one_path(r_cyl, l_cyl, wind_angle, r_axis, fiber_width, starting_offset_angle)
t_size = 1000;
t = linspace(0, l_cyl, t_size);
r = zeros (size(t));
theta_up = zeros(size(t));
theta_down = zeros(size(t));
x_up = t;
x_down = t;

wind_angle = calc_wind_angle(r_cyl, l_cyl, wind_angle, fiber_width);
wind_angle_down = pi - wind_angle;
ota = pi - (2 * wind_angle);

%populate wind-up vectors
for i = 1 : t_size
    r(i) = r_cyl;
    theta_up(i) = starting_offset_angle + (t(i) * pi) / (2 * r_cyl * cot(wind_angle));
    x_up(i) = t(i);
end

r_coord_ota = zeros(2);
theta_coord_ota = zeros(2);
x_coord_ota = zeros(2);

%populate OTA-top vectors
r_coord_ota(1) = r_cyl;
theta_coord_ota(1) = theta_up(t_size);
x_coord_ota(1) = x_up(t_size);

r_coord_ota(2) = r_cyl;
theta_coord_ota(2) = thet_coord_ota(1) + ota;
x_coord_ota = x_coord_ota(1);

%populate wind-down vectors
for i = 1 : t_size
    theta_down(i) = (t(i) * pi) / (2 * r_cyl * cot(wind_angle_down)) + theta_coord_ota(2);
    x_down(i) = l_cyl - x_up(i);
end

r_coord_ota_2(1) = r_cyl;
theta_coord_ota_2(1) = theta_down(t_size);
x_coord_ota_2(1) = x_down(t_size);

r_coord_ota_2(2) = r_cyl;
theta_coord_ota_2(2) = theta_coord_ota_2(1) + ota;
x_coord_ota_2(2) = x_coord_ota(2);

path_coordinates_r = [r, r_coord_ota, r, r_coord_ota_2];
path_coordinates_theta = [theta_up, theta_coord_ota, theta_down, theta_coord_ota_2];
path_coordinates_x = [x_up, x_coord_ota_1, x_down, x_coord_ota_2];

% first column r's, second column thetas, third column x's
path_coordinates = [path_coordinates_r;path_coordinates_theta;path_coordinates_x];


end
