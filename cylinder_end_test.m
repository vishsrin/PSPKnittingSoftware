r_cyl = 2;
r_axis = 0.25;
l_cyl = 15;
f_width = 0.1;
alpha = 79;

ota = calc_ota(r_cyl, r_axis, f_width, alpha);
draw_cyl_end(r_cyl, r_axis, ota, 0)

function ota = calc_ota(r_cyl, r_axis, f_width, alpha)
ota_max = 2 * acos(r_axis / r_cyl);
fiber_count = floor(2 * pi * r_cyl * cos(deg2rad(alpha)) / f_width);
fiber_step = floor(fiber_count * (deg2rad(ota_max) + (pi / (2 * r_cyl * cot(deg2rad(alpha))))) / pi);
i = 0.1;
ota = ota_max;
while (gcd(fiber_step, fiber_count) ~= 1)
    ota = ota - i;
    fiber_step = floor(fiber_count * (deg2rad(ota) + (pi / (2 * r_cyl * cot(deg2rad(alpha))))) / pi);
end
ota = rad2deg(ota);
end

function draw_cyl_end(r_cyl, r_axis, ota, alpha1)
t_size = 100;
ota = deg2rad(ota);
theta = linspace(0, 2 * pi, t_size);
t = linspace(0, 1, t_size);

x1 = r_cyl * cos(alpha1);
y1 = r_cyl * sin(alpha1);
x2 = r_cyl * cos(alpha1 + ota);
y2 = r_cyl * sin(alpha1 + ota);
dist_x = x2 - x1;
dist_y = y2 - y1;
x = zeros(size(t));
y = zeros(size(t));

for i = 1:0.01:size(t)
    x(i) = x1 + (i * dist_x);
    y(i) = y1 + (i * dist_y);
end

plot(x, y)

axis equal
hold on
des_x = r_cyl .* cos(theta);
des_y = r_cyl .* sin(theta);
plot(des_x, des_y)
hold on
des_x = r_axis .* cos(theta);
des_y = r_axis .* sin(theta);
plot(des_x, des_y)
xlabel("x");
ylabel("y");
end