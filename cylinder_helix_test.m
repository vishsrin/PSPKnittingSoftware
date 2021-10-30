%draw_cyl_with_helix(2,15,69)


function draw_cyl_with_helix(r_cyl,l_cyl,wind_angle)
alpha = deg2rad(wind_angle); % wind angle, radians
t_size = 1000;

%phi0 = 0; % starting phi
%zeta = 0.5 * h / r;
%theta0 = atan(1 / zeta); % starting theta
%theta_span = theta0:0.01:(pi / 2); % theta span vectors

t = linspace(0,l_cyl,t_size);
r = zeros(size(t));
theta = zeros(size(t));
z = zeros(size(t));

for i = 1:length(t)
    r(i) = r_cyl;
    theta(i) = (t(i) * pi) / (2 * r_cyl * cot(alpha));
    x(i) = t(i);
end

y = r .* cos(theta);
z = r .* sin(theta);

plot3(x,y,z);
axis equal
hold on
[y_cyl,z_cyl,x_cyl] = cylinder(r_cyl);
x_cyl = x_cyl * l_cyl;
mesh(x_cyl,y_cyl,z_cyl);
end

function ota = calc_ota(r_cyl, l_cyl, r_boss, wind_angle, fiber_width)
width_relative = acos(wind_angle) * fiber_width; %step 1: Calculate W_r
num_fibers = 2 * pi * r_cyl / width_relative; %step 2: Calculate F
num_fibers = ceil(num_fibers); %step 3: ceil F
ota = 2 * acos(r_boss / r_cyl); %step 4: Initial angle guess at max possible
fiber_step = floor((num_fibers / pi) * mod(((l_cyl * pi) / (2 * r_cyl * cot(wind_angle))), 2 * pi) + ota); %step 5: generate K value from ota guess
fiber_step = find_fiber_step(fiber_step, num_fibers);
ota = (fiber_step / num_fibers) * 2 * pi;
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
        fiber_step_inc = fiber_step_inc + 1
        GCD_inc = gcd(fiber_step_inc, num_fibers);
        if GCD_inc == 1
            fiber_step = fiber_step_inc;
            return
        end
        
        fiber_step_dec = fiber_step_dec - 1
        GCD_dec = gcd(fiber_step_dec, num_fibers);
        if GCD_dec == 1
            fiber_step = fiber_step_dec;
            return
        end
    end
end