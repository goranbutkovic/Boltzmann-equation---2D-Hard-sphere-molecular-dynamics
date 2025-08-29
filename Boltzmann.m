% 2D Hard-sphere molecular dynamics in a box
clear; clc;

% Parameters
N  = 20;        % number of particles
L  = 1.0;       % box side length
r  = 0.03;      % particle radius
dt = 0.005;     % time step
Nt = 3000;      % number of steps

% Initial positions: cluster in one corner
x = 0.2*rand(1,N)*L;
y = 0.2*rand(1,N)*L;

% Initial velocities: random
vx = 1.0*randn(1,N);
vy = 1.0*randn(1,N);

% Colors for particles
colors = lines(N);

% Figure setup
figure('Color','w');
h = scatter(x,y,200,'filled');
xlim([0 L]); ylim([0 L]); axis square;
title('2D Hard-sphere Molecular Dynamics');

for it = 1:Nt
    % --- Update positions ---
    x = x + vx*dt;
    y = y + vy*dt;
    
    % --- Wall collisions ---
    hitLeft   = x < r; hitRight = x > L-r;
    hitBottom = y < r; hitTop   = y > L-r;
    vx(hitLeft | hitRight) = -vx(hitLeft | hitRight);
    vy(hitBottom | hitTop) = -vy(hitBottom | hitTop);
    x(hitLeft)   = r;  x(hitRight) = L-r;
    y(hitBottom) = r;  y(hitTop)   = L-r;
    
    % --- Particle-particle collisions ---
    for i = 1:N
        for j = i+1:N
            dx = x(i)-x(j); dy = y(i)-y(j);
            dist = sqrt(dx^2+dy^2);
            if dist < 2*r
                % Normal unit vector
                nx = dx/dist; ny = dy/dist;
                % Relative velocity
                dvx = vx(i)-vx(j);
                dvy = vy(i)-vy(j);
                vn = dvx*nx + dvy*ny;
                if vn < 0 % approaching
                    % Impulse along normal (equal masses)
                    vx(i) = vx(i) - vn*nx;
                    vy(i) = vy(i) - vn*ny;
                    vx(j) = vx(j) + vn*nx;
                    vy(j) = vy(j) + vn*ny;
                end
                % Push them apart (avoid overlap sticking)
                overlap = 2*r - dist;
                x(i) = x(i) + nx*overlap/2;
                y(i) = y(i) + ny*overlap/2;
                x(j) = x(j) - nx*overlap/2;
                y(j) = y(j) - ny*overlap/2;
            end
        end
    end
    
    % --- Update plot ---
    set(h,'XData',x,'YData',y,'CData',colors);
    drawnow;
end
