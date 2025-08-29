%% Sticky Hard-Sphere / MD with Temperature Control
clear; clc; close all;

% --- Parameters ---
N = 80;           % number of particles
L = 10;           % box size
r = 0.15;         % particle radius (repulsive core)
epsilon = 1.0;    % attractive potential depth
dt = 0.003;       % time step
steps_per_frame = 10;
kB = 1.0;         % Boltzmann constant
targetT = 1.0;    % initial temperature

% --- Initial positions and velocities ---
x = 0.3*rand(1,N)*L;
y = 0.3*rand(1,N)*L;
vx = randn(1,N);
vy = randn(1,N);

% remove net momentum
vx = vx - mean(vx);
vy = vy - mean(vy);

% --- figure ---
fig = figure('Color','w');
h = scatter(x,y,80,'filled');
axis([0 L 0 L]); axis square; box on;
colormap(parula);
title('Sticky Hard-Sphere MD: Drag T to see droplet/gas');

% slider for temperature
sld = uicontrol('Style','slider','Units','normalized',...
    'Position',[0.18 0.02 0.6 0.04],'Min',0.1,'Max',5,'Value',targetT);
txt = uicontrol('Style','text','Units','normalized','Position',[0.81 0.02 0.16 0.05],...
    'String',sprintf('T = %.2f',targetT), 'FontSize',10);

% --- Simulation loop ---
while ishandle(fig)
    targetT = sld.Value;
    txt.String = sprintf('T = %.2f',targetT);
    
    for step = 1:steps_per_frame
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
        
        % --- Particle-particle collisions & attraction ---
        for i = 1:N-1
            for j = i+1:N
                dx = x(i)-x(j);
                dy = y(i)-y(j);
                dist = sqrt(dx^2+dy^2);
                nx = dx/dist; ny = dy/dist;
                dvx = vx(i)-vx(j);
                dvy = vy(i)-vy(j);
                
                % --- Hard-core collision ---
                if dist < 2*r
                    vn = dvx*nx + dvy*ny;
                    if vn < 0
                        vx(i) = vx(i) - vn*nx;
                        vy(i) = vy(i) - vn*ny;
                        vx(j) = vx(j) + vn*nx;
                        vy(j) = vy(j) + vn*ny;
                    end
                    % push apart
                    overlap = 2*r - dist;
                    x(i) = x(i) + nx*overlap/2;
                    y(i) = y(i) + ny*overlap/2;
                    x(j) = x(j) - nx*overlap/2;
                    y(j) = y(j) - ny*overlap/2;
                end
                
                % --- Short-range attraction ---
                if dist > 2*r && dist < 4*r
                    % simple linear attractive impulse
                    fmag = epsilon*(6 - (dist-2*r)/(2*r)); % linearly decreases to zero
                    vx(i) = vx(i) - fmag*nx*dt;
                    vy(i) = vy(i) - fmag*ny*dt;
                    vx(j) = vx(j) + fmag*nx*dt;
                    vy(j) = vy(j) + fmag*ny*dt;
                end
            end
        end
        
        % --- Thermostat: rescale velocities toward target T ---
        KE = 0.5*sum(vx.^2 + vy.^2);
        T_inst = KE/(N*kB);
        lambda = sqrt(targetT/T_inst);
        vx = vx * lambda;
        vy = vy * lambda;
    end
    
    % --- Update plot ---
    set(h,'XData',x,'YData',y);
    drawnow;
end