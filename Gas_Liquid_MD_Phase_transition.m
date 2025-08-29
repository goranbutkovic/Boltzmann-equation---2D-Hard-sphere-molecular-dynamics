%% Sticky Hard-Sphere MD with Temperature Control and GIF (Frozen droplet)
clear; clc; close all;

% --- Parameters ---
N = 280;               % number of particles
L = 10;                % box size
r = 0.15;              % particle radius (repulsive core)
epsilon = 1.0;         % attractive potential depth
dt = 0.003;            % time step
steps_per_frame = 10;
kB = 1.0;              % Boltzmann constant
targetT = 5.0;         % initial temperature
targetT_min = 0.01;    % minimum temperature to stop cooling
T_freeze = 0.05;       % below this, freeze droplet
extra_time = 5;         % seconds to continue after reaching min T
gif_filename = fullfile(pwd,'phase_transition.gif'); % save in MATLAB folder

% --- Initial positions: cluster in center ---
cluster_size = 2.0;
x = (L/2 - cluster_size/2) + cluster_size*rand(1,N);
y = (L/2 - cluster_size/2) + cluster_size*rand(1,N);

% --- Initial velocities ---
vx = 0.5*randn(1,N);
vy = 0.5*randn(1,N);

% remove net momentum
vx = vx - mean(vx);
vy = vy - mean(vy);

% --- Figure setup ---
fig = figure('Color','w');
h = scatter(x,y,80,'filled');
axis([0 L 0 L]); axis square; box on;
colormap(parula);
title('Gas–Liquid Phase Transition Using MD');

% Temperature display
txt = uicontrol('Style','text','Units','normalized','Position',[0.81 0.02 0.16 0.05],...
    'String',sprintf('T = %.2f',targetT), 'FontSize',10);

frame_count = 0;

%% --- Simulation loop: cooling ---
while ishandle(fig) && targetT > targetT_min
    targetT = targetT - 0.005;   % decrease temperature
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
                dx = x(i)-x(j); dy = y(i)-y(j);
                dist = sqrt(dx^2 + dy^2);
                nx = dx/dist; ny = dy/dist;
                dvx = vx(i)-vx(j); dvy = vy(i)-vy(j);
                
                % Hard-core collision
                if dist < 2*r
                    vn = dvx*nx + dvy*ny;
                    if vn < 0
                        vx(i) = vx(i) - vn*nx; vy(i) = vy(i) - vn*ny;
                        vx(j) = vx(j) + vn*nx; vy(j) = vy(j) + vn*ny;
                    end
                    overlap = 2*r - dist;
                    x(i) = x(i) + nx*overlap/2; y(i) = y(i) + ny*overlap/2;
                    x(j) = x(j) - nx*overlap/2; y(j) = y(j) - ny*overlap/2;
                end
                
                % Short-range attraction
                if dist > 2*r && dist < 4*r
                    fmag = epsilon*(5 - (dist-2*r)/(2*r)); % stronger attraction
                    vx(i) = vx(i) - fmag*nx*dt; vy(i) = vy(i) - fmag*ny*dt;
                    vx(j) = vx(j) + fmag*nx*dt; vy(j) = vy(j) + fmag*ny*dt;
                end
            end
        end
        
        % --- Thermostat (stop below T_freeze) ---
        if targetT > T_freeze
            KE = 0.5*sum(vx.^2 + vy.^2);
            T_inst = KE/(N*kB);
            lambda = sqrt(max(targetT/T_inst,0));
            vx = vx * lambda;
            vy = vy * lambda;
        else
            % small damping to freeze cluster
            vx = vx * 0.995;
            vy = vy * 0.995;
        end
    end
    
    % --- Update plot ---
    set(h,'XData',x,'YData',y);
    drawnow;
    
    % --- Capture frame for GIF ---
    frame_count = frame_count + 1;
    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if frame_count == 1
        imwrite(A,map,gif_filename,'gif','LoopCount',Inf,'DelayTime',0.05);
    else
        imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',0.05);
    end
end

%% --- Extra frames at minimum temperature (~5 seconds) ---
extra_frames = round(extra_time / (dt*steps_per_frame));
for ef = 1:extra_frames
    % --- Update positions ---
    x = x + vx*dt*steps_per_frame;
    y = y + vy*dt*steps_per_frame;
    
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
            dx = x(i)-x(j); dy = y(i)-y(j);
            dist = sqrt(dx^2 + dy^2);
            nx = dx/dist; ny = dy/dist;
            dvx = vx(i)-vx(j); dvy = vy(i)-vy(j);
            
            % Hard-core collision
            if dist < 2*r
                vn = dvx*nx + dvy*ny;
                if vn < 0
                    vx(i) = vx(i) - vn*nx; vy(i) = vy(i) - vn*ny;
                    vx(j) = vx(j) + vn*nx; vy(j) = vy(j) + vn*ny;
                end
                overlap = 2*r - dist;
                x(i) = x(i) + nx*overlap/2; y(i) = y(i) + ny*overlap/2;
                x(j) = x(j) - nx*overlap/2; y(j) = y(j) - ny*overlap/2;
            end
            
            % Short-range attraction
            if dist > 2*r && dist < 4*r
                fmag = epsilon*(5 - (dist-2*r)/(2*r));
                vx(i) = vx(i) - fmag*nx*dt; vy(i) = vy(i) - fmag*ny*dt;
                vx(j) = vx(j) + fmag*nx*dt; vy(j) = vy(j) + fmag*ny*dt;
            end
        end
    end
    
    % --- Small damping at frozen state ---
    vx = vx * 0.995;
    vy = vy * 0.995;
    
    % --- Update plot & GIF ---
    set(h,'XData',x,'YData',y);
    drawnow;
    
    frame_count = frame_count + 1;
    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,gif_filename,'gif','WriteMode','append','DelayTime',0.05);
end

disp(['GIF saved as ', gif_filename]);
