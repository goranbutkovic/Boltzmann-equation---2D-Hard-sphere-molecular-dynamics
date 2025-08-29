% 2D Hard-sphere molecular dynamics in a box (smooth GIF)
clear; clc;

% --- Parameters ---
N  = 20;        % number of particles
L  = 1.0;       % box side length
r  = 0.03;      % particle radius
dt = 0.005;     % time step
Nt = 3000;      % total simulation steps
frameStep = 3;  % save every 3rd frame for GIF smoothness

% Initial positions: cluster in one corner
x = 0.2*rand(1,N)*L/4;
y = 0.2*rand(1,N)*L/4;

% Initial velocities: random
vx = 1.0*randn(1,N);
vy = 1.0*randn(1,N);

% Colors for particles
colors = lines(N);

% --- Figure setup ---
figure('Color','w');
h = scatter(x,y,200,'filled');
xlim([0 L]); ylim([0 L]); axis square;
title('2D Hard-sphere Molecular Dynamics');

% --- GIF setup ---
gifFileName = 'HardSphereMD.gif';
delayTime = dt*frameStep; % delay between frames in seconds

% --- Main loop ---
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
                nx = dx/dist; ny = dy/dist;
                dvx = vx(i)-vx(j); dvy = vy(i)-vy(j);
                vn = dvx*nx + dvy*ny;
                if vn < 0
                    vx(i) = vx(i) - vn*nx;
                    vy(i) = vy(i) - vn*ny;
                    vx(j) = vx(j) + vn*nx;
                    vy(j) = vy(j) + vn*ny;
                end
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
    
    % --- Save to GIF every frameStep steps ---
    if mod(it, frameStep) == 0
        frame = getframe(gcf);
        im = frame2im(frame);
        [A,map] = rgb2ind(im,256);
        if it == frameStep
            imwrite(A,map,gifFileName,'gif','LoopCount',Inf,'DelayTime',delayTime);
        else
            imwrite(A,map,gifFileName,'gif','WriteMode','append','DelayTime',delayTime);
        end
    end
end

disp(['GIF saved as ', gifFileName]);
