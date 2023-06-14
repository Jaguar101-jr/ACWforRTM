
%% A simple Reverse-Time Migration (RTM) Acoustic Wavefield Extrapolator, 
%% A 2nd order Finite Difference Method program to solve/propagate (i.e., forward modelling) the 2D constant-density,
%% wave equation on a Cartesian mesh, wth 2D for recievers' signals, 3D plots of wavefield and a movie; Acoustic-Cartesian-Waves ACW for RTM. 
%% Please cite this code as: Hussein Abduelhaleim Hussein Muhammed (2022), Least-
%% Squares ReverseTime Migration in Pseudodepth Domain and Its Application.
%% China University of Petroleum (East China), Master Thesis,
%% School of Geosciences, Dept. of Geophysics. Library press©.

clc;
clear;
close all;

%% Prepare the movie file.
    vidObj = VideoWriter('ACWRTM.avi');
    open(vidObj);

%% Grid parameters
nx = 50;                % Number of grid points in x-direction
ny = 50;                % Number of grid points in y-direction
nt = 500;               % Number of time steps
dx = 20;                % Spatial step size in x-direction (m)
dy = 20;                % Spatial step size in y-direction (m)
dt = 0.001;             % Time step-size (in sec.) or sampling rate
c0 = 1500;              % Velocity of the medium (m/s)
srcx = 25;              % Source x-coordinate
srcy = 25;              % Source y-coordinate
recx = 30:5:40;         % Receiver x-coordinates
recy = 10:5:30;         % Receiver y-coordinates

%% Initialization
c = c0 * ones(nx, ny);  % Velocity model
u = zeros(nx, ny);      % Wavefield at current time step
u1 = zeros(nx, ny);     % Wavefield at previous time step
u2 = zeros(nx, ny);     % Wavefield at two time steps ago

% Source time function (Ricker source wavelet)
f0 = 30;                % Dominant frequency (Hz)
t0 = 1 / f0;            % Delay (s)
t = (0:nt-1) * dt;      % Time vector


%% Time-stepping solution loop

for n = 1:nt
    % Update source wavelet position
    srcx = srcx + 1;  % Example: Move source in positive x-direction

    % Ensure source stays within the grid
    if srcx > nx
        srcx = nx;
    end

    % Compute source term at current time step
    src = (1 - 2 * pi^2 * f0^2 * (t(n) - t0)^2) * exp(-pi^2 * f0^2 * (t(n) - t0)^2);

    % Inject source wavelet
    u(srcx, srcy) = u(srcx, srcy) + src * dt^2 / (dx * dy);

    % Apply the finite-difference stencil coeff. loop
    u(2:end-1, 2:end-1) = (2 - 2 * (c(2:end-1, 2:end-1) * dt / dx)^2 - 2 * (c(2:end-1, 2:end-1) * dt / dy)^2) .* u1(2:end-1, 2:end-1) ...
        - u2(2:end-1, 2:end-1) ...
        + (c(2:end-1, 2:end-1) * dt / dx)^2 .* (u1(3:end, 2:end-1) + u1(1:end-2, 2:end-1)) ...
        + (c(2:end-1, 2:end-1) * dt / dy)^2 .* (u1(2:end-1, 3:end) + u1(2:end-1, 1:end-2));

    % Record the wavefield at receivers
    for i = 1:length(recx)
        u_rec(n, i) = u(recx(i), recy(i));
    end

    % Update wavefields for next time step
    u2 = u1;
    u1 = u;
   
    % Apply thresholding to only show positive values of the wavefield u
    % u(u < 0) = 0;  % Set any negative values to zero

    % Plot wavefield at each time step
    figure(1);
    surf(u);
    title(sprintf('Wavefield at Time Step %d', n));
    xlabel('x');
    ylabel('y');
    zlabel('Amplitude');
    shading interp;
    colormap jet;
    drawnow;

    % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);

end

% Plot recorded waveforms at receivers
figure(2);
t_rec = (0:nt-1) * dt;

for i = 1:length(recx)
    subplot(length(recx), 1, i);
    plot(t_rec, u_rec(:, i));
    title(sprintf('Recorded Waveform at Receiver (%d, %d)', recx(i), recy(i)));
    xlabel('Time (s)');
    ylabel('Amplitude');
end

%% Close the file.
close(vidObj);
