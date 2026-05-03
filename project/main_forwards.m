%% Forward Problem Implementation
clear; close all; clc;
addpath(genpath(pwd));

% Frequency grid
f = logspace(1, 4, 100);

% Transducer intrinsic parameters
p = struct();
p.Re  = 4.7;
p.Le  = 0.34e-3;
p.Sd  = 0.0132;
p.Bl  = 9.46;
p.Rms = 0.93;
p.Mms = 0.019;
p.Cms = 294e-6;

% Geometry (source manifold)
Nrad   = 32;
Mtheta = 32;
geom = build_circle_mesh(p.Sd, Nrad, Mtheta);

% Environment / field
env = struct();
env.c = 343;
env.rho = 1.2000;
env.r_ref = 1.0;
env.off_axis = 0;
env.azimuth = 0;

% enclosure (optional boundary condition)
env.enc = struct( ...
    'type','sealed', ...
    'Vb', 20e-3, ...
    'Rloss', 0, ...
    'Rleak', 0);

% drive (normalized)
V_in = sqrt(p.Re);

% simulate
results = run_electroacoustic_engine(f, p, geom, env);

% post
metrics = speaker_power_metrics(results, V_in);

make_plots(results, metrics, V_in, env);
%plot_Ze(results);