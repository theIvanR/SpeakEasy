% Forwards Problem Implementation
clear; close all; clc;
addpath("functions\");

% Frequency grid
f = logspace(1, 4, 2000);

% Speaker parameters
p = struct();
    % Fixed
    p.c   = 343;     % m/s
    p.rho = 1.2;     % kg/m^3
    p.Sd  = 132e-4;  % m^2
    
    % Electrical
    p.Re  = 6;       % Ohm
    p.Le  = 0.29e-3;  % H
    
    % Mechanical
    p.Bl  = 10;      % N/A
    p.Rms = 2.0;     % N*s/m
    p.Mms = 10.4e-3; % kg
    p.Cms = 77e-6;   % m/N

% Auxiliary / model configuration
opts = struct();
    opts.Vb             = 0.02;      % sealed-box volume [m^3]
    opts.alpha          = 1;         % sealed-box damping factor
    opts.radiationModel = 'bem';     % 'constant' | 'bessel' | 'bem'
    opts.r_ref          = 1.0;       % reference distance [m]
    opts.plot           = 'reduced'; % 'none' | 'reduced' | 'full'

results = simulate_speaker_full(f, p, opts);