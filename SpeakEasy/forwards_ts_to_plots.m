% Forwards Problem Implementation
clear; close all; clc;
addpath("functions\");

% Frequency grid
f = logspace(1, 4, 2000);

% Speaker parameters
p = struct();
    % Fixed
    p.c   = 343;       % m/s
    p.rho = 1.2000;    % kg/m^3
    p.Sd  = 0.0132;    % m^2
    
    % Electrical
    p.Re  = 6.1007;    % Ohm
    p.Le  = 4.2008e-04; % H
    
    % Mechanical
    p.Bl  = 5.0213;    % N/A
    p.Rms = 3.2688;    % N*s/m
    p.Mms = 0.0287;    % kg
    p.Cms = 1.5811e-04; % m/N

% Auxiliary / model configuration
opts = struct();
    opts.Vb             = 0.02;      % sealed-box volume [m^3]
    opts.alpha          = 0.1;         % sealed-box damping factor
    opts.radiationModel = 'bem';     % 'constant' | 'bessel' | 'bem'
    opts.r_ref          = 1.0;       % reference distance [m]
    opts.plot           = 'reduced'; % 'none' | 'reduced' | 'full'

results = simulate_speaker_full(f, p, opts);