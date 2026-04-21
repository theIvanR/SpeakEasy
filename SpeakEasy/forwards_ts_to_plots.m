% Forwards Problem Implementation
clear; close all; clc;
addpath("functions\");

% Frequency grid
f = logspace(1, 4, 2000);

% Speaker parameters
p = struct( ...
    'c',    343,      ... % m/s
    'rho',  1.2,      ... % kg/m^3
    'Sd',   132e-4,   ... % m^2
    'Re',   6,        ... % Ohm
    'Le',   0.29e-3,  ... % H
    'Bl',   10,       ... % N/A
    'Rms',  2.0,      ... % N*s/m
    'Mms',  10.4e-3,  ... % kg
    'Cms',  77e-6    ... % m/N
);

% Auxilliary Parameters
opts = struct( ...
    'cfg',            'single', ... % 'single' | 'series' | 'parallel' | 'series_parallel' , (series = "ideal isobaric")
    'Vb',             0.02,      ... % sealed-box volume [m^3]
    'alpha',          1,     ... % sealed-box damping factor
    'radiationModel', 'bem',    ... % 'constant' | 'bessel' | 'bem'
    'r_ref',          1.0,      ... % reference distance for pressure proxy [m]
    'plot',           'reduced' );    % 'none' | 'reduced' | 'full'

results = simulate_speaker_full(f, p, opts);
