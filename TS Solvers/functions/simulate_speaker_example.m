%% Minimal code to simulate speaker

clear; clc; close all; 
addpath("functions\")

% Initialize Speaker
f = logspace(1,4,2000);

p = struct( ...
    'Re',   6,        ... % Ohm
    'Le',   0.29e-3,    ... % H
    'Bl',   10,       ... % N/A
    'Rms',  2.0,      ... % N*s/m
    'Mms',  10.4e-3,    ... % kg
    'Cms',  77e-6,     ... % m/N
    'Sd',   132e-4,     ... % m^2
    'rho',  1.2,        ... % kg/m^3
    'c',    343         ... % m/s
);

opts = struct( ...
    'cfg',   'single',  ... % 'single' | 'series' | 'parallel' | 'series_parallel'
    'box',   false,      ... % true = sealed box enabled
    'Vb',    0.2,       ... % m^3 (box volume)
    'alpha', 1.25       ... % dimensionless damping factor
);


%Plot speaker
Ze = simulate_speaker(f, p, opts);
figure;
  
    subplot(2,1,1)
    semilogx(f, abs(Ze));
    grid on;
    ylabel('|Z_e| (\Omega)');
    title('Impedance Magnitude');
    
    subplot(2,1,2)
    semilogx(f, angle(Ze) * 180/pi);
    grid on;
    xlabel('Frequency (Hz)');
    ylabel('Phase (deg)');
    title('Impedance Phase');
