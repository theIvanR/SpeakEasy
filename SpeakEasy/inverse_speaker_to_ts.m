% Inverse Problem Implementation
clear; clc; close all;
addpath('functions')


% A: Configure (TS operates in 0.5 - 5 F_res, set accordingly)
% -------------------------
cfg = struct();
cfg.csvFile      = "example_data/impedance_bode_plot.csv";
cfg.cohThreshold = 0.9;
cfg.stride       = 1;
cfg.fMin         = 30;
cfg.fMax         = 300; 

% TS Parameters
p0 = struct();
    % Fixed
    p0.c   = 343;     % m/s
    p0.rho = 1.2;     % kg/m^3
    p0.Sd  = 132e-4;  % m^2
    
    % Electrical
    p0.Re  = 6;       % Ohm
    p0.Le  = 0.29e-3;  % H
    
    % Mechanical
    p0.Bl  = 10;      % N/A
    p0.Rms = 2.0;     % N*s/m
    p0.Mms = 10.4e-3; % kg
    p0.Cms = 77e-6;   % m/N

% Extended Parameters
opts = struct();
    % Enclosure
    opts.box   = false;
    opts.Vb    = 0.2;
    opts.alpha = 1.25;
    
    % Leakage / measurement chain
    opts.leakage = struct();
    opts.leakage.R = 0.61;
    opts.leakage.L = 50e-6;


% B: Run fit
% -------------------------

% 1: Clean data and fit parameters
    df_clean = filter_impedance_csv(cfg.csvFile, cfg.cohThreshold, cfg.stride, cfg.fMin, cfg.fMax);
    result   = fit_speaker_params(df_clean, p0, opts);
    
    clc; % clear the verbose log
    disp('--- Extracted Features from CSV ---');
    disp(result.features);

    disp('--- Computed TS Parameters ---');
    disp(result.p);
    
    disp('--- Auxilliary Derived ---');
    fprintf('Fs (Hz): %.6e\n', result.fs);
    fprintf('Qms      : %.6f\n', result.Qms);
    fprintf('Qes      : %.6f\n', result.Qes);
    fprintf('Qts      : %.6f\n', result.Qts_fit);
    fprintf('Vas (m^3): %.6e\n', result.vas);
    
% Plot Results    
    plot_fit_result(result);
    grid on;


%% Local Functions

% 1: CSV Handling
function df_valid = filter_impedance_csv(filename, coh_threshold, stride, f_min, f_max)
    if nargin < 2 || isempty(coh_threshold), coh_threshold = 0.8; end
    if nargin < 3 || isempty(stride),       stride = 1;          end
    if nargin < 4 || isempty(f_min),         f_min = -inf;        end
    if nargin < 5 || isempty(f_max),         f_max = inf;         end

    df = readtable(filename);
    validate_columns(df, {'frequency_hz','Z_mag_ohm','Z_phase_deg','Coherence'});

    mask = df.Coherence >= coh_threshold;
    mask = mask & df.frequency_hz >= f_min & df.frequency_hz <= f_max;

    df_valid = sortrows(df(mask, :), 'frequency_hz');

    if stride > 1
        df_valid = df_valid(1:stride:end, :);
    end

    fprintf('Total rows: %d\n', height(df));
    fprintf('After coherence + frequency filter: %d\n', sum(mask));
    fprintf('After downsampling (stride=%d): %d\n', stride, height(df_valid));
end

function validate_columns(df, req)
    missing = req(~ismember(req, df.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing required column(s): %s', strjoin(missing, ', '));
    end
end

% 2: Curve Fit
function result = fit_speaker_params(df_clean, p0, opts)

    validate_columns(df_clean, {'frequency_hz','Z_mag_ohm','Z_phase_deg','Coherence'});

    f     = df_clean.frequency_hz(:);
    Zmag  = df_clean.Z_mag_ohm(:);
    Zph   = deg2rad(df_clean.Z_phase_deg(:));
    Zmeas = Zmag .* exp(1j * Zph);
    coh   = df_clean.Coherence(:);

    % Leakage defaults
    if ~isfield(opts, 'leakage') || isempty(opts.leakage)
        opts.leakage = struct('R', 0, 'L', 0);
    else
        if ~isfield(opts.leakage, 'R') || isempty(opts.leakage.R)
            opts.leakage.R = 0;
        end
        if ~isfield(opts.leakage, 'L') || isempty(opts.leakage.L)
            opts.leakage.L = 0;
        end
    end

    % Seed extraction: use only electrical HF info and resonance frequency
    features = estimate_impedance_features(f, Zmag, Zmeas, opts.leakage);

    Re0  = max(features.Re0, 1e-3);
    Le0  = max(features.Le0, 1e-7);
    Mms0 = max(p0.Mms, 1e-6);

    fs0 = features.fs;
    if ~isfinite(fs0) || fs0 <= 0
        fs0 = f(round(numel(f)/2));
    end

    % No Q-based crutch: use user-provided mechanical seeds
    Bl0  = max(p0.Bl, 1e-2);
    Rms0 = max(p0.Rms, 1e-4);

    % Optimization vector
    x0 = log([Re0, Le0, Bl0, Rms0, Mms0]);
    lb = log([1e-2, 1e-5, 1e-2, 1e-4, 1e-4]);
    ub = log([1e2,  1e-2,  1e2,  1e2,   1e-1]);

    fun = @(x) speaker_residual(x, f, Zmeas, coh, p0, opts, fs0);

    lsq_opts = optimoptions('lsqnonlin', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 200, ...
        'MaxIterations', 20, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);

    [x_fit, resnorm, residual, exitflag, output] = lsqnonlin(fun, x0, lb, ub, lsq_opts);

    % Reconstruct fitted parameters
    p_fit = p0;
    p_fit.Re  = exp(x_fit(1));
    p_fit.Le  = exp(x_fit(2));
    p_fit.Bl  = exp(x_fit(3));
    p_fit.Rms = exp(x_fit(4));
    p_fit.Mms = exp(x_fit(5));
    p_fit.Cms = 1 / ((2*pi*fs0)^2 * p_fit.Mms);

    % Model-consistent Q values from fitted parameters
    omega_s = 2*pi*fs0;
    Qms = omega_s * p_fit.Mms / p_fit.Rms;
    Qes = omega_s * p_fit.Mms * p_fit.Re / (p_fit.Bl^2);
    Qts_fit = (Qms * Qes) / (Qms + Qes);
    
    % Find Vas
    Vas = p0.rho * p0.c^2 * (p0.Sd^2) * p_fit.Cms;

    % Find Fs
    fs = 1/(2*pi*sqrt(p_fit.Mms * p_fit.Cms));

    % Make Result Struct
    result = struct();
        result.p         = p_fit;
        
        result.Qms       = Qms;
        result.Qes       = Qes;
        result.Qts_fit   = Qts_fit;
        result.vas       = Vas;
        result.fs        = fs;

        result.resnorm   = resnorm;
        result.residual  = residual;
        result.exitflag  = exitflag;
        result.output    = output;
        result.features  = features;
        result.f         = f;
        result.Zmeas     = Zmeas;
        result.coherence = coh;
    
        result.Z_fit     = simulate_speaker(f, p_fit, opts);

end

% 2.1 Estimate Key Features and helpers
function features = estimate_impedance_features(f, Zmag, Zmeas, leakage, f_max_res)

    if nargin < 4 || isempty(leakage)
        leakage.R = 0;
        leakage.L = 0;
    end
    if nargin < 5 || isempty(f_max_res)
        f_max_res = 5000;
    end

    f = f(:); Zmag = Zmag(:); Zmeas = Zmeas(:);

    w = 2*pi*f;
    Zcorr = Zmeas - leakage.R - 1j .* w .* leakage.L;
    Zmag_corr = abs(Zcorr);

    % smoothing for peak detection
    win = min(11, numel(Zmag_corr));
    Zs = movmean(Zmag_corr, win);

    idx_band = f <= min(f_max_res, max(f)/2);
    if ~any(idx_band)
        error('No data points found below the resonance search limit.');
    end

    f_band  = f(idx_band);
    Zs_band = Zs(idx_band);
    Zr_band = Zmag_corr(idx_band);

    [~, idx_peak] = max(Zs_band);
    fs    = f_band(idx_peak);
    Zpeak = Zr_band(idx_peak);

    % corrected HF estimate
    idx_hf = max(round(0.8 * numel(f)), 1):numel(f);
    omega  = 2*pi*f(idx_hf);
    Zhf    = Zcorr(idx_hf);

    Re_est = max(real(mean(tail(Zhf, min(20,numel(Zhf))))), 1e-3);

    target = Re_est + (Zpeak - Re_est) / sqrt(2);

    f1 = find_crossing_left(f_band, Zr_band, idx_peak, target);
    f2 = find_crossing_right(f_band, Zr_band, idx_peak, target);

    if isfinite(f1) && isfinite(f2) && f2 > f1
        Qts = fs / (f2 - f1);
    else
        Qts = NaN;
    end

    Re0 = median(real(Zhf), 'omitnan');
    Le0 = median(imag(Zhf) ./ max(omega, eps), 'omitnan');

    if ~isfinite(Re0) || Re0 <= 0
        Re0 = max(min(Zmag_corr), 1e-3);
    end
    if ~isfinite(Le0) || Le0 <= 0
        Le0 = 1e-4;
    end

    % --- Output ---
    features = struct( ...
        'fs', fs, ...
        'Zpeak', Zpeak, ...
        'Qts', Qts, ...
        'bandwidth_f1', f1, ...
        'bandwidth_f2', f2, ...
        'Re0', Re0, ...
        'Le0', Le0 ...
    );
end

function y = tail(x, n)
    n = min(n, numel(x));
    y = x(end-n+1:end);
end

function f_cross = find_crossing_left(f, y, idx_peak, target)
    f_cross = NaN;
    for i = idx_peak:-1:2
        if (y(i) >= target && y(i-1) < target) || (y(i) <= target && y(i-1) > target)
            f_cross = interp1([y(i-1), y(i)], [f(i-1), f(i)], target, 'linear', 'extrap');
            return;
        end
    end
end

function f_cross = find_crossing_right(f, y, idx_peak, target)
    f_cross = NaN;
    for i = idx_peak:numel(y)-1
        if (y(i) >= target && y(i+1) < target) || (y(i) <= target && y(i+1) > target)
            f_cross = interp1([y(i), y(i+1)], [f(i), f(i+1)], target, 'linear', 'extrap');
            return;
        end
    end
end

% 2.2 Compute residual (loss) 
function r = speaker_residual(x, f, Zmeas, coh, p_fixed, opts, fs0)
    p = p_fixed;
    p.Re  = exp(x(1));
    p.Le  = exp(x(2));
    p.Bl  = exp(x(3));
    p.Rms = exp(x(4));
    p.Mms = exp(x(5));
    p.Cms = 1 / ((2*pi*fs0)^2 * p.Mms);

    Zmodel = simulate_speaker(f, p, opts);
    dZ = Zmodel(:) - Zmeas(:);

    scale = max(abs(Zmeas(:)), 0.25 * p.Re);
    scale = max(scale, 1e-6);
    w = max(coh(:), 0.05);

    r = [real(dZ)./scale .* w; imag(dZ)./scale .* w];
end

% 3: Plot Results
function plot_fit_result(result)
    f = result.f;
    Zm = result.Zmeas;
    Zf = result.Z_fit;

    figure;

    subplot(2,1,1)
    plot(f, abs(Zm)); hold on;
    plot(f, abs(Zf), '--');
    ylabel('|Z| (Ohm)');
    legend('Measured','Fitted');
    grid on;

    subplot(2,1,2)
    plot(f, unwrap(angle(Zm))*180/pi); hold on;
    plot(f, unwrap(angle(Zf))*180/pi, '--');
    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    legend('Measured','Fitted');
    grid on;
end

