% Inverse Problem Implementation
clear; clc; close all;
addpath('functions')


% A: Configure
% -------------------------
cfg = struct();
cfg.csvFile      = "impedance_bode_plot.csv";
cfg.cohThreshold = 0.9;
cfg.stride       = 1;
cfg.fMin         = 30;
cfg.fMax         = 300;

% Initial guess of parameters
p0 = struct( ...
    'c',    343,      ... % m/s
    'rho',  1.2,      ... % kg/m^3
    'Sd',   132e-4,   ... % m^2
    'Re',   6,        ... % Ohm
    'Le',   0.29e-3,  ... % H
    'Bl',   10,       ... % N/A
    'Rms',  2.0,      ... % N*s/m
    'Mms',  10.4e-3,  ... % kg
    'Cms',  77e-6    ... % m/N
); %n R,L_wire is also contained here (and in csv!) ~1.7 Ohm and 50uH 

opts = struct('cfg','single','box',false);


% B: Run fit
% -------------------------

% 1: Clean data and fit parameters
    df_clean = filter_impedance_csv(cfg.csvFile, cfg.cohThreshold, cfg.stride, cfg.fMin, cfg.fMax);
    result   = fit_speaker_params(df_clean, p0, opts);
    
    disp(result.p);
    disp(result.features);
    Vas = result.p.rho * result.p.c^2 * result.p.Sd^2 * result.p.Cms;
    disp(Vas);
    
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

    features = estimate_impedance_features(f, Zmag, Zmeas);

    Re0 = max(features.Re0, 1e-3);
    Le0 = max(features.Le0, 1e-7);
    Mms0 = max(p0.Mms, 1e-6);

    fs0 = features.fs;
    if ~isfinite(fs0) || fs0 <= 0
        fs0 = f(round(numel(f)/2));
    end

    if isfinite(features.Qts) && features.Qts > 0
        Qts0 = features.Qts;
    else
        Qts0 = 0.4;
    end

    alpha = max(features.Zpeak / max(Re0, eps) - 1, 0.2);
    Qms0  = max(Qts0 * (1 + alpha), 0.1);
    Qes0  = max(Qms0 / alpha, 0.1);

    Rms0 = max((2*pi*fs0*Mms0) / Qms0, 1e-4);
    Bl0  = max(sqrt((2*pi*fs0*Mms0*Re0) / Qes0), 1e-2);

    % x = [log(Re), log(Le), log(Bl), log(Rms), log(Mms)]
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

    p_fit = p0;
    p_fit.Re  = exp(x_fit(1));
    p_fit.Le  = exp(x_fit(2));
    p_fit.Bl  = exp(x_fit(3));
    p_fit.Rms = exp(x_fit(4));
    p_fit.Mms = exp(x_fit(5));
    p_fit.Cms = 1 / ((2*pi*fs0)^2 * p_fit.Mms);

    result = struct();
    result.p         = p_fit;
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
function features = estimate_impedance_features(f, Zmag, Zmeas, f_max_res)

    if nargin < 4 || isempty(f_max_res)
        f_max_res = 5000;
    end

    f = f(:); Zmag = Zmag(:); Zmeas = Zmeas(:);

    if numel(f) ~= numel(Zmag) || numel(f) ~= numel(Zmeas)
        error('Inputs must have the same length.');
    end

    % --- Light smoothing ONLY for peak detection ---
    win = min(11, numel(Zmag));
    Zs = movmean(Zmag, win);

    % --- Resonance search band ---
    idx_band = f <= min(f_max_res, max(f)/2);
    if ~any(idx_band)
        error('No data points found below the resonance search limit.');
    end

    f_band   = f(idx_band);
    Zs_band  = Zs(idx_band);     % smoothed (for detection)
    Zr_band  = Zmag(idx_band);   % raw (for measurement)

    % --- Peak detection (robust) ---
    [~, idx_peak] = max(Zs_band);
    fs    = f_band(idx_peak);
    Zpeak = Zr_band(idx_peak);   % TRUE peak from raw data

    % --- Re estimate (HF region) ---
    Re_est = max(real(mean(tail(Zmeas, min(20,numel(Zmeas))))), 1e-3);

    % --- Q estimation (use RAW data!) ---
    target = Re_est + (Zpeak - Re_est) / sqrt(2);

    f1 = find_crossing_left(f_band, Zr_band, idx_peak, target);
    f2 = find_crossing_right(f_band, Zr_band, idx_peak, target);

    if isfinite(f1) && isfinite(f2) && f2 > f1
        Qts = fs / (f2 - f1);
    else
        Qts = NaN;
    end

    % --- High-frequency electrical seed ---
    idx_hf = max(round(0.8 * numel(f)), 1):numel(f);
    omega  = 2*pi*f(idx_hf);
    Zhf    = Zmeas(idx_hf);

    Re0 = median(real(Zhf), 'omitnan');
    Le0 = median(imag(Zhf) ./ max(omega, eps), 'omitnan');

    % --- Fallbacks ---
    if ~isfinite(Re0) || Re0 <= 0
        Re0 = max(min(Zmag), 1e-3);
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

