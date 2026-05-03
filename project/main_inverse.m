%% Inverse Problem Implementation (with confidence intervals)
clear; clc; close all;
addpath(genpath(pwd));

% ============================================================
% A) Configuration
% ============================================================
cfg = struct();
cfg.csvFile      = "/auxilliary/impedance_bode_plot.csv";
cfg.cohThreshold = 0.90;
cfg.stride       = 1;
cfg.bandLow      = 0.5;
cfg.bandHigh     = 5.0;

% --- SMART fitter controls ---
cfg.regLambda    = 1.0;      % overall MAP prior strength
cfg.freqFocusStrength = 0.5;  % emphasise resonance region
cfg.freqFocusWidth    = 0.7;

% --- Resonance anchor ---
cfg.fsAnchorLambda = 100;     % strong soft anchor to measured fs0

% --- Prior shaping ---
cfg.priorSigmaMin = 0.05;
cfg.priorSigmaMax  = 10.0;
cfg.priorSharpness = 1.5;
cfg.priorUsePass1SVD = true;
cfg.priorCenterMode = 'pass1';   % 'x0' or 'pass1'
cfg.includePriorInFisher = true;

% --- Solver / FD Jacobian ---
cfg.maxIter = 100;
cfg.maxEval = 1200;
cfg.fdStep  = 1e-5;
cfg.verbose = true;

% ============================================================
% B) Initial transducer parameters [GUESS AS CLOSE AS POSSIBLE]
% ============================================================
p0 = struct();

% Fixed
p0.c   = 343;
p0.rho = 1.2;
p0.Sd  = 140e-4;

% Initial guesses
p0.Re  = 6.0;
p0.Le  = 0.4e-3;
p0.Bl  = 8.0;
p0.Rms = 1.0;
p0.Mms = 25e-3;
p0.Cms = 300e-6;

% ============================================================
% C) Geometry / field / enclosure
% ============================================================
Nrad   = 32;
Mtheta = 32;
geom   = build_circle_mesh(p0.Sd, Nrad, Mtheta);

env = struct();
env.c   = p0.c;
env.rho = p0.rho;
env.r_ref       = 1.0;
env.off_axis    = 0;
env.azimuth     = 0;
env.green_scale = 2;
env.enc = struct('type','ib','Vb',[],'Rloss',[],'Rleak',[]);

% ============================================================
% D) Load + clean data
% ============================================================
df_clean = filter_impedance_csv(cfg.csvFile, cfg.cohThreshold);

% ============================================================
% E) SMART FIT
% ============================================================
result = fit_speaker_params_smart(df_clean, p0, geom, env, cfg);

clc;

disp('====================================================');
disp('SMART FIT RESULTS');
disp('====================================================');

disp('--- Fitted Parameters ---');
disp(result.p);

% ------------------------------------------------------------
% Derived quantities
% ------------------------------------------------------------
disp('--- Derived ---');
fprintf('Fs (Hz): %.6f\n', result.derived.fs);
fprintf('Qms    : %.6f\n', result.derived.Qms);
fprintf('Qes    : %.6f\n', result.derived.Qes);
fprintf('Qts    : %.6f\n', result.derived.Qts);
fprintf('Vas (m^3): %.6e\n', result.derived.Vas);

% ------------------------------------------------------------
% Confidence intervals
% ------------------------------------------------------------
disp('--- 95% Confidence Intervals (Parameters) ---');
for i = 1:numel(result.linearNames)
    name = result.linearNames{i};
    ci = result.ci95_linear(i,:);
    fprintf('%s: [%.4g, %.4g]\n', name, ci(1), ci(2));
end

disp('--- 95% Confidence Intervals (Derived) ---');
disp(result.derivedCI);

% ------------------------------------------------------------
% Conditioning
% ------------------------------------------------------------
fprintf('cond(J_data) = %.3e\n', result.condJ);

% ============================================================
% F) Plot
% ============================================================
plot_fit_result(result, cfg);
grid on;

% ============================================================
% Local helper functions
% ============================================================
function df_valid = filter_impedance_csv(filename, coh_threshold)
    if nargin < 2 || isempty(coh_threshold), coh_threshold = 0.8; end

    df = readtable(filename);
    validate_columns(df, {'frequency_hz','Z_mag_ohm','Z_phase_deg','Coherence'});

    mask = df.Coherence >= coh_threshold;
    df_valid = sortrows(df(mask, :), 'frequency_hz');

    fprintf('After coherence filter: %d\n', height(df_valid));
end

function validate_columns(df, req)
    missing = req(~ismember(req, df.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing required column(s): %s', strjoin(missing, ', '));
    end
end

function plot_fit_result(result, cfg)

    f  = result.f_plot;
    Zm = result.Zmeas_plot;
    Zf = result.Z_fit_plot;

    fs = result.fs_refined;

    figure('Name','Inverse Fit Result','Position',[350 250 1100 700]);

    % =========================
    % Magnitude
    % =========================
    subplot(2,1,1)
    semilogx(f, abs(Zm), 'LineWidth', 1.1); hold on;
    semilogx(f, abs(Zf), '--', 'LineWidth', 1.2);

    xline(fs, ':', 'Fs');

    if isfield(cfg, 'bandLow') && isfield(cfg, 'bandHigh')
        xline(cfg.bandLow * fs, '--', '0.7Fs');
        xline(cfg.bandHigh * fs, '--', '3Fs');
    end

    ylabel('|Z| (Ohm)');
    legend('Measured','Fitted');
    grid on;

    % =========================
    % Phase
    % =========================
    subplot(2,1,2)
    semilogx(f, unwrap(angle(Zm))*180/pi, 'LineWidth', 1.1); hold on;
    semilogx(f, unwrap(angle(Zf))*180/pi, '--', 'LineWidth', 1.2);

    xline(fs, ':', 'Fs');

    if isfield(cfg, 'bandLow') && isfield(cfg, 'bandHigh')
        xline(cfg.bandLow * fs, '--', '0.7Fs');
        xline(cfg.bandHigh * fs, '--', '3Fs');
    end

    ylabel('Phase (deg)');
    xlabel('Frequency (Hz)');
    legend('Measured','Fitted');
    grid on;
end