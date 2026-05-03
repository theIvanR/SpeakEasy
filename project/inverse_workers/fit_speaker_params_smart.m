function result = fit_speaker_params_smart(df_clean, p0, geom, env, cfg)
%FIT_SPEAKER_PARAMS_SMART  Stable MAP inverse fitter for loudspeaker TS parameters.
%
% Fitted log-parameter vector:
%   x = [log(Re), log(Le), log(Bl), log(Rms), log(Mms), log(Cms)]
%
% Main fixes:
%   - soft fs anchor to measured fs0
%   - pass-1 broad physical bounds
%   - pass-2 Fisher-informed adaptive trust bounds
%   - finite-difference Jacobian
%   - Fisher/SVD prior for sloppy directions
%
% Inputs
%   df_clean : table with columns frequency_hz, Z_mag_ohm, Z_phase_deg, Coherence
%   p0       : initial parameter struct
%   geom     : geometry struct
%   env      : environment struct
%   cfg      : configuration struct
%
% Output
%   result   : fit, covariance, derived quantities, and diagnostics

    if nargin < 5 || isempty(cfg)
        cfg = struct();
    end

    cfg = apply_defaults(cfg, struct( ...
        'bandLow', 0.7, ...
        'bandHigh', 3.0, ...
        'stride', 1, ...
        'cohThreshold', [], ...
        'freqFocusStrength', 0.45, ...
        'freqFocusWidth', 0.7, ...
        'fsAnchorLambda', 100, ...
        'maxIter', 100, ...
        'maxEval', 1200, ...
        'fdStep', 1e-5, ...
        'verbose', true, ...
        'regLambda', 1.0, ...
        'priorCenterMode', 'pass1', ...
        'priorSigmaMin', 0.05, ...
        'priorSigmaMax', 10.0, ...
        'priorSharpness', 1.5, ...
        'priorUsePass1SVD', true, ...
        'includePriorInFisher', true, ...
        'trustMinLogWidth', log(1.10), ...
        'trustMaxLogWidth', log(2.0), ...
        'fisherTrustSigma', 2.0));

    validate_columns(df_clean, {'frequency_hz','Z_mag_ohm','Z_phase_deg','Coherence'});

    % ---------------------------------------------------------------------
    % Data unpacking
    % ---------------------------------------------------------------------
    f_all   = df_clean.frequency_hz(:).';
    Zmag    = df_clean.Z_mag_ohm(:);
    Zphase  = deg2rad(df_clean.Z_phase_deg(:));
    Zmeas   = Zmag .* exp(1j * Zphase);
    coh_all = df_clean.Coherence(:);

    % ---------------------------------------------------------------------
    % Initial feature extraction
    % ---------------------------------------------------------------------
    features = estimate_initial_features(f_all, Zmeas);
    fs0 = features.fs;

    if ~isfinite(fs0) || fs0 <= 0
        error('Could not estimate an initial resonance frequency.');
    end

    cfg.fs0 = fs0;

    % ---------------------------------------------------------------------
    % Fit band selection around fs0
    % ---------------------------------------------------------------------
    fLow  = cfg.bandLow  * fs0;
    fHigh = cfg.bandHigh * fs0;
    bandMask = f_all >= fLow & f_all <= fHigh;

    if ~any(bandMask)
        warning('Automatic fit band is empty. Falling back to full cleaned trace.');
        bandMask = true(size(f_all));
        fLow = min(f_all);
        fHigh = max(f_all);
    end

    idxBand = find(bandMask);
    if cfg.stride > 1
        idxBand = idxBand(1:cfg.stride:end);
    end

    f_fit      = f_all(idxBand);
    Z_fit_meas = Zmeas(idxBand);
    coh_fit    = coh_all(idxBand);

    if numel(f_fit) < 8
        warning('Very few points in the selected fit band; fit quality may be limited.');
    end

    if cfg.verbose
        fprintf('[SPEAKER FIT] Band: %.6g Hz to %.6g Hz | points = %d\n', fLow, fHigh, numel(f_fit));
        fprintf('[SPEAKER FIT] Initial fs0 = %.6g Hz\n', fs0);
    end

    % ---------------------------------------------------------------------
    % Initial parameter vector in log-space
    %   x = [log(Re), log(Le), log(Bl), log(Rms), log(Mms), log(Cms)]
    % ---------------------------------------------------------------------
    Re0  = max(features.Re0, 1e-3);
    Le0  = max(features.Le0, 1e-7);
    Bl0  = max(getfield_with_default(p0, 'Bl', 1.0), 1e-3); %#ok<GFLD>
    Rms0 = max(getfield_with_default(p0, 'Rms', 1.0), 1e-5); %#ok<GFLD>
    Mms0 = max(getfield_with_default(p0, 'Mms', 1e-3), 1e-6); %#ok<GFLD>
    Cms0 = max(getfield_with_default(p0, 'Cms', 1e-6), 1e-9); %#ok<GFLD>

    x0 = log([Re0, Le0, Bl0, Rms0, Mms0, Cms0]);

    % Pass 1 broad physical bounds
    lb1 = log([1.0,  1e-6, 0.5, 0.1, 1e-4, 1e-7]);
    ub1 = log([32,   5e-3, 30,  50,  0.5,  1e-3]);

    opts = optimoptions('lsqnonlin', ...
        'Display', ternary(cfg.verbose, 'iter', 'off'), ...
        'MaxFunctionEvaluations', cfg.maxEval, ...
        'MaxIterations', cfg.maxIter, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12, ...
        'OptimalityTolerance', 1e-12);

    % =====================================================================
    % PASS 1: data + fs anchor, no prior
    % =====================================================================
    resid1 = @(x) speaker_fit_residual(x, f_fit, Z_fit_meas, coh_fit, geom, env, cfg);
    [x1, resnorm1, residual1, exitflag1, output1] = lsqnonlin(resid1, x0, lb1, ub1, opts); %#ok<ASGLU>

    p1 = x_to_params(x1, p0);
    [J1_data, ~, meta1] = build_jacobian_fd(x1, f_fit, Z_fit_meas, coh_fit, geom, env, cfg, []);
    [U1, S1, V1] = svd(J1_data, 'econ');
    singular1 = diag(S1);
    condJ1 = singular1(1) / max(singular1(end), eps);

    % ---------------------------------------------------------------------
    % Build Bayesian prior aligned with pass-1 SVD basis
    % ---------------------------------------------------------------------
    if strcmpi(cfg.priorCenterMode, 'pass1')
        x_prior = x1;
    else
        x_prior = x0;
    end

    if cfg.priorUsePass1SVD
        Vprior = V1;
        s = singular1(:);
        smax = max(s);
        if ~isfinite(smax) || smax <= 0
            sNorm = ones(size(s));
        else
            sNorm = s ./ smax;
        end

        % Stiff directions -> weaker penalty
        % Sloppy directions -> stronger penalty
        sigmaPrior = cfg.priorSigmaMin + ...
            (cfg.priorSigmaMax - cfg.priorSigmaMin) .* (sNorm .^ cfg.priorSharpness);

        sigmaPrior = max(min(sigmaPrior, cfg.priorSigmaMax), cfg.priorSigmaMin);
    else
        Vprior = eye(numel(x0));
        sigmaPrior = cfg.priorSigmaMax * ones(numel(x0),1);
    end

    prior = struct();
    prior.center = x_prior(:);
    prior.V = Vprior;
    prior.sigma = sigmaPrior(:);
    prior.weight = sqrt(max(cfg.regLambda, eps));
    prior.description = 'Anisotropic Gaussian prior aligned with pass-1 SVD directions';

    % ---------------------------------------------------------------------
    % Fisher-informed adaptive trust bounds for PASS 2
    % ---------------------------------------------------------------------
    [lb2, ub2, trustInfo] = fisher_trust_bounds(x_prior, J1_data, cfg);

    % =====================================================================
    % PASS 2: MAP fit with fs anchor + SVD prior
    % =====================================================================
    resid2 = @(x) speaker_residual(x, f_fit, Z_fit_meas, coh_fit, geom, env, cfg, prior);
    [xhat, resnorm, residual, exitflag, output] = lsqnonlin(resid2, x_prior, lb2, ub2, opts);

    phat = x_to_params(xhat, p0);
    sim  = run_electroacoustic_engine(f_fit, phat, geom, env);

    % ---------------------------------------------------------------------
    % Jacobians at the optimum
    % ---------------------------------------------------------------------
    [J_data, J_total, meta] = build_jacobian_fd(xhat, f_fit, Z_fit_meas, coh_fit, geom, env, cfg, prior);

    F_data = J_data.' * J_data;
    F_map  = J_total.' * J_total;

    if cfg.includePriorInFisher
        F_use = F_map;
    else
        F_use = F_data;
    end

    [U, S, V] = svd(J_data, 'econ');
    singular = diag(S);
    condJ = singular(1) / max(singular(end), eps);
    condF = cond(F_use);

    % Posterior covariance in log-parameter space
    if condF < 1e12
        Cov_log = inv(F_use);
    else
        Cov_log = pinv(F_use);
    end

    se_log = sqrt(max(diag(Cov_log), 0));
    corr = fisher_correlation(F_use);

    % 95% CIs in log-space and linear-space
    paramNames = {'log(Re)','log(Le)','log(Bl)','log(Rms)','log(Mms)','log(Cms)'};
    logVals = xhat(:);
    ci95_log = [logVals - 1.96 * se_log, logVals + 1.96 * se_log];

    linearNames = {'Re','Le','Bl','Rms','Mms','Cms'};
    pVec = exp(logVals);
    ci95_lin = [exp(ci95_log(:,1)), exp(ci95_log(:,2))];

    % Derived parameters and propagated CI
    derived = compute_derived_quantities(phat, env);
    [derivedCI, derivedCovLog] = propagate_derived_ci(xhat, Cov_log, phat, env);

    % ---------------------------------------------------------------------
    % Plotting data
    % ---------------------------------------------------------------------
    idxPlot = find(bandMask);
    if cfg.stride > 1
        idxPlot = idxPlot(1:cfg.stride:end);
    end
    f_plot     = f_all(idxPlot);
    Zmeas_plot = Zmeas(idxPlot);
    Zfit_plot  = run_electroacoustic_engine(f_plot, phat, geom, env).Ze;

    % ---------------------------------------------------------------------
    % Pack result
    % ---------------------------------------------------------------------
    result = struct();

    result.p = phat;
    result.pVec = pVec;
    result.f = f_all;
    result.Zmeas = Zmeas;
    result.coherence = coh_all;
    result.fitMask = bandMask;
    result.f_fit = f_fit;
    result.f_plot = f_plot;
    result.Zmeas_plot = Zmeas_plot;
    result.Z_fit_plot = Zfit_plot;
    result.sim = sim;

    result.fs0 = fs0;
    result.fs_refined = phat.fs;

    result.pass1 = struct( ...
        'x', x1, 'p', p1, 'resnorm', resnorm1, 'residual', residual1, ...
        'exitflag', exitflag1, 'output', output1, 'J_data', J1_data, ...
        'singular', singular1, 'sloppyVec', V1(:,end), ...
        'meta', meta1);

    result.fit = struct( ...
        'x', xhat, 'resnorm', resnorm, 'residual', residual, ...
        'exitflag', exitflag, 'output', output, 'meta', meta);

    result.prior = prior;

    result.J_data = J_data;
    result.J_total = J_total;
    result.F_data = F_data;
    result.F_map = F_map;
    result.F_used = F_use;
    result.U = U;
    result.S = S;
    result.V = V;
    result.singular = singular;
    result.condJ = condJ;
    result.condF = condF;
    result.corr = corr;
    result.paramNames = paramNames;
    result.linearNames = linearNames;

    result.Cov_log = Cov_log;
    result.se_log = se_log;
    result.ci95_log = ci95_log;
    result.ci95_linear = ci95_lin;

    result.derived = derived;
    result.derivedCI = derivedCI;
    result.derivedCovLog = derivedCovLog;

    result.trust = trustInfo;
    result.trust.lb2 = lb2;
    result.trust.ub2 = ub2;

    result.notes = struct();
    result.notes.parameterization = 'Re, Le, Bl, Rms, Mms, Cms with fs/Qts/Qes/Qms derived';
    result.notes.prior = prior.description;
    result.notes.priorCenterMode = cfg.priorCenterMode;
    result.notes.priorSigmaMin = cfg.priorSigmaMin;
    result.notes.priorSigmaMax = cfg.priorSigmaMax;
    result.notes.priorSharpness = cfg.priorSharpness;
    result.notes.includePriorInFisher = cfg.includePriorInFisher;
    result.notes.fsAnchorLambda = cfg.fsAnchorLambda;
    result.notes.fs0 = fs0;
    result.notes.trustMinLogWidth = cfg.trustMinLogWidth;
    result.notes.trustMaxLogWidth = cfg.trustMaxLogWidth;
    result.notes.fisherTrustSigma = cfg.fisherTrustSigma;

    if cfg.verbose
        fprintf('[SPEAKER FIT] Pass 1 fs0 = %.6g Hz | refined fs = %.6g Hz\n', fs0, phat.fs);
        fprintf('[SPEAKER FIT] cond(J_data) = %.3e | cond(F_used) = %.3e\n', condJ, condF);
        fprintf('[SPEAKER FIT] Fit complete. Exitflag = %d\n', exitflag);
    end
end

% ========================================================================
% Residuals
% ========================================================================
function r = speaker_residual(x, f, Zmeas, coh, geom, env, cfg, prior)
    r = speaker_fit_residual(x, f, Zmeas, coh, geom, env, cfg);

    % Bayesian prior residual in SVD basis:
    %   z = V'*(x - mu)
    %   r_prior = weight * z ./ sigma
    if nargin >= 8 && ~isempty(prior)
        delta = x(:) - prior.center(:);
        z = prior.V' * delta;
        r_prior = prior.weight * (z ./ max(prior.sigma, eps));
        r = [r; r_prior(:)];
    end
end

function r = speaker_fit_residual(x, f, Zmeas, coh, geom, env, cfg)
    p = x_to_params(x, struct());
    sim = run_electroacoustic_engine(f, p, geom, env);

    if ~isfield(sim, 'Ze')
        error('run_electroacoustic_engine must return sim.Ze');
    end

    dZ = sim.Ze(:) - Zmeas(:);

    % Mild resonance emphasis using derived fs
    w_freq = resonance_window(f(:), p.fs, cfg.freqFocusStrength, cfg.freqFocusWidth);

    % Scale by magnitude, avoid division by tiny numbers
    scale = max(abs(Zmeas(:)), 0.25 * max(abs(p.Re), 1e-6));
    scale = max(scale, 1e-6);

    w = max(coh(:), 0.05) .* w_freq;

    r_data = [real(dZ)./scale .* w; imag(dZ)./scale .* w];
    r = r_data;

    % Soft resonance anchor: keep fs near the measured zero-crossing estimate
    if isfield(cfg, 'fsAnchorLambda') && cfg.fsAnchorLambda > 0 && ...
            isfield(cfg, 'fs0') && isfinite(cfg.fs0) && cfg.fs0 > 0

        fs_model = p.fs;
        r_fs = sqrt(cfg.fsAnchorLambda) * (log(fs_model) - log(cfg.fs0));
        r = [r; r_fs];
    end
end

% ========================================================================
% Jacobian of the residual (finite differences for the fit residual;
% anchor row added analytically for stability)
% ========================================================================
function [J_data, J_total, meta] = build_jacobian_fd(x, f, Zmeas, coh, geom, env, cfg, prior)
    x = x(:);
    n = numel(x);

    r0 = speaker_fit_residual(x, f, Zmeas, coh, geom, env, cfg);
    m = numel(r0);
    J_data = zeros(m, n);

    for k = 1:n
        h = cfg.fdStep * max(1, abs(x(k)));
        dx = zeros(n, 1);
        dx(k) = h;

        rp = speaker_fit_residual(x + dx, f, Zmeas, coh, geom, env, cfg);
        rm = speaker_fit_residual(x - dx, f, Zmeas, coh, geom, env, cfg);

        J_data(:, k) = (rp - rm) / (2*h);
    end

    % Add prior Jacobian analytically
    if nargin >= 8 && ~isempty(prior)
        J_prior = prior.weight * (diag(1 ./ max(prior.sigma, eps)) * prior.V');
        J_total = [J_data; J_prior];
    else
        J_total = J_data;
    end

    p = x_to_params(x, struct());
    meta = struct('fs', p.fs, 'r0', r0);
end

% ========================================================================
% Initial feature estimation
% ========================================================================
function features = estimate_initial_features(f, Zmeas)
    f = f(:);
    Zmeas = Zmeas(:);

    ph = unwrap(angle(Zmeas));
    ph_deg = ph * 180/pi;
    fs = estimate_zero_crossing_frequency(f, ph_deg);

    Zmag = abs(Zmeas);
    idx_hf_start = max(round(0.8 * numel(f)), 1);
    idx_hf = idx_hf_start:numel(f);

    omega = 2*pi*f(idx_hf);
    Zhf = Zmeas(idx_hf);

    Re0 = median(real(Zhf), 'omitnan');
    Le0 = median(imag(Zhf) ./ max(omega, eps), 'omitnan');

    if ~isfinite(Re0) || Re0 <= 0
        Re0 = max(min(Zmag), 1e-3);
    end
    if ~isfinite(Le0) || Le0 <= 0
        Le0 = 1e-4;
    end

    features = struct('fs', fs, 'Re0', Re0, 'Le0', Le0);
end

function fs = estimate_zero_crossing_frequency(f, ph_deg)
    f = f(:);
    ph_deg = ph_deg(:);

    s = sign(ph_deg);
    s(s == 0) = NaN;

    idx = find(~isnan(s(1:end-1)) & ~isnan(s(2:end)) & (s(1:end-1) ~= s(2:end)), 1, 'first');

    if isempty(idx)
        [~, idx0] = min(abs(ph_deg));
        fs = f(idx0);
        return;
    end

    fs = interp1(ph_deg(idx:idx+1), f(idx:idx+1), 0, 'linear', 'extrap');
end

% ========================================================================
% Parameter mapping
% ========================================================================
function p = x_to_params(x, p0)
    x = x(:);
    p = p0;

    p.Re  = exp(x(1));
    p.Le  = exp(x(2));
    p.Bl  = exp(x(3));
    p.Rms = exp(x(4));
    p.Mms = exp(x(5));
    p.Cms = exp(x(6));

    p.fs = 1 / (2*pi*sqrt(max(p.Mms * p.Cms, eps)));
    ws = 2*pi*p.fs;

    p.Qms = ws * p.Mms / max(p.Rms, eps);
    p.Qes = ws * p.Mms * p.Re / max(p.Bl^2, eps);
    p.Qts = (p.Qms * p.Qes) / max(p.Qms + p.Qes, eps);
end

function q = compute_derived_quantities(p, env)
    ws = 2*pi*p.fs;
    Cms = p.Cms;
    Qms = ws * p.Mms / p.Rms;
    Qes = ws * p.Mms * p.Re / (p.Bl^2);
    Qts = (Qms * Qes) / max(Qms + Qes, eps);
    Vas = env.rho * env.c^2 * (p.Sd^2) * Cms;

    q = struct();
    q.fs = p.fs;
    q.Cms = Cms;
    q.Qms = Qms;
    q.Qes = Qes;
    q.Qts = Qts;
    q.Vas = Vas;
    q.Rms = p.Rms;
end

function [derivedCI, derivedCovLog] = propagate_derived_ci(xhat, Cov_log, p, env)
    q0 = compute_derived_quantities(p, env);
    qNames = {'fs','Cms','Qms','Qes','Qts','Vas','Rms'};

    epsStep = 1e-6;
    J = zeros(numel(qNames), numel(xhat));

    for k = 1:numel(xhat)
        dx = zeros(size(xhat));
        dx(k) = epsStep;

        p_plus = x_to_params(xhat + dx, p);
        p_minus = x_to_params(xhat - dx, p);

        q_plus = compute_derived_quantities(p_plus, env);
        q_minus = compute_derived_quantities(p_minus, env);

        for i = 1:numel(qNames)
            name = qNames{i};
            vp = q_plus.(name);
            vm = q_minus.(name);

            if vp > 0 && vm > 0
                J(i,k) = (log(vp) - log(vm)) / (2*epsStep);
            else
                J(i,k) = (vp - vm) / (2*epsStep);
            end
        end
    end

    derivedCovLog = J * Cov_log * J.';
    se = sqrt(max(diag(derivedCovLog), 0));

    derivedCI = struct();
    for i = 1:numel(qNames)
        name = qNames{i};
        mu = q0.(name);
        if mu > 0
            lo = exp(log(mu) - 1.96 * se(i));
            hi = exp(log(mu) + 1.96 * se(i));
        else
            lo = mu - 1.96 * se(i);
            hi = mu + 1.96 * se(i);
        end
        derivedCI.(name) = [lo, hi];
    end
end

% ========================================================================
% Fisher-informed adaptive trust bounds
% ========================================================================
function [lb, ub, info] = fisher_trust_bounds(x_center, J, cfg)
%FISHER_TRUST_BOUNDS  Axis-aligned trust bounds learned from Fisher information.
%
% This is a first-step approximation to a Fisher ellipsoid.
% Since lsqnonlin only accepts box bounds, we convert the local covariance
% into per-parameter half-widths and cap them to keep the search local.

    x_center = x_center(:);

    F = J.' * J;
    F = 0.5 * (F + F.');

    % Robust covariance estimate in log-parameter space
    if any(~isfinite(F(:))) || rcond(F) < 1e-14
        Cov = pinv(F + 1e-12 * eye(size(F)));
    else
        Cov = pinv(F);
    end

    sigma = sqrt(max(diag(Cov), eps));

    % Fisher-scaled half-widths
    rawHalfWidth = cfg.fisherTrustSigma * sigma;

    % Clamp to a sensible local neighborhood:
    % - not tighter than trustMinLogWidth
    % - not wider than trustMaxLogWidth
    halfWidth = min(max(rawHalfWidth, cfg.trustMinLogWidth), cfg.trustMaxLogWidth);

    lb = x_center - halfWidth(:);
    ub = x_center + halfWidth(:);

    info = struct();
    info.F = F;
    info.Cov = Cov;
    info.sigma = sigma;
    info.rawHalfWidth = rawHalfWidth;
    info.halfWidth = halfWidth;
end

% ========================================================================
% Resonance weighting helper
% ========================================================================
function w = resonance_window(f, fs, strength, widthFactor)
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        w = ones(size(f));
        return;
    end

    strength = min(max(strength, 0), 1);
    widthFactor = max(widthFactor, 0.05);

    g = exp(-((f - fs) ./ (widthFactor * fs)).^2);
    w = (1 - strength) + strength * g;
    w = max(w, 0.1);
end

% ========================================================================
% Fisher helpers
% ========================================================================
function corr = fisher_correlation(F)
    d = sqrt(max(diag(F), eps));
    corr = F ./ (d * d.');
end

% ========================================================================
% Small utilities
% ========================================================================
function s = getfield_with_default(st, field, default)
    if isstruct(st) && isfield(st, field) && ~isempty(st.(field))
        s = st.(field);
    else
        s = default;
    end
end

function out = apply_defaults(in, defaults)
    out = in;
    f = fieldnames(defaults);
    for i = 1:numel(f)
        if ~isfield(out, f{i}) || isempty(out.(f{i}))
            out.(f{i}) = defaults.(f{i});
        end
    end
end

function validate_columns(df, req)
    missing = req(~ismember(req, df.Properties.VariableNames));
    if ~isempty(missing)
        error('Missing required column(s): %s', strjoin(missing, ', '));
    end
end

function y = ternary(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end
