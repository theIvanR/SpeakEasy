function results = simulate_speaker_full(f, p, opts)
%SIMULATE_SPEAKER_FULL Full loudspeaker comparison model.
%
% Always computes and compares:
%   - Infinite baffle (IB)
%   - Sealed enclosure
%
% Radiation model:
%   opts.radiationModel = 'constant' | 'bessel' | 'bem'
%
% Returns a structured result and optionally plots.

    if nargin < 3
        opts = struct();
    end

    % -------------------- defaults --------------------
    if ~isfield(opts, 'cfg') || isempty(opts.cfg), opts.cfg = 'single'; end
    if ~isfield(opts, 'Vb') || isempty(opts.Vb), opts.Vb = 0.2; end
    if ~isfield(opts, 'alpha') || isempty(opts.alpha), opts.alpha = 1.25; end
    if ~isfield(opts, 'radiationModel') || isempty(opts.radiationModel), opts.radiationModel = 'bem'; end
    if ~isfield(opts, 'r_ref') || isempty(opts.r_ref), opts.r_ref = 1.0; end
    if ~isfield(opts, 'plot') || isempty(opts.plot), opts.plot = 'full';  end

    % -------------------- validate --------------------
    validate_inputs(f, p, opts);

    % -------------------- unpack --------------------
    f = f(:);
    w = 2*pi*f;

    Re  = p.Re;
    Le  = p.Le;
    Bl  = p.Bl;
    Rms = p.Rms;
    Mms = p.Mms;
    Cms = p.Cms;
    Sd  = p.Sd;
    rho = p.rho;
    c   = p.c;

    % -------------------- configuration --------------------
    [Cms_free, Sd_eff] = configure_geometry(Cms, Sd, opts.cfg);
    a = sqrt(Sd_eff/pi);

    % -------------------- common radiation impedance --------------------
    [R_rad_raw, X_rad_raw] = compute_radiation_impedance( ...
        f, w, a, rho, c, opts.radiationModel);

    R_rad_w = R_rad_raw(:);
    X_rad_w = X_rad_raw(:);

    % keep raw model output untouched; use a sanitized copy where needed
    R_rad_for_power = max(R_rad_w, 0);

    % pressure proxy at reference distance
    H_rad = (1i .* w .* rho .* Sd_eff) ./ (2*pi*opts.r_ref);

    % -------------------- solve scenarios --------------------
    common = struct();
    common.f = f;
    common.w = w;
    common.Re = Re;
    common.Le = Le;
    common.Bl = Bl;
    common.Rms = Rms;
    common.Mms = Mms;
    common.Cms_free = Cms_free;
    common.Sd_eff = Sd_eff;
    common.rho = rho;
    common.c = c;
    common.H_rad = H_rad;
    common.R_rad = R_rad_w;
    common.X_rad = X_rad_w;
    common.R_rad_for_power = R_rad_for_power;

    freeAir = solve_scenario(common, struct( ...
        'name', 'IB', ...
        'useBox', false, ...
        'Vb', NaN, ...
        'alpha', NaN));

    sealedBox = solve_scenario(common, struct( ...
        'name', 'Sealed', ...
        'useBox', true, ...
        'Vb', opts.Vb, ...
        'alpha', opts.alpha));

    % -------------------- summary prints --------------------
    fprintf('\n=== COMPARISON SUMMARY ===\n');
    fprintf('cfg: %s | radiationModel: %s | Sd_eff = %.6g m^2\n', opts.cfg, opts.radiationModel, Sd_eff);

    fprintf('\n-- Infinite Baffle (IB) --\n');
    fprintf('analytic no-Ma: f_res = %.4f Hz\n', freeAir.res_noMa.f_res_analytic);
    fprintf('numeric no-Ma:  f_res = %.4f Hz\n', freeAir.res_noMa.f_res_num);
    fprintf('numeric with Ma: f_res = %.4f Hz\n', freeAir.res_full.f_res_num);
    fprintf('Ma_at_res = %.6g kg | M_eff = %.6g kg | Re[Zm]@res = %.6g | Q = %.4f\n', ...
        freeAir.res_full.Ma_at_res, freeAir.res_full.M_eff_res, freeAir.res_full.ReZm_at_res, freeAir.res_full.Q_from_Z);

    fprintf('\n-- Sealed box --\n');
    fprintf('Vb = %.4f m^3 | alpha = %.3f | Cms_box = %.6e m/N\n', sealedBox.Vb, sealedBox.alpha, sealedBox.Cms_box);
    fprintf('analytic sealed (no-Ma): f_c = %.4f Hz\n', sealedBox.f_c_analytic);
    fprintf('self-consistent sealed with Ma: f_c = %.4f Hz\n', sealedBox.f_c);
    fprintf('numeric no-Ma:  f_res = %.4f Hz\n', sealedBox.res_noMa.f_res_num);
    fprintf('numeric with Ma: f_res = %.4f Hz\n', sealedBox.res_full.f_res_num);
    fprintf('Ma_at_res = %.6g kg | M_eff = %.6g kg | Re[Zm]@res = %.6g | Q = %.4f\n', ...
        sealedBox.res_full.Ma_at_res, sealedBox.res_full.M_eff_res, sealedBox.res_full.ReZm_at_res, sealedBox.res_full.Q_from_Z);

    % -------------------- plots --------------------
    plot_results(freeAir, sealedBox, f, opts);


    % -------------------- outputs --------------------
    results = struct();
    results.f = f;
    results.w = w;
    results.cfg = opts.cfg;
    results.radiationModel = lower(opts.radiationModel);
    results.Sd_eff = Sd_eff;
    results.rho = rho;
    results.c = c;
    results.R_rad = R_rad_w;
    results.X_rad = X_rad_w;

    results.IB = freeAir;
    results.Sealed = sealedBox;

    results.summary = struct( ...
        'Cms_free', Cms_free, ...
        'f_res_free', 1/(2*pi*sqrt(Mms * Cms_free)), ...
        'V_as', rho * c^2 * Cms_free * Sd_eff^2);
end

% =====================================================================
% Helpers
% =====================================================================

function validate_inputs(f, p, opts)
    if ~isvector(f) || isempty(f)
        error('f must be a nonempty vector of positive frequencies.');
    end
    if any(~isfinite(f)) || any(f <= 0)
        error('f must contain only positive finite values.');
    end
    if any(diff(f(:)) <= 0)
        error('f must be strictly increasing.');
    end

    required = {'Re','Le','Bl','Rms','Mms','Cms','Sd','rho','c'};
    for k = 1:numel(required)
        if ~isfield(p, required{k})
            error('Missing field p.%s.', required{k});
        end
        if ~isscalar(p.(required{k})) || ~isfinite(p.(required{k}))
            error('p.%s must be a finite scalar.', required{k});
        end
    end

    if ~isfield(opts, 'cfg') || isempty(opts.cfg)
        return;
    end

    if ~ismember(lower(opts.cfg), {'single','series','parallel','series_parallel'})
        error('Unknown cfg. Use ''single'', ''series'', ''parallel'', or ''series_parallel''.');
    end

    if ~ismember(lower(opts.radiationModel), {'constant','bessel','bem'})
        error('Unknown radiationModel. Use ''constant'', ''bessel'', or ''bem''.');
    end

    if ~isscalar(opts.r_ref) || ~isfinite(opts.r_ref) || opts.r_ref <= 0
        error('opts.r_ref must be a positive finite scalar.');
    end

    if ~isscalar(opts.Vb) || ~isfinite(opts.Vb) || opts.Vb <= 0
        error('opts.Vb must be a positive finite scalar.');
    end

    if ~isscalar(opts.alpha) || ~isfinite(opts.alpha) || opts.alpha < 0
        error('opts.alpha must be a nonnegative finite scalar.');
    end
end

function [Cms_eff, Sd_eff] = configure_geometry(Cms, Sd, cfg)
    switch lower(cfg)
        case 'single'
            Cms_eff = Cms;
            Sd_eff  = Sd;
        case 'series'
            Cms_eff = Cms/2;
            Sd_eff  = Sd;
        case 'parallel'
            Cms_eff = 2*Cms;
            Sd_eff  = 2*Sd;
        case 'series_parallel'
            Cms_eff = Cms;
            Sd_eff  = 2*Sd;
        otherwise
            error('Unknown cfg.');
    end
end

function [R_rad, X_rad] = compute_radiation_impedance(f, w, a, rho, c, radiationModel)
    switch lower(radiationModel)
        case 'constant'
            ka = w*a/c;
            R_rad = rho * c * (pi*a^2) .* (ka.^2 / 2);
            X_rad = rho * c * (pi*a^2) .* (8 * ka / (3*pi));

        case 'bessel'
            ka = w*a/c;
            ka_safe = max(ka, 1e-12);
            J1_term = besselj(1, 2*ka_safe) ./ ka_safe;
            Y1_term = bessely(1, 2*ka_safe) ./ ka_safe;
            R_rad = rho * c * (pi*a^2) .* (1 - J1_term);
            X_rad = rho * c * (pi*a^2) .* (2/pi) .* Y1_term;

        case 'bem'
            [R_rad, X_rad] = pistonBEM(a, f, rho, c, 100);

        otherwise
            error('Unknown radiationModel.');
    end
end

function scenario = solve_scenario(common, cfg)
    f = common.f;
    w = common.w;

    Re  = common.Re;
    Le  = common.Le;
    Bl  = common.Bl;
    Rms = common.Rms;
    Mms = common.Mms;

    Cms_free = common.Cms_free;
    R_rad    = common.R_rad;
    X_rad    = common.X_rad;
    H_rad    = common.H_rad;

    if cfg.useBox
        Vb = cfg.Vb;
        alpha = cfg.alpha;

        C_box   = Vb / (common.rho * common.c^2 * common.Sd_eff^2);
        Cms_box = (Cms_free * C_box) / (Cms_free + C_box);
        R_box   = alpha * (common.rho * common.c * common.Sd_eff^2 / Vb);

        f_c_analytic = 1 / (2*pi*sqrt(Mms * Cms_box));
        f_c = estimate_self_consistent_box_resonance(f, X_rad, Mms, Cms_box, f_c_analytic);
    else
        Vb = NaN;
        alpha = NaN;
        C_box = Inf;
        Cms_box = Cms_free;
        R_box = 0;

        f_c_analytic = 1 / (2*pi*sqrt(Mms * Cms_free));
        f_c = NaN;
    end

    Ma_w = X_rad ./ (w + eps);

    Zm_full = (Rms + R_rad + R_box) + 1i .* w .* (Mms + Ma_w) + 1 ./ (1i .* w .* Cms_box);
    Zm_noMa = (Rms + R_rad + R_box) + 1i .* w .* Mms + 1 ./ (1i .* w .* Cms_box);

    Ze_full = Re + 1i .* w .* Le + (Bl^2) ./ Zm_full;
    Ze_noMa = Re + 1i .* w .* Le + (Bl^2) ./ Zm_noMa;

    H_v = Bl ./ (Zm_full .* Ze_full);
    H_p = H_rad .* H_v;

    % derived observables
    p_ref = 20e-6;
    spl = 20*log10(abs(H_p) / p_ref);

    H_x = H_v ./ (1i .* w);
    x_mm = abs(H_x) * 1e3;

    I = 1 ./ Ze_full;
    P_ac = 0.5 * abs(H_v).^2 .* common.R_rad_for_power;
    P_el = 0.5 * abs(I).^2 .* real(Ze_full);
    eta = 100 * (P_ac ./ max(P_el, eps));

    Zmag = abs(Ze_full);
    Zph  = angle(Ze_full) * 180/pi;

    [~, idx_spl] = max(spl);
    [~, idx_x] = max(x_mm);
    [~, idx_eta] = max(eta);

    res_full = extract_resonance(f, w, Zm_full, Cms_box, Mms, X_rad);
    res_noMa = extract_resonance(f, w, Zm_noMa, Cms_box, Mms, X_rad);

    scenario = struct();
    scenario.name = cfg.name;
    scenario.Vb = Vb;
    scenario.alpha = alpha;
    scenario.C_box = C_box;
    scenario.Cms_box = Cms_box;
    scenario.R_box = R_box;
    scenario.f_c_analytic = f_c_analytic;
    scenario.f_c = f_c;

    scenario.Zm_full = Zm_full;
    scenario.Zm_noMa = Zm_noMa;
    scenario.Ze_full = Ze_full;
    scenario.Ze_noMa = Ze_noMa;
    scenario.H_v = H_v;
    scenario.H_p = H_p;
    scenario.spl = spl;
    scenario.x_mm = x_mm;
    scenario.eta = eta;
    scenario.Zmag = Zmag;
    scenario.Zph = Zph;
    scenario.P_ac = P_ac;
    scenario.P_el = P_el;
    scenario.res_full = res_full;
    scenario.res_noMa = res_noMa;

    scenario.peaks = struct( ...
        'idx_spl', idx_spl, ...
        'idx_x', idx_x, ...
        'idx_eta', idx_eta);
end

function f_c = estimate_self_consistent_box_resonance(f, X_rad, Mms, Cms_box, f_init)
    f_old = f_init;
    tol = 1e-9;
    maxit = 100;

    for it = 1:maxit
        omega = 2*pi*max(f_old, 1e-9);
        X_rad_fc = interp1(f, X_rad, f_old, 'linear', 'extrap');
        Ma_fc = X_rad_fc / (omega + eps);
        M_eff = Mms + Ma_fc;
        f_new = 1 / (2*pi*sqrt(M_eff * Cms_box));

        if abs(f_new - f_old) < tol
            f_old = f_new;
            break;
        end
        f_old = f_new;
    end

    f_c = f_old;
end

function res = extract_resonance(f, w, Zm, Cms_local, Mms, X_rad)
    % Prefer a real zero-crossing of imag(Zm) if available; otherwise fall
    % back to the minimum of abs(imag(Zm)). This makes the estimate less
    % sensitive to the frequency grid.

    y = imag(Zm);

    f_res_num = find_resonance_frequency(f, y);
    omega_res = 2*pi*f_res_num;

    f_res_analytic = 1/(2*pi*sqrt(Mms * Cms_local));

    Zm_res = interp1(f, Zm, f_res_num, 'linear', 'extrap');
    X_res  = interp1(f, X_rad, f_res_num, 'linear', 'extrap');

    Ma_at_res = X_res / (omega_res + eps);
    M_eff_res = Mms + Ma_at_res;
    ReZm_at_res = real(Zm_res);

    Q_from_Z = omega_res * M_eff_res / max(ReZm_at_res, eps);

    res = struct( ...
        'f_res_num', f_res_num, ...
        'f_res_analytic', f_res_analytic, ...
        'omega_res', omega_res, ...
        'Ma_at_res', Ma_at_res, ...
        'M_eff_res', M_eff_res, ...
        'ReZm_at_res', ReZm_at_res, ...
        'Q_from_Z', Q_from_Z);
end

function f_res = find_resonance_frequency(f, y)
    % Locate a sign change near the minimum |imag(Zm)|, if possible.
    % Otherwise return the grid point where |y| is smallest.

    [~, idx_min] = min(abs(y));

    % exact zero on grid
    if abs(y(idx_min)) < 1e-14
        f_res = f(idx_min);
        return;
    end

    % detect sign changes
    s = sign(y);
    s(s == 0) = 1;
    zc = find(s(1:end-1) .* s(2:end) < 0);

    if isempty(zc)
        f_res = f(idx_min);
        return;
    end

    % choose the sign-change interval closest to the minimum-|y| point
    mid_f = 0.5 * (f(zc) + f(zc+1));
    [~, j] = min(abs(mid_f - f(idx_min)));
    i = zc(j);

    % linear interpolation of root in that interval
    y1 = y(i);
    y2 = y(i+1);
    f1 = f(i);
    f2 = f(i+1);

    if abs(y2 - y1) < eps
        f_res = 0.5 * (f1 + f2);
    else
        f_res = f1 - y1 * (f2 - f1) / (y2 - y1);
    end
end

function plot_results(freeAir, sealedBox, f, opts)
    mode = lower(string(opts.plot));

    switch mode
        case "none"
            return;

        case "reduced"
            plot_combined_panel(freeAir, sealedBox, f);

        case "full"
            % SPL
            figure('Name','SPL: IB vs Sealed','Position',[100 100 900 400]);
            semilogx(f, freeAir.spl, 'b', 'LineWidth', 1.6); hold on;
            semilogx(f, sealedBox.spl, 'g', 'LineWidth', 1.6);
            grid on; xlim([10 2e4]);
            xlabel('Frequency (Hz)');
            ylabel('SPL proxy (dB re 20 \muPa)');
            title('SPL: Infinite Baffle vs Sealed');
            xline(freeAir.res_full.f_res_num, 'b--', 'Label', 'IB f_{res}', 'LabelVerticalAlignment', 'bottom');
            xline(sealedBox.res_full.f_res_num, 'g--', 'Label', 'Sealed f_{res}', 'LabelVerticalAlignment', 'bottom');
            legend('IB SPL','Sealed SPL','Location','best');

            % displacement
            figure('Name','Cone displacement','Position',[150 150 900 400]);
            semilogx(f, freeAir.x_mm, 'b', 'LineWidth', 1.4); hold on;
            semilogx(f, sealedBox.x_mm, 'g', 'LineWidth', 1.4);
            grid on; xlim([10 2e4]);
            xlabel('Frequency (Hz)');
            ylabel('Displacement (mm/V)');
            title('Cone displacement: IB vs Sealed');
            xline(freeAir.res_full.f_res_num, 'b--', 'Label', 'IB f_{res}', 'LabelVerticalAlignment', 'bottom');
            xline(sealedBox.res_full.f_res_num, 'g--', 'Label', 'Sealed f_{res}', 'LabelVerticalAlignment', 'bottom');
            legend('IB displacement','Sealed displacement','Location','best');

            % efficiency
            figure('Name','Efficiency comparison','Position',[200 200 900 350]);
            semilogx(f, freeAir.eta, 'b', 'LineWidth', 1.4); hold on;
            semilogx(f, sealedBox.eta, 'g', 'LineWidth', 1.4);
            grid on; xlim([10 2e4]);
            xlabel('Frequency (Hz)');
            ylabel('Efficiency (%)');
            title('Reference efficiency: IB vs Sealed');
            xline(f(idx_of_max(freeAir.eta)), 'b--', 'Label', sprintf('IB peak %.2f%%', max(freeAir.eta)));
            xline(f(idx_of_max(sealedBox.eta)), 'g--', 'Label', sprintf('Sealed peak %.2f%%', max(sealedBox.eta)));
            legend('IB \eta','Sealed \eta','Location','best');

            % impedance
            figure('Name','Electrical impedance','Position',[250 250 900 500]);
            subplot(2,1,1);
            semilogx(f, freeAir.Zmag, 'b', 'LineWidth', 1.4); hold on;
            semilogx(f, sealedBox.Zmag, 'g', 'LineWidth', 1.4);
            grid on; xlim([10 2e4]);
            ylabel('|Z_e| (\Omega)');
            title('Electrical impedance magnitude');
            legend('IB |Z_e|','Sealed |Z_e|','Location','best');

            subplot(2,1,2);
            semilogx(f, freeAir.Zph, 'b', 'LineWidth', 1.1); hold on;
            semilogx(f, sealedBox.Zph, 'g', 'LineWidth', 1.1);
            grid on; xlim([10 2e4]);
            xlabel('Frequency (Hz)');
            ylabel('Phase (deg)');
            title('Electrical impedance phase');
            legend('IB \angle Z_e','Sealed \angle Z_e','Location','best');

            % powers
            figure('Name','Powers','Position',[300 300 900 400]);
            semilogx(f, freeAir.P_ac, 'b', 'LineWidth', 1.2); hold on;
            semilogx(f, freeAir.P_el, 'b--', 'LineWidth', 1.0);
            semilogx(f, sealedBox.P_ac, 'g', 'LineWidth', 1.2);
            semilogx(f, sealedBox.P_el, 'g--', 'LineWidth', 1.0);
            grid on; xlim([10 2e4]);
            set(gca, 'YScale', 'log');
            xlabel('Frequency (Hz)');
            ylabel('Power (W)');
            title('Radiated acoustic power (solid) and electrical input power (dashed)');
            legend('IB P_{ac}','IB P_{el}','Sealed P_{ac}','Sealed P_{el}','Location','best');

            % combined panel
            plot_combined_panel(freeAir, sealedBox, f);

        otherwise
            error('Unknown plot mode. Use ''none'', ''reduced'', or ''full''.');
    end
end

function plot_combined_panel(freeAir, sealedBox, f)
    figure('Name','Combined comparison','Position',[350 350 1100 600]);

    subplot(3,1,1);
    semilogx(f, freeAir.x_mm, 'b', 'LineWidth', 1.2); hold on;
    semilogx(f, sealedBox.x_mm, 'g', 'LineWidth', 1.2);
    grid on; xlim([10 2e4]);
    ylabel('x (mm/V)');
    title('Displacement comparison');
    legend('IB','Sealed');

    subplot(3,1,2);
    semilogx(f, freeAir.spl, 'b', 'LineWidth', 1.2); hold on;
    semilogx(f, sealedBox.spl, 'g', 'LineWidth', 1.2);
    grid on; xlim([10 2e4]);
    ylabel('SPL proxy (dB)');
    title('SPL comparison');

    subplot(3,1,3);
    semilogx(f, freeAir.eta, 'b', 'LineWidth', 1.2); hold on;
    semilogx(f, sealedBox.eta, 'g', 'LineWidth', 1.2);
    grid on; xlim([10 2e4]);
    xlabel('Frequency (Hz)');
    ylabel('\eta (%)');
    title('Efficiency comparison');
    legend('IB','Sealed');
end

function idx = idx_of_max(x)
    [~, idx] = max(x);
end