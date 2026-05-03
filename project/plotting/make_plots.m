function make_plots(results, metrics, V_in, env)
%MAKE_PLOTS  3-panel plot: displacement, SPL, efficiency.

    if nargin < 4 || isempty(env)
        env = struct();
    end

    if ~isfield(env, 'p_ref') || isempty(env.p_ref)
        env.p_ref = 20e-6;
    end

    f = results.f(:).';

    x_mm = abs(results.Hx .* V_in) * 1e3;
    p_out = abs(results.Hp .* V_in);
    SPL = 20 * log10(max(p_out, eps) / env.p_ref);
    eta_pct = 100 * abs(metrics.eta);

    figure('Name','Speaker Physics (True Power Model)','Position',[350 350 1100 600]);

    subplot(3,1,1);
    semilogx(f, x_mm, 'LineWidth', 1.2);
    grid on; xlim([min(f) max(f)]);
    ylabel('x (mm)');
    title('Cone Displacement');

    subplot(3,1,2);
    semilogx(f, SPL, 'LineWidth', 1.2);
    grid on; xlim([min(f) max(f)]);
    ylabel('SPL (dB re 20 \muPa)');
    title('Radiated SPL');

    subplot(3,1,3);
    semilogx(f, eta_pct, 'LineWidth', 1.2);
    grid on; xlim([min(f) max(f)]);
    xlabel('Frequency (Hz)');
    ylabel('\eta (%)');
    title('True Electroacoustic Efficiency (Pa / Pe)');
end