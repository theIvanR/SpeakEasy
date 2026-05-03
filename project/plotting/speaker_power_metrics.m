function metrics = speaker_power_metrics(results, V_in)
%SPEAKER_POWER_METRICS  True electrical + acoustic power efficiency.

    if nargin < 2 || isempty(V_in)
        V_in = 1;
    end

    f = results.f(:).';

    Ze = results.Ze;
    I  = V_in ./ Ze;

    % Time-averaged electrical input power
    Pe = 0.5 * real(V_in .* conj(I));

    % Cone velocity
    v = results.Hv .* V_in;

    % Acoustic power
    Pa = 0.5 * abs(v).^2 .* real(results.Zrad_bem);

    eta = Pa ./ max(Pe, eps);

    metrics = struct();
    metrics.f   = f;
    metrics.Pe  = Pe;
    metrics.Pa  = Pa;
    metrics.eta = eta;
    metrics.I   = I;
    metrics.v   = v;
end