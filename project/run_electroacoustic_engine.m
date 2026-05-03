function results = run_electroacoustic_engine(f, p, geom, env)
%RUN_ELECTROACOUSTIC_ENGINE  Orchestrates the full electroacoustic solve.
%
% Required p fields:
%   Re, Le, Bl, Rms, Mms, Cms
%
% Required env fields:
%   geom        : mesh structure with fields xc, A, aeq, Np
%   c           : sound speed [m/s]
%   rho         : density [kg/m^3]
%   r_ref       : reference distance [m]
%   off_axis    : polar angle [rad]
%   azimuth     : azimuth angle [rad]
%
% Optional env fields:
%   green_scale : 2 for half-space / IB, 1 for free-space
%   enc         : enclosure struct for z_box()

    f = f(:).';                 % row vector
    w = 2*pi*f;
    s = 1i * w;

    validate_transducer(p);

    assert(isfield(env, 'c')    && ~isempty(env.c),    'env.c is required');
    assert(isfield(env, 'rho')  && ~isempty(env.rho),  'env.rho is required');


    if ~isfield(env, 'r_ref') || isempty(env.r_ref)
        env.r_ref = 1.0;
    end
    if ~isfield(env, 'off_axis') || isempty(env.off_axis)
        env.off_axis = 0;
    end
    if ~isfield(env, 'azimuth') || isempty(env.azimuth)
        env.azimuth = 0;
    end
    if ~isfield(env, 'green_scale') || isempty(env.green_scale)
        env.green_scale = 2;   % IB / half-space default
    end

    % Build reusable acoustic operator once
    Aac = build_acoustic_operator(f, env, geom);

    % Mechanical subsystem
    Zdriver = z_driver(s, p);
    Zbox    = z_box(s, p, env);

    % Rigid-piston radiation impedance assumption for now
    v_unit = ones(geom.Np, 1);

    Zrad_bem = zeros(1, numel(f));
    G_BEM    = zeros(1, numel(f));

    xobs = [env.r_ref * sin(env.off_axis) * cos(env.azimuth), ...
            env.r_ref * sin(env.off_axis) * sin(env.azimuth), ...
            env.r_ref * cos(env.off_axis)];

    for fi = 1:numel(f)
        Zrad_bem(fi) = Aac.impedance(fi, v_unit);
        G_BEM(fi)    = Aac.pressure_at(fi, xobs, v_unit);
    end

    Zm_tot = Zdriver + Zrad_bem + Zbox;

    % Electrical impedance
    Ze = p.Re + s*p.Le + (p.Bl.^2) ./ Zm_tot;

    % Voltage-to-velocity
    Hv = p.Bl ./ (Zm_tot .* Ze);

    % Displacement
    Hx = Hv ./ s;
    Hx(f == 0) = 0;

    % Voltage-to-pressure
    Hp = G_BEM .* Hv;

    % Pack
    results = struct();
    results.f        = f;
    results.s        = s;
    results.Zdriver  = Zdriver;
    results.Zrad_bem = Zrad_bem;
    results.Zbox     = Zbox;
    results.Zm_tot   = Zm_tot;
    results.Ze       = Ze;
    results.Hv       = Hv;
    results.Hx       = Hx;
    results.G_BEM    = G_BEM;
    results.Hp       = Hp;
    results.geom     = geom;
    results.env      = env;
    results.p        = p;
    results.Aac      = Aac;
end