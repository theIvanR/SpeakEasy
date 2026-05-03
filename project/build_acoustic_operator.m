function A = build_acoustic_operator(f, env, geom)
%BUILD_ACOUSTIC_OPERATOR  Reusable Green operator for radiation.
%
% env.green_scale:
%   2 -> infinite baffle / half-space (default)
%   1 -> free space
%
% The operator can compute:
%   - pressure on the surface
%   - total force
%   - radiation impedance
%   - pressure at an observation point

    if ~isfield(env, 'green_scale') || isempty(env.green_scale)
        env.green_scale = 2;
    end

    f = f(:).';
    w = 2*pi*f;
    k = w / env.c;

    A = struct();
    A.f = f;
    A.w = w;
    A.k = k;
    A.env = env;
    A.geom = geom;
    A.green_scale = env.green_scale;

    % Precompute source-source distance matrix once
    xc = geom.xc;
    dx = xc(:,1) - xc(:,1).';
    dy = xc(:,2) - xc(:,2).';
    A.Rss = sqrt(dx.^2 + dy.^2);

    % Methods
    A.pressure_surface = @(fi, v) pressure_surface(A, fi, v);
    A.force             = @(fi, v) surface_force(A, fi, v);
    A.impedance         = @(fi, v) surface_impedance(A, fi, v);
    A.pressure_at       = @(fi, xobs, v) pressure_at(A, fi, xobs, v);
end

function psurf = pressure_surface(A, fi, v)
% Pressure on the radiating surface for a given velocity field v.

    v = v(:);
    if isscalar(v)
        v = repmat(v, A.geom.Np, 1);
    end

    G = source_kernel(A, fi);
    psurf = 1i * A.w(fi) * A.env.rho * (G * v);
end

function F = surface_force(A, fi, v)
% Total radiated force from surface pressure.

    psurf = pressure_surface(A, fi, v);
    F = sum(psurf .* A.geom.A);
end

function Z = surface_impedance(A, fi, v)
% Radiation impedance = force / area-weighted mean velocity.

    v = v(:);
    if isscalar(v)
        v = repmat(v, A.geom.Np, 1);
    end

    F = surface_force(A, fi, v);
    vavg = sum(v .* A.geom.A) / sum(A.geom.A);

    if abs(vavg) < eps
        Z = 0;
    else
        Z = F / vavg;
    end
end

function pobs = pressure_at(A, fi, xobs, v)
% Pressure at one observation point xobs = [x y z].

    v = v(:);
    if isscalar(v)
        v = repmat(v, A.geom.Np, 1);
    end

    xobs = xobs(:).';

    dx = A.geom.xc(:,1) - xobs(1);
    dy = A.geom.xc(:,2) - xobs(2);
    dz = A.geom.xc(:,3) - xobs(3);
    R  = sqrt(dx.^2 + dy.^2 + dz.^2);

    G = A.green_scale * exp(-1i * A.k(fi) * R) ./ (4*pi*max(R, eps));

    % Source area weighting is included here
    G = G .* A.geom.A;

    pobs = 1i * A.w(fi) * A.env.rho * sum(G .* v);
end

function G = source_kernel(A, fi)
% Source-to-source Green kernel matrix with panel area weighting.

    R = A.Rss;
    k = A.k(fi);

    G = A.green_scale * exp(-1i * k * max(R, eps)) ./ (4*pi*max(R, eps));

    % Multiply each source column by source panel area
    G = G .* (A.geom.A.');

    % Diagonal self-term:
    % Half-space (green_scale=2) -> ~ a_eq
    % Free-space (green_scale=1) -> ~ a_eq/2
    diag_term = 0.5 * A.green_scale * A.geom.aeq(:) .* exp(-1i * k * 0.5 * A.geom.aeq(:));

    idx = 1:size(R,1)+1:numel(R);
    G(idx) = diag_term;
end