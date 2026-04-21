function [Rrad, Xrad] = pistonBEM(a, f, rho, c, N)
% PISTONBEM  Fast ring‐BEM radiation impedance of a baffled piston
%   [Rrad,Xrad] = pistonBEM(a,f,rho,c,N)
%   N = # rings (default 100).
    if nargin<5, N = 100; end

    % Constants & discretization
    omega = 2*pi*f(:).';         % 1×M
    k     = omega / c;           % 1×M
    dr    = a/N;
    r     = (dr/2 : dr : a-dr/2)';  % N×1 ring radii
    A     = 2*pi*r*dr;              % N×1 ring element areas
    [Ri, Rj] = meshgrid(r, r);      % N×N
    AiAj = A * A.';                 % N×N

    % Theta quadrature (uniform)
    M     = 32;
    thetas = linspace(2*pi/M, 2*pi-2*pi/M, M);
    cosT   = cos(thetas);

    % Pre‐allocate
    Zrad = zeros(1, numel(f));

    % Loop only over frequency
    for fi = 1:numel(f)
        Ki = k(fi);

        % Build Gij(θ) for all θ at once: N×N×M
        %   Rij(:,:,m) = sqrt(Ri.^2 + Rj.^2 - 2 Ri.*Rj*cosT(m))
        Rij = sqrt( Ri.^2 + Rj.^2 - 2*(Ri.*Rj).*reshape(cosT,1,1,[]) );
        Rij = max(Rij, dr*1e-3);   % floor to avoid singularity

        Gij = exp(-1i * Ki * Rij) ./ Rij;  % N×N×M

        Gavg = mean(Gij, 3);  % average over θ → N×N

        % Sum all mutual impedances in one dot‐product INTO HALF SPACE!!!
        Zrad(fi) = 2 * (1i*omega(fi)*rho/(4*pi)) * sum( AiAj(:) .* Gavg(:) ); %2 factor for half space
    end

    % Split real/imag
    Rrad = real(Zrad);
    Xrad = imag(Zrad);
end

