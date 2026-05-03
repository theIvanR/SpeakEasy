function geom = build_square_mesh(Lx, Ly, Nx, Ny)
%BUILD_SQUARE_MESH  Uniform rectangular panel mesh in the z=0 plane.

    if nargin < 4 || isempty(Ny)
        Ny = Nx;
    end

    x_edges = linspace(-Lx/2, Lx/2, Nx + 1);
    y_edges = linspace(-Ly/2, Ly/2, Ny + 1);

    Np  = Nx * Ny;
    xc  = zeros(Np, 3);
    A   = zeros(Np, 1);
    aeq = zeros(Np, 1);

    idx = 1;
    for ix = 1:Nx
        x1 = x_edges(ix);
        x2 = x_edges(ix + 1);
        xm = 0.5 * (x1 + x2);

        for iy = 1:Ny
            y1 = y_edges(iy);
            y2 = y_edges(iy + 1);
            ym = 0.5 * (y1 + y2);

            Ai = (x2 - x1) * (y2 - y1);

            xc(idx, :) = [xm, ym, 0];
            A(idx)     = Ai;
            aeq(idx)   = sqrt(Ai/pi);

            idx = idx + 1;
        end
    end

    geom = struct();
    geom.xc  = xc;
    geom.A   = A;
    geom.aeq = aeq;
    geom.Np  = Np;
end