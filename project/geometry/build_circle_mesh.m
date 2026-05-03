function geom = build_circle_mesh(Sd, Nrad, Mtheta)
%BUILD_CIRCLE_MESH  Polar panel mesh for a circular piston.

    a = sqrt(Sd / pi);

    r_edges = linspace(0, a, Nrad + 1);
    t_edges = linspace(0, 2*pi, Mtheta + 1);

    Np  = Nrad * Mtheta;
    xc  = zeros(Np, 3);
    A   = zeros(Np, 1);
    aeq = zeros(Np, 1);

    idx = 1;
    for it = 1:Mtheta
        t1 = t_edges(it);
        t2 = t_edges(it + 1);
        tm = 0.5 * (t1 + t2);

        for ir = 1:Nrad
            r1 = r_edges(ir);
            r2 = r_edges(ir + 1);

            Ai = 0.5 * (r2^2 - r1^2) * (t2 - t1);

            if abs(r2^2 - r1^2) < eps
                rc = 0;
            else
                rc = (2/3) * (r2^3 - r1^3) / (r2^2 - r1^2);
            end

            xc(idx, :) = [rc*cos(tm), rc*sin(tm), 0];
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