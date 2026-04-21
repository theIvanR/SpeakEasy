function Ze = simulate_speaker(f, p, opts)

% --- unpack ---
w   = 2*pi*f(:);

Re  = p.Re;
Le  = p.Le;
Bl  = p.Bl;
Rms = p.Rms;
Mms = p.Mms;
Cms = p.Cms;
Sd  = p.Sd;
rho = p.rho;
c   = p.c;

% --- leakage / measurement chain ---
R_leak = 0;
L_leak = 0;
if isfield(opts, 'leakage')
    if isfield(opts.leakage, 'R') && ~isempty(opts.leakage.R)
        R_leak = opts.leakage.R;
    end
    if isfield(opts.leakage, 'L') && ~isempty(opts.leakage.L)
        L_leak = opts.leakage.L;
    end
end

% --- radiation radius ---
a = sqrt(Sd/pi);

% --- radiation impedance (external function) ---
[R_rad, X_rad] = pistonBEM(a, f, rho, c, 100);
R_rad = R_rad(:);
X_rad = X_rad(:);

% --- added mass ---
Ma = X_rad ./ (w + eps);

% --- enclosure ---
R_box = 0;

if isfield(opts, 'box') && opts.box
    if ~isfield(opts, 'alpha')
        error('opts.alpha must be provided when opts.box is true.');
    end
    if ~isfield(opts, 'Vb')
        error('opts.Vb must be provided when opts.box is true.');
    end

    alpha = opts.alpha;
    Vb    = opts.Vb;

    % Mechanical compliance of sealed box seen by the cone
    C_box = Vb / (rho * c^2 * Sd^2);

    % Simple box loss term
    R_box = alpha * (rho * c * Sd^2 / Vb);

    % Suspension compliance in parallel with box compliance
    Cms = (Cms * C_box) / (Cms + C_box);
end

% --- mechanical impedance ---
Zm = (Rms + R_rad + R_box) + 1i .* w .* (Mms + Ma) + 1 ./ (1i .* w .* Cms);

% --- electrical impedance of the driver core ---
Ze_driver = Re + 1i .* w .* Le + (Bl^2) ./ Zm;

% --- total measured impedance including leakage ---
Ze = Ze_driver + R_leak + 1i .* w .* L_leak;

end