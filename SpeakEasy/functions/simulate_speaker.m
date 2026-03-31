function Ze = simulate_speaker(f, p, opts)

% --- unpack ---
w = 2*pi*f(:);
Re = p.Re; Le = p.Le; Bl = p.Bl;
Rms = p.Rms; Mms = p.Mms; Cms = p.Cms;
Sd = p.Sd; rho = p.rho; c = p.c;

% --- configuration (isobaric etc.) ---
if ~isfield(opts,'cfg') || isempty(opts.cfg)
    cfg = 'single';
else
    cfg = opts.cfg;
end

switch cfg
    case 'single'
        Cms_eff = Cms; Sd_eff = Sd;
    case 'series'
        Cms_eff = Cms/2; Sd_eff = Sd;
    case 'parallel'
        Cms_eff = 2*Cms; Sd_eff = 2*Sd;
    case 'series_parallel'
        Cms_eff = Cms; Sd_eff = 2*Sd;
    otherwise
        error('Unknown cfg');
end

a = sqrt(Sd_eff/pi);

% --- radiation impedance (external function) ---
[R_rad, X_rad] = pistonBEM(a, f, rho, c, 100);
R_rad = R_rad(:); 
X_rad = X_rad(:);

% --- added mass ---
Ma = X_rad ./ (w + eps);

% --- enclosure ---
R_box = 0;

if isfield(opts,'box') && opts.box
    if ~isfield(opts,'alpha')
        error('opts.alpha must be provided when opts.box is true.');
    end
    if ~isfield(opts,'Vb')
        error('opts.Vb must be provided when opts.box is true.');
    end

    alpha = opts.alpha;
    Vb = opts.Vb;

    C_box = Vb / (rho * c^2 * Sd_eff^2);
    R_box = alpha * (rho * c * Sd_eff^2 / Vb);

    Cms_eff = (Cms_eff * C_box) / (Cms_eff + C_box);
end

% --- mechanical impedance ---
Zm = (Rms + R_rad + R_box) + 1i .* w .* (Mms + Ma) + 1 ./ (1i .* w .* Cms_eff);

% --- electrical impedance ---
Ze = Re + 1i .* w .* Le + (Bl^2) ./ Zm;

end