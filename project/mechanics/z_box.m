function Z = z_box(s, p, env)
%Z_BOX  Enclosure impedance model (mechanical domain, strict mode)
%
% Supports:
%   - IB (infinite baffle)
%   - sealed box
%
% Future:
%   - ported (Helmholtz)
%   - CAD (full wave / BEM / modal)

    Z = zeros(size(s));  % default: IB / no enclosure

    % ---- Guard: env + enc ----
    if ~isstruct(env) || ~isfield(env, 'enc') || isempty(env.enc)
        return;
    end

    enc = env.enc;

    if ~isfield(enc, 'type') || isempty(enc.type)
        return;
    end

    type = lower(string(enc.type));

    switch type

        % =========================================================
        % 1) INFINITE BAFFLE (IB)
        % =========================================================
        case {"ib","infinitebaffle","none","open"}

            assert_no_fields(enc, {'Vb','Rleak','Rloss','Fb','port'}, ...
                'IB must be parameter-free');

            Z = zeros(size(s));
            return;

        % =========================================================
        % 2) SEALED ENCLOSURE
        % =========================================================
        case "sealed"
        
            % ---- HARD REQUIREMENTS (no defaults allowed) ----
            required = {'Vb','Rloss','Rleak'};
            for i = 1:numel(required)
                f = required{i};
                if ~isfield(enc, f) || isempty(enc.(f))
                    error("sealed enclosure requires enc.%s (no defaults allowed)", f);
                end
            end
        
            % ---- constants (allowed to come from env) ----
            rho = env.rho;
            c   = env.c;
        
            % ---- physics ----
            Cab = enc.Vb / (rho * c^2);
        
            % total acoustic loss MUST be explicitly provided
            Rac = enc.Rloss + enc.Rleak;
        
            % acoustic impedance
            Zac = 1 ./ (s * Cab) + Rac;
        
            % mechanical transform
            Z = (p.Sd^2) .* Zac;
        
            return;

        % =========================================================
        % 3) PORTED (not implemented)
        % =========================================================
        case "ported"
            error("PORTED enclosure not implemented (needs Helmholtz + radiation coupling)");

        % =========================================================
        % 4) CAD / MODAL (not implemented)
        % =========================================================
        case "cad"
            error("CAD enclosure not implemented (requires wave/BEM modal solver)");
        
        % Something like this and we make another folder for FEM
        % if ~isfield(enc, 'fem')
        %     error("CAD enclosure missing compiled FEM model");
        % end
        % 
        % fem = enc.fem;
        % 
        % Z = fem_eval(fem, s);

        % =========================================================
        % UNKNOWN TYPE
        % =========================================================
        otherwise
            error("Unknown enclosure type: %s", type);
    end
end

% ================================================================
% Helper: assert forbidden fields absent
% ================================================================
function assert_no_fields(s, fields, msg)
    for i = 1:numel(fields)
        if isfield(s, fields{i}) && ~isempty(s.(fields{i}))
            error("%s (forbidden field: %s)", msg, fields{i});
        end
    end
end