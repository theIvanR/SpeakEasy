% ============================================================
% ELECTROACOUSTIC MODEL CONTRACT (READ FIRST)
% ============================================================
%
% This solver operates under the following physical assumptions:
%
% 1. ACOUSTIC MEDIUM
%    - Linear, inviscid, isentropic fluid (Helmholtz regime)
%    - Constant sound speed c and density rho
%    - No thermal or viscous boundary layer losses
%
% 2. RADIATION MODEL
%    - Surface Boundary Element Method (BEM)
%    - Arbitrary geometry supported via discretized mesh:
%        geom = {xc, A, aeq, Np}
%    - Radiation is fully geometry-dependent
%
% 3. STRUCTURAL MODEL
%    - Single-degree-of-freedom rigid piston (current version)
%    - No cone breakup modes included (future extension planned)
%
% 4. ENCLOSURE MODEL
%    - Lumped mechanical impedance ONLY
%    - IB (infinite baffle) = Zbox = 0
%    - Sealed box = compliance + optional damping
%    - No spatial cavity modes included
%
% 5. COUPLING
%    - Electromechanical coupling via Bl product
%    - Acoustic loading enters as mechanical impedance
%
% 6. IMPORTANT DESIGN RULE
%    - Geometry affects ONLY radiation operators
%    - Enclosure affects ONLY lumped mechanical impedance
%    - No hidden coupling paths are allowed
%
% ============================================================