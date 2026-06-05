function p = cyclo_params()
% CYCLO_PARAMS  Fixed constants for the K-H-V cycloid speed reducer model.
% Baseline values follow the XW3-type example in Wang, Luo & Su (2016).
% Items marked "-- verify" depend on your material / bearing / geometry and
% should be set to match your actual design before trusting the optimum.

    % ---- application (FIXED, not design variables) ----
    p.i   = 43;             % transmission ratio
    p.zc  = 43;             % cycloid gear teeth  (i = zc)
    p.zp  = p.zc + 1;       % pin teeth           (zp = zc + 1)
    p.zw  = 8;              % number of output pins
    p.P   = 0.75;           % input power [kW]
    p.n   = 1440;           % input speed [rpm]
    p.M   = 9550*p.P/(p.n/p.i)*1e3;   % output torque [N.mm]

    % ---- material / contact ----
    p.Ee       = 2.06e5;    % equivalent elastic modulus [MPa] (GCr15) -- verify
    p.sHP      = 1000;      % allowable contact stress, gear  [MPa]
    p.sFP      = 150;       % allowable bending stress, pin wheel [MPa]
    p.sBP      = 200;       % allowable bending stress, pin [MPa]
    p.sZHP     = 1000;      % allowable contact stress, pin/hole [MPa] -- verify
    p.rho_emin = 5.0;       % min equivalent curvature radius [mm] -- PLACEHOLDER
                            %   (really a function of the profile; see notes)

    % ---- friction / efficiency ----
    p.mu     = 0.06;        % pin<->cycloid friction   (0.05..0.10)
    p.fw     = 0.04;        % output-part friction     (0.008..0.08)
    p.eta_zx = 0.99;        % pivoted-arm bearing efficiency
    p.eta_gx = 0.995;       % rolling bearing efficiency

    % ---- geometry helpers ----
    p.Delta2 = 2.0;         % pin-sleeve wall thickness [mm] -- verify
    p.Deltac = 1.0;         % space-ring thickness [mm]
    p.Kw     = 1.4;         % load distribution factor (pin bending)
    p.L      = 30;          % pin effective length [mm] -- verify

    % ---- coefficient / size ranges (handbook tables, zc=43, zp=44) ----
    p.K1_lo = 0.65; p.K1_hi = 0.90;   % short width coefficient
    p.K2_lo = 1.0;  p.K2_hi = 1.6;    % pin-diameter coefficient
    p.Dp_lo = 140;  p.Dp_hi = 155;    % pin centre-circle diameter [mm]

    % ---- pivoted-arm bearing life ----
    p.C_bearing = 3.0e4;    % rated dynamic load [N] -- from bearing catalog
    p.Lh_min    = 5000;     % required life [h]

    % ---- box bounds (Wang et al. Table 4): X = [Dp drp B D K1 Dw dsw] ----
    p.lb = [140  7.0   7   50  0.65   88  11];
    p.ub = [155 10.4  12   55  0.90  104  14];

    % ---- normalisation baselines for the weighted scalar objective ----
    p.V0   = 1.0e5;         % ~ baseline volume [mm^3]
    p.eta0 = 0.86;          % ~ baseline efficiency
end
