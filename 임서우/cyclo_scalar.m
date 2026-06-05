function f = cyclo_scalar(x, p, w)
% CYCLO_SCALAR  Weighted single objective so the Golinski-style single-objective
% solvers (fmincon/ga/particleswarm/simulannealbnd) can be reused directly.
%   w = 1  -> minimise inefficiency only (maximise efficiency)
%   w = 0  -> minimise volume only
% Both terms are normalised to ~O(1) at the baseline.
    eta = cyclo_efficiency(x, p);
    V   = cyclo_volume(x, p);
    f = w*(1 - eta)/(1 - p.eta0) + (1 - w)*V/p.V0;
end
