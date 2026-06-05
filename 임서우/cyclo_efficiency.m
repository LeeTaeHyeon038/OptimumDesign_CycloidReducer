function eta = cyclo_efficiency(x, p)
% CYCLO_EFFICIENCY  Total efficiency of the cycloid reducer (Wang Eq.6).
%   x = [Dp drp B D K1 Dw dsw]
    Dp = x(1); drp = x(2); dsw = x(7); K1 = x(5); Dw = x(6);

    % meshing efficiency, pin wheel <-> cycloid gear (Eq.4)
    a = (Dp - drp)*4*p.mu/(K1*p.zc*Dp*pi);
    b = (Dp - drp)*4*p.mu/(K1*Dp*pi);
    eta_x = (1 - a)/(1 + b);

    % efficiency of the output part (Eq.5)
    eta_sx = 1 - 4*p.fw*K1*dsw*Dp/(pi*Dw*(dsw + 2*p.Delta2));

    % total (Eq.6): meshing * output * arm-bearing * (rolling-bearing)^2
    eta = eta_x*eta_sx*p.eta_zx*p.eta_gx^2;
end
