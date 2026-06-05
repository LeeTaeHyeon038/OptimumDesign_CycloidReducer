function [c, ceq] = cyclo_nonlcon(x, p)
% CYCLO_NONLCON  16 inequality constraints (Wang et al. y1..y16), all <= 0.
% Same [c,ceq] interface that fmincon / ga / gamultiobj expect, exactly like
% the Golinski Speed_Reducer_con file -- only the equations are replaced.
%   x = [Dp drp B D K1 Dw dsw]
    Dp=x(1); drp=x(2); B=x(3); D=x(4); K1=x(5); Dw=x(6); dsw=x(7);
    zp=p.zp; zc=p.zc; zw=p.zw;

    c = zeros(16,1);

    % y1: no undercut / sharp tooth tip   ->  drp/2 < |rho0|min
    thr = (zc-1)/(2*zc+1);
    if K1 > thr
        rho0min = (Dp/2)*sqrt(27*(1-K1^2)*(zp-1)/(zp+1)^3);
    else
        rho0min = Dp*(1-K1)^2/(2*(zp*K1+1));
    end
    c(1) = drp/2 - rho0min;

    % y2,y3: short-width-coefficient range
    c(2) = p.K1_lo - K1;
    c(3) = K1 - p.K1_hi;

    % y4,y5: pin-diameter-coefficient range   K2 = (Dp/drp) sin(pi/zp)
    K2 = (Dp/drp)*sin(pi/zp);
    c(4) = p.K2_lo - K2;
    c(5) = K2 - p.K2_hi;

    % y6: contact strength of the cycloid gear (Eq.22)
    c(6) = 0.418*sqrt( p.Ee/(B*K1*zc*Dp) * (4.4*p.M/p.rho_emin) ) - p.sHP;

    % y7: bending strength of the pin wheel (Eq.25); 2-pivot if Dp<390
    if Dp < 390, k = 1.41; else, k = 0.48; end
    c(7) = k*4.4*9550*p.P*p.L/(K1*Dp*p.n*drp^2) - p.sFP;

    % y8: contact strength between pin and pin-hole (Eq.27)
    rsw   = dsw/2;
    inner = (rsw + p.Delta2)^2*zp + 0.5*K1*Dp*(rsw + p.Delta2);
    c(8)  = 300*sqrt( K1*p.M*Dp/(zw*Dw*B*inner) ) - p.sZHP;

    % y9: bending strength of the pin (Eq.29)
    c(9) = 4.4*p.Kw*p.M*(1.5*B + p.Deltac)/(0.1*zw*Dw*dsw^3) - p.sBP;

    % y10,y11: pin centre-circle diameter range
    c(10) = p.Dp_lo - Dp;
    c(11) = Dp - p.Dp_hi;

    % y12,y13: minimum wall thickness (centre hole / between pin holes)
    et = K1*Dp/zp;
    c(12) = 0.06*Dp - Dw + D + dsw + 2*p.Delta2 + et + 0.15;
    c(13) = 0.03*Dp - Dw*sin(pi/zw) + dsw + 2*p.Delta2 + et + 0.15;

    % y14,y15: cycloid gear width range  (0.05 Dp <= B <= 0.1 Dp)
    c(14) = 0.05*Dp - B;
    c(15) = B - 0.1*Dp;

    % y16: pivoted-arm bearing life >= Lh_min
    R  = 0.825*p.M*zp/(K1*Dp*zc);   % equivalent radial load [N]
    pp = 1.25*R;                    % dynamic load
    n1 = p.n*zp/(zp-1);             % arm-bearing speed [rpm]
    Lh = (1e6/(60*n1))*(p.C_bearing/pp)^(10/3);
    c(16) = p.Lh_min - Lh;

    ceq = [];                       % no equality constraints
end
