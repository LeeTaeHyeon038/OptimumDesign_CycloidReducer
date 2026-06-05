% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyqminJac(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
    else
        min_alpha=alp+pi/zs;
    end      
    if q > abs(GearRho(e, zk, zs, m, q, min_alpha))
       %q=2*(q-abs(GearRho(e, zk, zs, m, q, min_alpha)));
       q=2*(GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi) ...
           -GearVolume(r, e, zk, zs, abs(GearRho(e, zk, zs, m, q, min_alpha)), m, N, h, Rs, Rh, zi)) ...
           *GearVolumedq(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);       
    else
       q=0;
    end
    pen=q;
end