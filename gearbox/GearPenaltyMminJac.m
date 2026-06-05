% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyMminJac(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    M_MIN=0.5;
    if m<M_MIN
       %m=2*(M_MIN-m);
       m=2*(GearVolume(r, e, zk, zs, q, M_MIN, N, h, Rs, Rh, zi) ...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)) ...
           *GearVolumedm(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);       
    else
       m=0;
    end
    pen=m;
end