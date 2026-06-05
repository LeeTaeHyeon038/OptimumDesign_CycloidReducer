% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyMmaxJac(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    M_MAX=0.85;
    if m>M_MAX
       %m=2*(m-M_MAX);
       m=2*(GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi) ...
           -GearVolume(r, e, zk, zs, q, M_MAX, N, h, Rs, Rh, zi)) ...
           *GearVolumedm(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);       
    else
       m=0;
    end
    pen=m;
end