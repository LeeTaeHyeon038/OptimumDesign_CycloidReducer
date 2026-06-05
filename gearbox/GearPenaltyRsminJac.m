% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRsminJac(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    RS_MIN=3;
    if Rs<RS_MIN
       %Rs=2*(RS_MIN-Rs);
       Rs=2*(GearVolume(r, e, zk, zs, q, m, N, h, RS_MIN, Rh, zi) ...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)) ...
           *GearVolumedRs(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);       
    else
       Rs=0;
    end
    pen=Rs;
end