% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRsmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    RS_MIN=3;
    if Rs<RS_MIN
       %Rs=(RS_MIN-Rs)^2;
       Rs=(GearVolume(r, e, zk, zs, q, m, N, h, RS_MIN, Rh, zi)...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi))^2;       
    else
       Rs=0;
    end
    pen=Rs;
end