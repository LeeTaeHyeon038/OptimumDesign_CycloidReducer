% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyHmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    H_MIN=0.2;
    if h<H_MIN
       %h=(H_MIN-h)^2;
       h=(GearVolume(r, e, zk, zs, q, m, N, H_MIN, Rs, Rh, zi)...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi))^2;
    else
       h=0;
    end
    pen=h;
end