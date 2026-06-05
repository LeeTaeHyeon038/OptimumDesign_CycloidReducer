% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyMmax(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    M_MAX=0.85;
    if m>M_MAX
       %m=(m-M_MAX)^2;
       m=(GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)...
           -GearVolume(r, e, zk, zs, q, M_MAX, N, h, Rs, Rh, zi))^2;       
    else
       m=0;
    end
    pen=m;
end