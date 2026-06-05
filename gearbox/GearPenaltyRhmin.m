% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    if Rh<Rs
       %Rh=(Rs-Rh)^2;
       Rh=(GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rs, zi)...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi))^2;       
    else
       Rh=0;
    end
    pen=Rh;
end