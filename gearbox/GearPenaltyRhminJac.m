% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhminJac(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    if Rh<Rs
       %Rh=2*(Rs-Rh);
       Rh=2*(GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rs, zi) ...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)) ...
           *GearVolumedRh(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);       
    else
       Rh=0;
    end
    pen=Rh;
end