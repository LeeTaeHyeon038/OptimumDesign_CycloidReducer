% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function vol=GearVolumedRs(R, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    vol=(GearVolume(R, e, zk, zs, q, m, N, h, Rs+0.001, Rh, zi) ...
        -GearVolume(R, e, zk, zs, q, m, N, h, Rs, Rh, zi))/0.001;
end