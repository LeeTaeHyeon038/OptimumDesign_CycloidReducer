% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function vol=GearVolumedh(R, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    vol=(GearVolume(R, e, zk, zs, q, m, N, h+0.001, Rs, Rh, zi) ...
        -GearVolume(R, e, zk, zs, q, m, N, h, Rs, Rh, zi))/0.001;
end