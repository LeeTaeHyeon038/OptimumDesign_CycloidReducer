% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function vol=GearVolumedq(R, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    vol=(GearVolume(R, e, zk, zs, q+0.001, m, N, h, Rs, Rh, zi) ...
        -GearVolume(R, e, zk, zs, q, m, N, h, Rs, Rh, zi))/0.001;
end