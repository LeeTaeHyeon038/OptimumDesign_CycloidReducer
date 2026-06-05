% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyVmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    MIN_VOL=GearGetConstant(6);
    vol=GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
    if vol < MIN_VOL
            vol=(MIN_VOL-vol)^2;
    else
            vol=0;
    end
    pen=vol;
end