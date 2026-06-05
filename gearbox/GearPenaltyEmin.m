% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyEmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    E_MIN=q*(zk+1)/(3*sqrt(3)*zk)*sqrt((zk+1)/(zk-1))*sqrt(m^2/(1-m^2));
    if e<E_MIN
       %e=(E_MIN-e)^2;
       e=(GearVolume(r, E_MIN, zk, zs, q, m, N, h, Rs, Rh, zi)...
           -GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi))^2;
    else
       e=0;
    end
    pen=e;
end