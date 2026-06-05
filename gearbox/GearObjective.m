% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function [vsobj]=GearObjective(x)
    zs=15;
    zk=16;
    r=x(1)*(zs+1)/x(3);
    e=x(1);
    q=x(2);
    m=x(3);
    N=1000;
    h=x(4);
    Rs=x(5);
    Rh=x(6);
    zi=8;
    Mh=1550;
    nu1=0.3;
    nu2=0.3;
    Emod1=200000;
    Emod2=200000;
    Rw=x(7);
    
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
        max_alpha=alp+pi/zs;
    else
        max_alpha=alp;
        min_alpha=alp+pi/zs;
    end    
    
    vsobj=[GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, max_alpha);
           GearStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, min_alpha);
           GearStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2);
           GearPenaltyVmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyEmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyMmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyMmax(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyHmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyRsmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyRhmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyqmin(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);
           GearPenaltyRwmin(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2);
           GearPenaltyRwmax(r, zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2);
           GearPenaltyRhopMin(e, zk, zs, m, q);
           GearPenaltyRhopMax(e, zk, zs, m, q);
           GearPenaltyRhonMin(e, zk, zs, m, q);
           GearPenaltyRhonMax(e, zk, zs, m, q);
           GearPenaltyStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2);
           GearPenaltyStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2);
           GearPenaltyStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2);
           ];
end