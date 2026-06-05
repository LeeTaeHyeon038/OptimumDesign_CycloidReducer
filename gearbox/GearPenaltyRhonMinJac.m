% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhonMinJac(e, zk, zs, m, q)
    NMIN_RHO=GearGetConstant(4);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
    else
        min_alpha=alp+pi/zs;
    end  
    nminrho=GearRho(e, zk, zs, m, q, min_alpha);
    if nminrho>NMIN_RHO
        nminrho=2*(nminrho-NMIN_RHO);
    else
        nminrho=0;
    end 
    pen=nminrho;
end