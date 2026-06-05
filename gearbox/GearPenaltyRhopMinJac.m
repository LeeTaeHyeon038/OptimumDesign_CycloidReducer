% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhopMinJac(e, zk, zs, m, q)
    PMIN_RHO=GearGetConstant(2);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        max_alpha=alp+pi/zs;
    else
        max_alpha=alp;
    end  
    pminrho=GearRho(e, zk, zs, m, q, max_alpha);
    if pminrho<PMIN_RHO
        pminrho=2*(PMIN_RHO-abs(pminrho));
    else
        pminrho=0;
    end 
    pen=pminrho;
end