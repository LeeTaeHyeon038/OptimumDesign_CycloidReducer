% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhopMax(e, zk, zs, m, q)
    PMAX_RHO=GearGetConstant(3);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        max_alpha=alp+pi/zs;
    else
        max_alpha=alp;
    end  
    pmaxrho=GearRho(e, zk, zs, m, q, max_alpha);
    if pmaxrho>PMAX_RHO
        pmaxrho=(abs(pmaxrho)-PMAX_RHO)^2;
    else
        pmaxrho=0;
    end 
    pen=pmaxrho;
end