% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRhonMaxJac(e, zk, zs, m, q)
    NMAX_RHO=GearGetConstant(5);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
    else
        min_alpha=alp+pi/zs;
    end  
    nmaxrho=GearRho(e, zk, zs, m, q, min_alpha);
    if nmaxrho<NMAX_RHO
        nmaxrho=2*(NMAX_RHO-nmaxrho);
    else
        nmaxrho=0;
    end 
    pen=nmaxrho;
end