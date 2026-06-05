% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2)   
    MAX_CONTACT_STRESS=GearGetConstant(1);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
    else
        min_alpha=alp+pi/zs;
    end    
    gstressexnegpen=GearStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, min_alpha);
    if gstressexnegpen > MAX_CONTACT_STRESS
        gstressexnegpen=(gstressexnegpen-MAX_CONTACT_STRESS)^2;
    else
        gstressexnegpen=0;
    end 
    pen=gstressexnegpen;
end