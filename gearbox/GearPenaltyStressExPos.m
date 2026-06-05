% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2)   
    MAX_CONTACT_STRESS=GearGetConstant(1);
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        max_alpha=alp+pi/zs;
    else
        max_alpha=alp;
    end      
    gstressexpospen=GearStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, max_alpha);
    if gstressexpospen > MAX_CONTACT_STRESS
        gstressexpospen=(gstressexpospen-MAX_CONTACT_STRESS)^2;
    else
        gstressexpospen=0;
    end
    pen=gstressexpospen;
end