% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2)   
    MAX_CONTACT_STRESS=GearGetConstant(1);
    gstressinpen=GearStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2);
    if gstressinpen > MAX_CONTACT_STRESS
        gstressinpen=(gstressinpen-MAX_CONTACT_STRESS)^2;
    else
        gstressinpen=0;
    end
    pen=gstressinpen;
end