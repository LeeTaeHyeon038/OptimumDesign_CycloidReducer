% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressExdqNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha)
    StressExdqNeg=(GearStressExNeg(e, zk, zs, m, q+0.001, h, Mh, nu1, Emod1, nu2, Emod2, alpha)...
        -GearStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha))/0.001;
    sp=StressExdqNeg;
end