% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressExdbNeg(e, zk, zs, m, q, h, Mh, b, nu1, Emod1, nu2, Emod2, alpha)
        StressExdbNeg=(GearStressExNeg(e, zk, zs, m, q, h, Mh, b+0.001, nu1, Emod1, nu2, Emod2, alpha)...
        -GearStressExNeg(e, zk, zs, m, q, h, Mh, b, nu1, Emod1, nu2, Emod2, alpha))/0.001;
        sp=StressExdbNeg;
end