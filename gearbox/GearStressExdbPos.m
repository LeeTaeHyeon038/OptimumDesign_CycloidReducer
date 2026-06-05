% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressExdbPos(e, zk, zs, m, q, h, Mh, b, nu1, Emod1, nu2, Emod2, alpha)
        StressExdbPos=(GearStressExPos(e, zk, zs, m, q, h, Mh, b+0.001, nu1, Emod1, nu2, Emod2, alpha)...
        -GearStressExPos(e, zk, zs, m, q, h, Mh, b, nu1, Emod1, nu2, Emod2, alpha))/0.001;
        sp=StressExdbPos;
end