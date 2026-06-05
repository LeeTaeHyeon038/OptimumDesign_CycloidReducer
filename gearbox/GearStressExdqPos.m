% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressExdqPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha)
    StressExdqPos=(GearStressExPos(e, zk, zs, m, q+0.001, h, Mh, nu1, Emod1, nu2, Emod2, alpha)...
        -GearStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha))/0.001;
    sp=StressExdqPos;
end