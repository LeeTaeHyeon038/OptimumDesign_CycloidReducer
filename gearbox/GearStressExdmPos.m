% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressExdmPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha)
    StressExdmPos=(GearStressExPos(e, zk, zs, m+0.001, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha)...
        -GearStressExPos(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha))/0.001;
    sp=StressExdmPos;
end