% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function grho=GearRhoBis(e, zk, zs, m, q, alpha)
    grho=(GearRhoPrim(e, zk, zs, m, q, alpha+0.001)-GearRhoPrim(e, zk, zs, m, q, alpha))/0.001;
end