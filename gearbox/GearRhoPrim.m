% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function grho=GearRhoPrim(e, zk, zs, m, q, alpha)
    grho=(GearRho(e, zk, zs, m, q, alpha+0.001)-GearRho(e, zk, zs, m, q, alpha))/0.001;
end