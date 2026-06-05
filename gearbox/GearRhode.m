% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function de=GearRhode(e, zk, zs, m, q, alpha)
de=(GearRho(e+0.001, zk, zs, m, q, alpha)-GearRho(e, zk, zs, m, q, alpha))/0.001;
end