% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function dm=GearRhodm(e, zk, zs, m, q, alpha)
dm=(GearRho(e, zk, zs, m+0.001, q, alpha)-GearRho(e, zk, zs, m, q, alpha))/0.001;
end