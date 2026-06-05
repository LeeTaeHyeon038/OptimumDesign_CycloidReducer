% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
x0=[0.9, 0.8]';
xopt=GearSimpleLevMar(3000, @GearRosenbrock, @GearRosenbrockJac, x0);
GearValuePrev=GearRosenbrock(x0);
GearValueNext=GearRosenbrock(xopt);