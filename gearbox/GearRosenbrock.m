% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function gr=GearRosenbrock(x)
    gr=(1-x(1))^2+100*(x(2)-x(1)^2)^2;
end