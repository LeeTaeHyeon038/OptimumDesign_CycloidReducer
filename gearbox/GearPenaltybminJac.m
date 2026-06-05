% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltybminJac(e, zk, zs, m, q, b)
    alp=GearAlpha(e, zk, zs, m, q, 0);
    rho=GearRho(e, zk, zs, m, q, alp);
    if rho<0
        min_alpha=alp;
    else
        min_alpha=alp+pi/zs;
    end      
    if b > abs(GearRho(e, zk, zs, m, q, min_alpha))
       b=2*(b-abs(GearRho(e, zk, zs, m, q, min_alpha)));
    else
       b=0;
    end
    pen=b;
end