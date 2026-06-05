% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function n=GearNewton(e, zk, zs, m, alpha)
    prev=alpha-GearRho(e, zk, zs, m, alpha)/GearSymbolicsDiff(e, zk, zs, m, alpha);
    next=prev-GearRho(e, zk, zs, m, prev)/GearSymbolicsDiff(e, zk, zs, m, prev);
    while(abs(prev-next)>1e-6)
        prev=next-GearRho(e, zk, zs, m, next)/GearSymbolicsDiff(e, zk, zs, m, next);
        next=prev-GearRho(e, zk, zs, m, prev)/GearSymbolicsDiff(e, zk, zs, m, prev);
    end
    n=next;
end