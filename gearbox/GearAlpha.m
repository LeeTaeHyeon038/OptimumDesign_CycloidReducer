% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018

function alpharet=GearAlpha(e, zk, zs, m, q, initalpha)
    prev=initalpha-GearRhoPrim(e, zk, zs, m, q, initalpha)/GearRhoBis(e, zk, zs, m, q, initalpha);
    next=prev-GearRhoPrim(e, zk, zs, m, q, prev)/GearRhoBis(e, zk, zs, m, q, prev);
    while(abs(prev-next)>1e-6)
        prev=next-GearRhoPrim(e, zk, zs, m, q, next)/GearRhoBis(e, zk, zs, m, q, next);
        next=prev-GearRhoPrim(e, zk, zs, m, q, prev)/GearRhoBis(e, zk, zs, m, q, prev);
    end
    alpharet=next;    
end