% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function vol=GearVolume(R, e, zk, zs, q, m, N, h, Rs, Rh, zi)
    dt=2*pi/N;
    vol=0;
    for i=0:(N-1)
        alpha=i*dt+dt/2;
        A=abs(q*cos(atan(sin(alpha*zs)...
            /(cos(alpha*zs)+1/m))+alpha)-e*cos(alpha*zk)...
            -R*cos(alpha))*abs(q*(((zs*sin(alpha*zs)^2)...
            /(cos(alpha*zs)+1/m)^2+(zs*cos(alpha*zs))/(cos(alpha*zs)+1/m))...
            /(sin(alpha*zs)^2/(cos(alpha*zs)+1/m)^2+1)+1)*cos(atan(sin(alpha*zs)...
            /(cos(alpha*zs)+1/m))+alpha)-e*zk*cos(alpha*zk)-R*cos(alpha));    
        vol=vol+A*dt;
    end
    vol=vol*h-pi*Rh*Rh*zi*h+pi*q*q*zk*h+pi*Rs*Rs*zi*h;
end