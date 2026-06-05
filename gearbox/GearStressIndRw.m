% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function dstress=GearStressIndRw(Rs, Mh, zs, Rw, Rh, h, zi, nu1, Emod1, nu2, Emod2)
        sleeveangle=2*pi/zi;
        Mc=Mh*zs; %Mh*i, i - gear transmission ratio 
        nn=(1-nu1*nu1)/Emod1+(1-nu2*nu2)/Emod2;
        dSdRw=-(0.5642*Mc*(Rh-Rs)*cos(sleeveangle))...
            /(sqrt(2)*Rh*Rs*Rw^2*h*nn*pi*sqrt((Mc*(Rh-Rs)...
            *cos(sleeveangle))/(Rh*Rs*Rw*h*nn*pi)));
        dstress=dSdRw;
end