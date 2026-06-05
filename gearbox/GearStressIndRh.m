% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function dstress=GearStressIndRh(Rs, Mh, zs, Rw, Rh, h, zi, nu1, Emod1, nu2, Emod2)
        sleeveangle=2*pi/zi;
        Mc=Mh*zs; %Mh*i, i - gear transmission ratio 
        nn=(1-nu1*nu1)/Emod1+(1-nu2*nu2)/Emod2;
        dSdRh=(0.5642*((Mc*cos(sleeveangle))/(Rh*Rs*Rw*h*nn*pi)-(Mc*(Rh-Rs)*cos(sleeveangle))/(Rh^2*Rs*Rw*h*nn*pi)))/(sqrt(2)*sqrt((Mc*(Rh-Rs)*cos(sleeveangle))/(Rh*Rs*Rw*h*nn*pi)));
        dstress=dSdRh;
end