% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function sp=GearStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2)
        sleeveangle=2*pi/zi;
        Mc=Mh*zs; %Mh*i, i - gear transmission ratio
        insleeveforce=2*Mc*(cos(sleeveangle))/(pi*Rw);
        nn=(1-nu1*nu1)/Emod1+(1-nu2*nu2)/Emod2;
        MaxStressIn=0.5642*sqrt((insleeveforce/h)*(Rh-Rs)/(nn*Rh*Rs));
        sp=MaxStressIn;
end