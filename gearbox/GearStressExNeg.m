% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018

function sp=GearStressExNeg(e, zk, zs, m, q, h, Mh, nu1, Emod1, nu2, Emod2, alpha)
    rho=GearRho(e, zk, zs, m, q, alpha);
    Mc=Mh*zs; %Mh*i, i - gear transmission ratio
    nn=(1-nu1*nu1)/Emod1+(1-nu2*nu2)/Emod2;
    exsleeveforce=4*Mc/(e*zs*zk);
    MaxStressExNeg=0.5642*sqrt((exsleeveforce/h)*(abs(rho)-q)/(nn*abs(rho)*q));
    sp=MaxStressExNeg;
end