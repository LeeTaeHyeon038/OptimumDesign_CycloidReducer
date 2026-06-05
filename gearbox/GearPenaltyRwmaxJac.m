% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function pen=GearPenaltyRwmaxJac(r, zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2)   
    if Rw > (r-Rs)
       %Rw=2*(Rw-(r-Rs));
       Rw=2*(GearStressIn(zs, zi, h, Mh, Rw, Rs, Rh, nu1, Emod1, nu2, Emod2) ...
           -GearStressIn(zs, zi, h, Mh, (r-Rs), Rs, Rh, nu1, Emod1, nu2, Emod2)) ...
           *GearStressIndRw(Rs, Mh, zs, Rw, Rh, h, zi, nu1, Emod1, nu2, Emod2);       
    else
       Rw=0;
    end
    pen=Rw;
end