% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
function cnst=GearGetConstant(n)
    switch(n)
        case 1 
            cnst=400; %MAX_CONTACT_STRESS
        case 2 
            cnst=9; %PMIN_RHO
        case 3
            cnst=100; %PMAX_RHO
        case 4
            cnst=-2; %NMIN_RHO
        case 5
            cnst=-100; %NMAX_RHO
        case 6
            cnst=1000; %MIN_VOL
        otherwise
            cnst=-1;
    end
end
