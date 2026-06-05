% Cycloidal Gearbox Optimization Toolkit
% (C) Roman Król
% University of Technology and Humanities in Radom, Poland
% December 2018
xp=[3.0, 3.3, 0.56, 10, 5, 8, 40]';
%xfin=GearLevMar(3000, @GearObjective, @GearJac, xp);
xfin=GearSimpleLevMar(3000, @GearObjective, @GearJac, xp);
GearDraw(xfin, xp);
VolumeFin=GearVolume(xfin(1)*(15+1)/xfin(3), xfin(1), 16, 15, xfin(2), xfin(3), 1000, xfin(4), xfin(5), xfin(6), 8);
VolumeXp=GearVolume(xp(1)*(15+1)/xp(3), xp(1), 16, 15, xp(2), xp(3), 1000, xp(4), xp(5), xp(6), 8);
% R=xfin(1)*(15+1)/xfin(3);
% e=xfin(1);
% q=xfin(2);
% m=xfin(3);
% h=xfin(4);
% Rs=xfin(5);
% Rh=xfin(6);
% Rw=xfin(7);
    alp=GearAlpha(xp(1), 16, 15, xp(3), xp(2), 0);
    rho=GearRho(xp(1), 16, 15, xp(3), xp(2), alp);
    if rho<0
        min_alpha=alp;
        max_alpha=alp+pi/15;
    else
        max_alpha=alp;
        min_alpha=alp+pi/15;
    end  
StressExPosPrev=GearStressExPos(xp(1), 16, 15, xp(3), xp(2), xp(4), 1550, 0.3, 200000, 0.3, 200000, max_alpha);
StressExNegPrev=GearStressExNeg(xp(1), 16, 15, xp(3), xp(2), xp(4), 1550, 0.3, 200000, 0.3, 200000, min_alpha);
StressInPrev=GearStressIn(15, 8, xp(4), 1550, xp(7), xp(5), xp(6), 0.3, 200000, 0.3, 200000);
VolumePrev=GearVolume(xp(1)*16/xp(3), xp(1), 16, 15, xp(2), xp(3), 1000, xp(4), xp(5), xp(6), 8);

    alp=GearAlpha(xfin(1), 16, 15, xfin(3), xfin(2), 0);
    rho=GearRho(xfin(1), 16, 15, xfin(3), xfin(2), alp);
    if rho<0
        min_alpha=alp;
        max_alpha=alp+pi/15;
    else
        max_alpha=alp;
        min_alpha=alp+pi/15;
    end  
StressExPosNext=GearStressExPos(xfin(1), 16, 15, xfin(3), xfin(2), xfin(4), 1550, 0.3, 200000, 0.3, 200000, max_alpha);
StressExNegNext=GearStressExNeg(xfin(1), 16, 15, xfin(3), xfin(2), xfin(4), 1550, 0.3, 200000, 0.3, 200000, min_alpha);
StressInNext=GearStressIn(15, 8, xfin(4), 1550, xfin(7), xfin(5), xfin(6), 0.3, 200000, 0.3, 200000);
VolumeNext=GearVolume(xfin(1)*16/xfin(3), xfin(1), 16, 15, xfin(2), xfin(3), 1000, xfin(4), xfin(5), xfin(6), 8);

xfin2=GearSteepestDescent(10000, @GearObjective, @GearJac, xp);
figure;
GearDraw(xfin2, xp);
VolumeFin2=GearVolume(xfin2(1)*(15+1)/xfin2(3), xfin2(1), 16, 15, xfin2(2), xfin2(3), 1000, xfin2(4), xfin2(5), xfin2(6), 8);
VolumeXp2=GearVolume(xp(1)*(15+1)/xp(3), xp(1), 16, 15, xp(2), xp(3), 1000, xp(4), xp(5), xp(6), 8);
% R=xfin(1)*(15+1)/xfin(3);
% e=xfin(1);
% q=xfin(2);
% m=xfin(3);
% h=xfin(4);
% Rs=xfin(5);
% Rh=xfin(6);
% Rw=xfin(7);
    alp=GearAlpha(xp(1), 16, 15, xp(3), xp(2), 0);
    rho=GearRho(xp(1), 16, 15, xp(3), xp(2), alp);
    if rho<0
        min_alpha=alp;
        max_alpha=alp+pi/15;
    else
        max_alpha=alp;
        min_alpha=alp+pi/15;
    end  
StressExPosPrev2=GearStressExPos(xp(1), 16, 15, xp(3), xp(2), xp(4), 1550, 0.3, 200000, 0.3, 200000, max_alpha);
StressExNegPrev2=GearStressExNeg(xp(1), 16, 15, xp(3), xp(2), xp(4), 1550, 0.3, 200000, 0.3, 200000, min_alpha);
StressInPrev2=GearStressIn(15, 8, xp(4), 1550, xp(7), xp(5), xp(6), 0.3, 200000, 0.3, 200000);
VolumePrev2=GearVolume(xp(1)*16/xp(3), xp(1), 16, 15, xp(2), xp(3), 1000, xp(4), xp(5), xp(6), 8);

    alp=GearAlpha(xfin2(1), 16, 15, xfin2(3), xfin2(2), 0);
    rho=GearRho(xfin2(1), 16, 15, xfin2(3), xfin2(2), alp);
    if rho<0
        min_alpha=alp;
        max_alpha=alp+pi/15;
    else
        max_alpha=alp;
        min_alpha=alp+pi/15;
    end  
StressExPosNext2=GearStressExPos(xfin2(1), 16, 15, xfin2(3), xfin2(2), xfin2(4), 1550, 0.3, 200000, 0.3, 200000, max_alpha);
StressExNegNext2=GearStressExNeg(xfin2(1), 16, 15, xfin2(3), xfin2(2), xfin2(4), 1550, 0.3, 200000, 0.3, 200000, min_alpha);
StressInNext2=GearStressIn(15, 8, xfin2(4), 1550, xfin2(7), xfin2(5), xfin2(6), 0.3, 200000, 0.3, 200000);
VolumeNext2=GearVolume(xfin2(1)*16/xfin2(3), xfin2(1), 16, 15, xfin2(2), xfin2(3), 1000, xfin2(4), xfin2(5), xfin2(6), 8);

