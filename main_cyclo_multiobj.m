%% main_cyclo_multiobj.m
% Multi-objective version: volume vs. inefficiency Pareto front via gamultiobj,
% using the SAME model functions as the single-objective study. This is the step
% the Golinski benchmark never takes (it is single-objective weight only), but
% the cycloid problem is naturally two-objective, so the trade-off matters.
%
% Requires: Global Optimization Toolbox (gamultiobj).
clear; clc; close all;

p     = cyclo_params();
lb    = p.lb;  ub = p.ub;  nvars = numel(lb);

% objective vector: [ volume ; 1 - efficiency ]  (both minimised)
fmulti = @(x) [ cyclo_volume(x,p), 1 - cyclo_efficiency(x,p) ];
con    = @(x) cyclo_nonlcon(x,p);

opt = optimoptions('gamultiobj','Display','iter','PopulationSize',150, ...
                   'PlotFcn',@gaplotpareto);
[X,F] = gamultiobj(fmulti,nvars,[],[],[],[],lb,ub,con,opt);

% F(:,1) = volume, F(:,2) = 1 - eta  ->  plot volume vs. efficiency
figure;
plot(F(:,1), 1 - F(:,2), 'o','MarkerFaceColor',[0 0.45 0.74]); grid on;
xlabel('Volume  [mm^3]'); ylabel('Efficiency  \eta');
title('Pareto front: volume vs. efficiency (cycloid reducer)');

% X holds the design vector [Dp drp B D K1 Dw dsw] for each Pareto point;
% pick the knee point (or your preferred trade-off) as the final design.
