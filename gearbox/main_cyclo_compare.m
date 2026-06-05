%% main_cyclo_compare.m
% Single-objective optimisation of a K-H-V cycloid speed reducer, reusing the
% Golinski-style code skeleton: obj handle + nonlcon + per-solver driver call.
% Compares fmincon (SQP), ga, particleswarm, simulannealbnd, exactly like the
% Aktemur & Gusseinov Golinski study -- only the physics (cyclo_* files) differ.
%
% Requires: Optimization Toolbox (fmincon) and Global Optimization Toolbox
%           (ga, particleswarm, simulannealbnd).
clear; clc; close all;

p     = cyclo_params();
lb    = p.lb;  ub = p.ub;  nvars = numel(lb);
w     = 0.5;          % objective weight: 1 = efficiency only, 0 = volume only
Pen   = 1e4;          % penalty factor for bounds-only solvers (PSO, SA)

obj  = @(x) cyclo_scalar(x, p, w);            % fmincon / ga objective
con  = @(x) cyclo_nonlcon(x, p);              % fmincon / ga constraints
% penalised objective for solvers that do NOT accept nonlcon (PSO, SA):
objP = @(x) cyclo_scalar(x, p, w) + Pen*sum( max(0, con(x)).^2 );

results = struct('name',{},'f',{},'eta',{},'V',{},'maxc',{},'t',{},'x',{});

%% 1) fmincon (SQP) -- gradient based, start-sensitive => multi-start
Nstart = 10;  bestf = inf;  bestx = [];
optF = optimoptions('fmincon','Algorithm','sqp','Display','off');
tt = tic;
for k = 1:Nstart
    x0 = lb + rand(1,nvars).*(ub-lb);
    [xk,fk,flag] = fmincon(obj,x0,[],[],[],[],lb,ub,con,optF);
    if flag > 0 && fk < bestf, bestf = fk; bestx = xk; end
end
results(end+1) = pack('fmincon (SQP, 10 starts)', bestx, p, w, toc(tt));

%% 2) ga  (handles nonlcon directly)
optG = optimoptions('ga','Display','off','PopulationSize',80, ...
                    'PlotFcn',@gaplotbestf);
tt = tic;
xg = ga(obj,nvars,[],[],[],[],lb,ub,con,optG);
results(end+1) = pack('ga', xg, p, w, toc(tt));

%% 3) particleswarm  (no nonlcon -> penalty objective)
optP = optimoptions('particleswarm','Display','off','SwarmSize',60, ...
                    'PlotFcn',@pswplotbestf);
tt = tic;
xp = particleswarm(objP,nvars,lb,ub,optP);
results(end+1) = pack('particleswarm', xp, p, w, toc(tt));

%% 4) simulannealbnd  (no nonlcon -> penalty objective)
optS = optimoptions('simulannealbnd','Display','off','PlotFcn',@saplotbestf);
x0  = lb + rand(1,nvars).*(ub-lb);
tt  = tic;
xs  = simulannealbnd(objP,x0,lb,ub,optS);
results(end+1) = pack('simulannealbnd', xs, p, w, toc(tt));

%% comparison table
fprintf('\n%-28s %9s %9s %12s %12s %8s\n', ...
        'Method','f','eta','V [mm^3]','max C','t [s]');
fprintf('%s\n', repmat('-',1,82));
for r = results
    fprintf('%-28s %9.4f %9.4f %12.0f %12.2e %8.2f\n', ...
            r.name, r.f, r.eta, r.V, r.maxc, r.t);
end
fprintf('\n(feasible if max C <= ~1e-6)\n');

%% ---- local helper: evaluate and package one result ----
function r = pack(name, x, p, w, t)
    [c,~]  = cyclo_nonlcon(x, p);
    r.name = name;
    r.f    = cyclo_scalar(x, p, w);
    r.eta  = cyclo_efficiency(x, p);
    r.V    = cyclo_volume(x, p);
    r.maxc = max(c);
    r.t    = t;
    r.x    = x;
end
