clear; clc; close all;

%% Cycloid reducer optimization
% x = [Dp, drp, B, D, K1, Dw, dsw]
% Dp  : diameter of pin wheel central circle
% drp : diameter of pin wheel
% B   : width of cycloid gear
% D   : diameter of cycloid gear center hole
% K1  : short width coefficient
% Dw  : diameter of pin central circle
% dsw : diameter of pin

nvars = 7;

lb = [140, 7.0, 7.0, 50, 0.65, 88, 11];
ub = [155, 10.4, 12, 55, 0.90, 104, 14];

x0 = [144, 10, 11, 53.5, 0.70, 90, 12];

rho = 1e8;

results = struct();

%% 1. SQP
opts_sqp = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display','iter', ...
    'MaxIterations',10000, ...
    'MaxFunctionEvaluations',50000, ...
    'OptimalityTolerance',1e-10, ...
    'ConstraintTolerance',1e-10, ...
    'StepTolerance',1e-10);

tic;
[x_sqp, f_sqp, exit_sqp] = fmincon(@single_obj, x0, [], [], [], [], lb, ub, @nonlcon, opts_sqp);
time_sqp = toc;

results.SQP.x = x_sqp;
results.SQP.fval = f_sqp;
results.SQP.volume = volume_fun(x_sqp);
results.SQP.efficiency = efficiency_fun(x_sqp);
results.SQP.max_constraint = max_constraints(x_sqp);
results.SQP.time = time_sqp;
results.SQP.exitflag = exit_sqp;

%% 2. GA
opts_ga = optimoptions('ga', ...
    'Display','iter', ...
    'PopulationSize',100, ...
    'MaxGenerations',500, ...
    'MaxStallGenerations',100, ...
    'FunctionTolerance',1e-10, ...
    'ConstraintTolerance',1e-10);

tic;
[x_ga, f_ga, exit_ga] = ga(@single_obj, nvars, [], [], [], [], lb, ub, @nonlcon, opts_ga);
time_ga = toc;

results.GA.x = x_ga;
results.GA.fval = f_ga;
results.GA.volume = volume_fun(x_ga);
results.GA.efficiency = efficiency_fun(x_ga);
results.GA.max_constraint = max_constraints(x_ga);
results.GA.time = time_ga;
results.GA.exitflag = exit_ga;

%% 3. SA with penalty
penalty_obj = @(x) single_obj(x) + rho*sum(max(0,constraint_values(x)).^2);

opts_sa = optimoptions('simulannealbnd', ...
    'Display','iter', ...
    'MaxIterations',5000, ...
    'FunctionTolerance',1e-10);

tic;
[x_sa, f_sa_penalty, exit_sa] = simulannealbnd(penalty_obj, x0, lb, ub, opts_sa);
time_sa = toc;

results.SA.x = x_sa;
results.SA.fval = single_obj(x_sa);
results.SA.penalty_fval = f_sa_penalty;
results.SA.volume = volume_fun(x_sa);
results.SA.efficiency = efficiency_fun(x_sa);
results.SA.max_constraint = max_constraints(x_sa);
results.SA.time = time_sa;
results.SA.exitflag = exit_sa;

%% 4. PSO with penalty
opts_pso = optimoptions('particleswarm', ...
    'Display','iter', ...
    'SwarmSize',100, ...
    'MaxIterations',500, ...
    'MaxStallIterations',100, ...
    'FunctionTolerance',1e-10);

tic;
[x_pso, f_pso_penalty, exit_pso] = particleswarm(penalty_obj, nvars, lb, ub, opts_pso);
time_pso = toc;

results.PSO.x = x_pso;
results.PSO.fval = single_obj(x_pso);
results.PSO.penalty_fval = f_pso_penalty;
results.PSO.volume = volume_fun(x_pso);
results.PSO.efficiency = efficiency_fun(x_pso);
results.PSO.max_constraint = max_constraints(x_pso);
results.PSO.time = time_pso;
results.PSO.exitflag = exit_pso;

%% 5. NSGA-II using gamultiobj
opts_nsga2 = optimoptions('gamultiobj', ...
    'Display','iter', ...
    'PopulationSize',150, ...
    'MaxGenerations',300, ...
    'FunctionTolerance',1e-10, ...
    'ConstraintTolerance',1e-10, ...
    'PlotFcn',{@gaplotpareto});

tic;
[x_nsga2, f_nsga2, exit_nsga2] = gamultiobj(@multi_obj, nvars, [], [], [], [], lb, ub, @nonlcon, opts_nsga2);
time_nsga2 = toc;

% Choose one representative solution from Pareto front
% Objective 1 = normalized volume
% Objective 2 = normalized efficiency loss
[~, best_idx] = min(sum(f_nsga2,2));

x_best_nsga2 = x_nsga2(best_idx,:);

results.NSGA2.x = x_best_nsga2;
results.NSGA2.fval = sum(f_nsga2(best_idx,:));
results.NSGA2.volume = volume_fun(x_best_nsga2);
results.NSGA2.efficiency = efficiency_fun(x_best_nsga2);
results.NSGA2.max_constraint = max_constraints(x_best_nsga2);
results.NSGA2.time = time_nsga2;
results.NSGA2.exitflag = exit_nsga2;
results.NSGA2.pareto_x = x_nsga2;
results.NSGA2.pareto_f = f_nsga2;

%% Summary table
algorithms = {'SQP'; 'GA'; 'SA'; 'PSO'; 'NSGA-II'};

X = [
    results.SQP.x;
    results.GA.x;
    results.SA.x;
    results.PSO.x;
    results.NSGA2.x
];

Fval = [
    results.SQP.fval;
    results.GA.fval;
    results.SA.fval;
    results.PSO.fval;
    results.NSGA2.fval
];

Volume = [
    results.SQP.volume;
    results.GA.volume;
    results.SA.volume;
    results.PSO.volume;
    results.NSGA2.volume
];

Efficiency = [
    results.SQP.efficiency;
    results.GA.efficiency;
    results.SA.efficiency;
    results.PSO.efficiency;
    results.NSGA2.efficiency
];

MaxConstraint = [
    results.SQP.max_constraint;
    results.GA.max_constraint;
    results.SA.max_constraint;
    results.PSO.max_constraint;
    results.NSGA2.max_constraint
];

Time = [
    results.SQP.time;
    results.GA.time;
    results.SA.time;
    results.PSO.time;
    results.NSGA2.time
];

Exitflag = [
    results.SQP.exitflag;
    results.GA.exitflag;
    results.SA.exitflag;
    results.PSO.exitflag;
    results.NSGA2.exitflag
];

T = table(algorithms, Fval, Volume, Efficiency, MaxConstraint, Time, Exitflag);
disp(T);

disp('Design variables [Dp drp B D K1 Dw dsw]');
disp(array2table(X, ...
    'VariableNames', {'Dp','drp','B','D','K1','Dw','dsw'}, ...
    'RowNames', algorithms));

%% Pareto front plot
figure;
scatter(results.NSGA2.pareto_f(:,1), results.NSGA2.pareto_f(:,2), 35, 'filled');
xlabel('Normalized volume objective');
ylabel('Normalized efficiency loss objective');
title('NSGA-II Pareto Front');
grid on;

figure;
scatter(Volume, Efficiency, 60, 'filled');
text(Volume, Efficiency, algorithms, 'VerticalAlignment','bottom', 'HorizontalAlignment','left');
xlabel('Volume');
ylabel('Efficiency');
title('Algorithm Comparison');
grid on;

%% Single scalar objective for SQP, GA, SA, PSO
function f = single_obj(x)
    V = volume_fun(x);
    eta = efficiency_fun(x);

    x_ref = [144, 10, 11, 53.5, 0.70, 90, 12];
    V0 = volume_fun(x_ref);
    eta0 = efficiency_fun(x_ref);

    wV = 0.5;
    wEta = 0.5;

    f = wV*(V/V0) + wEta*((1 - eta)/(1 - eta0));
end

%% Multi-objective function for NSGA-II
function f = multi_obj(x)
    V = volume_fun(x);
    eta = efficiency_fun(x);

    x_ref = [144, 10, 11, 53.5, 0.70, 90, 12];
    V0 = volume_fun(x_ref);
    eta0 = efficiency_fun(x_ref);

    f1 = V/V0;
    f2 = (1 - eta)/(1 - eta0);

    f = [f1, f2];
end

%% Volume function
function V = volume_fun(x)
    Dp  = x(1);
    drp = x(2);
    B   = x(3);
    D   = x(4);
    K1  = x(5);
    dsw = x(7);

    zc = 43;
    zp = 44;
    zw = 8;
    Delta2 = 1.0;

    term1 = Dp - K1*Dp/zp - drp;
    term2 = dsw + 2*Delta2 + K1*Dp/zp;

    V = 0.25*pi*B*(term1^2 - term2^2*zw - D^2) ...
        + K1*Dp/zp*zc*B;
end

%% Efficiency function
function eta = efficiency_fun(x)
    Dp  = x(1);
    drp = x(2);
    K1  = x(5);
    Dw  = x(6);
    dsw = x(7);

    zc = 43;

    mu = 0.07;
    fw = 0.03;
    eta_zx = 0.99;
    eta_gx = 0.995;
    Delta2 = 1.0;

    eta_x = (1 - ((Dp - drp)*4*mu)/(K1*zc*Dp*pi)) / ...
            (1 + ((Dp - drp)*4*mu)/(K1*Dp*pi));

    eta_sx = 1 - (4*fw*K1*dsw*Dp)/(pi*Dw*(dsw + 2*Delta2));

    eta = eta_x * eta_zx * eta_gx^2 * eta_sx;
end

%% Nonlinear constraints
function [c, ceq] = nonlcon(x)
    c = constraint_values(x);
    ceq = [];
end

function c = constraint_values(x)
    Dp  = x(1);
    drp = x(2);
    B   = x(3);
    D   = x(4);
    K1  = x(5);
    Dw  = x(6);
    dsw = x(7);

    zp = 44;
    zw = 8;
    Delta2 = 1.0;

    V = volume_fun(x);
    eta = efficiency_fun(x);

    pin_coeff = (Dp/drp)*sin(pi/zp);

    c = zeros(9,1);

    c(1) = 0.65 - K1;
    c(2) = K1 - 0.90;

    c(3) = 0.80 - pin_coeff;
    c(4) = pin_coeff - 1.50;

    c(5) = -V;

    c(6) = eta - 1;
    c(7) = 0.70 - eta;

    c(8) = D + 2*(dsw + 2*Delta2) - Dp;

    c(9) = dsw + 2*Delta2 - Dw/zw;
end

function m = max_constraints(x)
    m = max(constraint_values(x));
end