% WangRun_weighted.m
% =========================================================
% K-H-V형 XW3 사이클로이드 감속기
% 가중합법(Weighted Sum Method) 다목적 최적설계
%
% 목적함수: f = wV*(V/V0) + wEta*((1-eta)/(1-eta0))
% 알고리즘: fmincon SQP
% 제약조건: WangCon_fmincon.m (16개, 논문 기준)
%
% 실행 내용:
%   1. 단일 가중치(wV=0.5) SQP 최적화 + KKT 분석
%   2. 가중치 스윕(wV: 0.05~0.95) → 트레이드오프 곡선
%   3. 단목적 결과(체적만)와 비교 표
% =========================================================

clc; clear;

%% ── 0. 초기값 및 변수 범위 ────────────────────────────────
% x = [Dp, drp, B, D, K1, Dw, dsw]
x0 = [144,  10,   11,  53.5, 0.6069, 90,  12 ]';
lb = [140,   7,    7,  50,   0.60,   88,  11 ]';
ub = [155,  10.4, 12,  55,   0.90,  104,  14 ]';

%% ── 1. 초기 설계 성능 계산 ───────────────────────────────
V0   = WangObj_fmincon(x0);
eta0 = WangEff(x0);

fprintf('=== 초기 설계 ===\n');
fprintf('체적 V0   = %.0f mm³\n', V0);
fprintf('효율 eta0 = %.4f\n\n', eta0);

%% ── 2. fmincon 공통 옵션 ─────────────────────────────────
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'off', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-6);

%% ── 3. wV = 0.5 단일 실행 + KKT 분석 ────────────────────
wV   = 0.5;
wEta = 1 - wV;

obj_w = @(x) WangObj_weighted(x, wV, wEta, V0, eta0);

options_verbose = optimoptions(options, 'Display', 'iter');

fprintf('=== 가중합 최적화 (wV=%.1f, wEta=%.1f) ===\n', wV, wEta);
[xopt, fopt, exitflag, output, lambda] = fmincon( ...
    obj_w, x0, [], [], [], [], lb, ub, @WangCon_fmincon, options_verbose);

V_opt   = WangObj_fmincon(xopt);
eta_opt = WangEff(xopt);

fprintf('\n종료 조건 (exitflag): %d', exitflag);
switch exitflag
    case  1; fprintf(' → KKT 조건 만족 (정상 수렴)\n');
    case  2; fprintf(' → 스텝 크기 허용오차 수렴\n');
    case  0; fprintf(' → 최대 반복 초과\n');
    case -2; fprintf(' → 가용 영역 없음\n');
    otherwise; fprintf('\n');
end
fprintf('반복 횟수: %d, 함수 평가: %d\n\n', output.iterations, output.funcCount);

%% 설계변수 비교
names = {'Dp (핀중심원직경)','drp (핀직경)','B (기어폭)', ...
         'D (중심홀직경)','K1 (단폭계수)','Dw (출력핀중심원)','dsw (출력핀직경)'};
units = {'mm','mm','mm','mm','-','mm','mm'};

fprintf('=== 설계변수 비교 (wV=0.5) ===\n');
fprintf('  %-22s  %8s  %8s\n','변수','초기값','최적값');
fprintf('  %s\n', repmat('-',1,42));
for i = 1:7
    fprintf('  %-22s  %8.4f  %8.4f  %s\n', names{i}, x0(i), xopt(i), units{i});
end

fprintf('\n=== 목적함수 및 성능 비교 ===\n');
fprintf('  %-15s  %10s  %10s  %8s\n','항목','초기','최적','변화');
fprintf('  %s\n', repmat('-',1,48));
fprintf('  %-15s  %10.0f  %10.0f  %+7.1f%%\n','체적 (mm³)', V0,   V_opt,   100*(V_opt-V0)/V0);
fprintf('  %-15s  %10.4f  %10.4f  %+7.4f\n', '효율 η',     eta0, eta_opt, eta_opt-eta0);

%% KKT 라그랑지 승수 분석
fprintf('\n=== KKT 라그랑지 승수 분석 (wV=0.5) ===\n');
constraint_names = { ...
    'y1  언더컷/첨예 방지', ...
    'y2  단폭계수 하한 K1>=0.60', ...
    'y3  단폭계수 상한 K1<=0.9', ...
    'y4  핀직경계수 K2 하한', ...
    'y5  핀직경계수 K2 상한', ...
    'y6  사이클로이드 접촉강도', ...
    'y7  핀기어 굽힘강도', ...
    'y8  핀-핀홀 접촉강도', ...
    'y9  핀 굽힘강도', ...
    'y10 Dp 하한 >=140mm', ...
    'y11 Dp 상한 <=155mm', ...
    'y12 핀홀 직경 조건 1', ...
    'y13 핀홀 직경 조건 2', ...
    'y14 기어폭 하한 B>=0.05Dp', ...
    'y15 기어폭 상한 B<=0.1Dp', ...
    'y16 피벗 베어링 수명 >=5000h'};

lam = lambda.ineqnonlin;
fprintf('활성 제약 (λ > 1e-4):\n');
active_any = false;
for i = 1:length(lam)
    if abs(lam(i)) > 1e-4
        fprintf('  [활성] c(%2d) %s\n         λ = %.6f\n', i, constraint_names{i}, lam(i));
        active_any = true;
    end
end
if ~active_any
    fprintf('  활성 제약 없음\n');
end

fprintf('\n경계 활성 변수 (λ > 1e-4):\n');
lam_lb = lambda.lower;
lam_ub = lambda.upper;
bound_any = false;
for i = 1:7
    if abs(lam_lb(i)) > 1e-4
        fprintf('  [하한] %s  λ=%.6f\n', names{i}, lam_lb(i));
        bound_any = true;
    end
    if abs(lam_ub(i)) > 1e-4
        fprintf('  [상한] %s  λ=%.6f\n', names{i}, lam_ub(i));
        bound_any = true;
    end
end
if ~bound_any
    fprintf('  경계 활성 변수 없음\n');
end

%% 제약 만족 확인
fprintf('\n=== 최적해 제약 위반 확인 ===\n');
[c_opt, ~] = WangCon_fmincon(xopt);
any_viol = false;
for i = 1:16
    if c_opt(i) > options.ConstraintTolerance
        fprintf('  c(%2d) = %.6f  위반\n', i, c_opt(i));
        any_viol = true;
    end
end
if ~any_viol
    fprintf('  모든 제약 만족\n');
end

%% ── 4. 가중치 스윕 → 트레이드오프 곡선 ─────────────────
fprintf('\n=== 가중치 스윕 실행 중 ===\n');

wV_vec = 0.05 : 0.05 : 0.95;   % wV = 0.05, 0.10, ..., 0.95 (19개)
nW = length(wV_vec);

V_sweep   = zeros(nW, 1);
eta_sweep = zeros(nW, 1);
exit_sweep = zeros(nW, 1);
x_sweep   = zeros(nW, 7);

for k = 1:nW
    wVk   = wV_vec(k);
    wEtak = 1 - wVk;

    obj_k = @(x) WangObj_weighted(x, wVk, wEtak, V0, eta0);

    [xk, ~, exk] = fmincon(obj_k, x0, [], [], [], [], lb, ub, ...
                            @WangCon_fmincon, options);

    V_sweep(k)    = WangObj_fmincon(xk);
    eta_sweep(k)  = WangEff(xk);
    exit_sweep(k) = exk;
    x_sweep(k,:)  = xk;

    fprintf('  wV = %.2f  →  V = %7.0f mm³,  η = %.4f  (exit=%d)\n', ...
        wVk, V_sweep(k), eta_sweep(k), exk);
end

%% ── 5. 단목적 결과 참고값 (체적만 최소화) ───────────────
% 단목적 체적 최소화 재실행 (비교용)
obj_V = @(x) WangObj_fmincon(x) / V0;
[x_single, ~, ~] = fmincon(obj_V, x0, [], [], [], [], lb, ub, ...
                            @WangCon_fmincon, options);
V_single   = WangObj_fmincon(x_single);
eta_single = WangEff(x_single);

%% ── 6. 결과 비교 표 ─────────────────────────────────────
fprintf('\n=== 결과 비교 표 ===\n');
fprintf('  %-20s  %10s  %8s\n','방법','체적 (mm³)','효율 η');
fprintf('  %s\n', repmat('-',1,42));
fprintf('  %-20s  %10.0f  %8.4f\n', '초기 설계',             V0,        eta0);
fprintf('  %-20s  %10.0f  %8.4f\n', '단목적 (체적만)',        V_single,  eta_single);
fprintf('  %-20s  %10.0f  %8.4f\n', '가중합 wV=0.5',          V_opt,     eta_opt);

%% ── 7. 트레이드오프 곡선 플롯 ───────────────────────────
figure('Name','체적-효율 트레이드오프 곡선','NumberTitle','off');
hold on;

% 스윕 결과
scatter(V_sweep, eta_sweep, 50, wV_vec, 'filled');
colorbar; colormap(jet);
clim([0 1]);

% 가중치 레이블 (일부만 표시)
label_idx = [1, 5, 10, 15, 19];  % wV = 0.05, 0.25, 0.50, 0.75, 0.95
for k = label_idx
    text(V_sweep(k)+200, eta_sweep(k), sprintf('w_V=%.2f', wV_vec(k)), ...
        'FontSize', 8, 'Color', 'k');
end

% 단목적 및 초기점 표시
scatter(V0,       eta0,       100, 'k', 's', 'filled', 'DisplayName','초기 설계');
scatter(V_single, eta_single, 100, 'r', '^', 'filled', 'DisplayName','단목적(체적)');
scatter(V_opt,    eta_opt,    100, 'b', 'p', 'filled', 'DisplayName','가중합 w_V=0.5');

legend('Location','best');
xlabel('체적 V (mm³)');
ylabel('전달 효율 η');
title('체적-효율 트레이드오프 곡선 (가중치 스윕, SQP)');
grid on;
hold off;

%% ── 8. 설계변수 변화 플롯 (가중치에 따라) ───────────────
var_labels = {'D_p','d_{rp}','B','D','K_1','D_w','d_{sw}'};

figure('Name','가중치별 최적 설계변수','NumberTitle','off');
for i = 1:7
    subplot(2,4,i);
    plot(wV_vec, x_sweep(:,i), 'b-o', 'MarkerSize', 4);
    yline(lb(i), 'r--', 'LB');
    yline(ub(i), 'g--', 'UB');
    xlabel('w_V (체적 가중치)');
    ylabel(var_labels{i});
    title(var_labels{i});
    grid on;
end
sgtitle('가중치 w_V에 따른 최적 설계변수 변화');

%% ── Figure 4: 단목적 vs 가중합 설계변수 비교 ────────────
names_kor = {'D_p','d_{rp}','B','D','K_1','D_w','d_{sw}'};
lb_w = [140, 7, 7, 50, 0.60, 88, 11]';
ub_w = [155, 10.4, 12, 55, 0.90, 104, 14]';

% 단목적 최적값 재계산 (비교용)
options_fig = optimoptions('fmincon','Algorithm','sqp','Display','off', ...
    'OptimalityTolerance',1e-6,'ConstraintTolerance',1e-6);
obj_V_only = @(x) WangObj_fmincon(x) / V0;
[x_single_fig, ~] = fmincon(obj_V_only, x0, [], [], [], [], lb_w, ub_w, ...
    @WangCon_fmincon, options_fig);

figure('Name','Figure 4: 단목적 vs 가중합 설계변수 비교','NumberTitle','off');

x_ax = 1:7;
bw   = 0.25;

hold on;
b0 = bar(x_ax - bw,   x0',           bw, 'FaceColor',[0.75 0.75 0.75], 'EdgeColor','w');
b1 = bar(x_ax,        x_single_fig', bw, 'FaceColor',[0.40 0.65 0.95], 'EdgeColor','w');
b2 = bar(x_ax + bw,   xopt',         bw, 'FaceColor',[1.00 0.50 0.30], 'EdgeColor','w');

% 상·하한선
for i = 1:7
    plot([i-bw*1.8, i+bw*1.8], [lb_w(i), lb_w(i)], 'r--', 'LineWidth',1);
    plot([i-bw*1.8, i+bw*1.8], [ub_w(i), ub_w(i)], 'g--', 'LineWidth',1);
end
hold off;

set(gca, 'XTick', x_ax, 'XTickLabel', names_kor, 'FontSize',10);
ylabel('설계변수 값 (mm 또는 무차원)');
legend([b0, b1, b2], {'초기값','단목적 최적','가중합 최적 (wV=0.5)'}, 'Location','best');
title('Figure 4: 단목적 vs 가중합(wV=0.5) — 설계변수 비교');
grid on;
