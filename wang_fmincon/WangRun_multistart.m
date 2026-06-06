% WangRun_multistart.m
% =========================================================
% 다중 초기점(Multi-Start) SQP 실험
%
% 단목적(체적만)과 가중합(wV=0.5) 두 경우를 동일한
% 랜덤 초기점 집합으로 실험하여 결과를 비교한다.
%
% 필요 파일 (같은 폴더):
%   WangObj_fmincon.m, WangCon_fmincon.m,
%   WangObj_weighted.m, WangEff.m
% =========================================================

clc; clear;
rng(42);   % 재현성을 위한 난수 시드 고정

%% ── 0. 설계변수 범위 및 기준 초기점 ─────────────────────
lb = [140,   7,    7,  50,  0.60,  88,  11]';
ub = [155,  10.4, 12,  55,  0.90, 104,  14]';

x0_ref = [144, 10, 11, 53.5, 0.6069, 90, 12]';  % 논문 초기값

%% ── 1. 가중합 파라미터 및 기준값 ────────────────────────
wV   = 0.5;
wEta = 0.5;
V0   = WangObj_fmincon(x0_ref);
eta0 = WangEff(x0_ref);

obj_single = @(x) WangObj_fmincon(x);
obj_weighted = @(x) WangObj_weighted(x, wV, wEta, V0, eta0);

%% ── 2. fmincon 옵션 ──────────────────────────────────────
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'off', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-6);

%% ── 3. 기준 초기점 결과 (비교용) ────────────────────────
[x_s_ref, V_s_ref, ex_s_ref] = fmincon(obj_single, x0_ref, ...
    [], [], [], [], lb, ub, @WangCon_fmincon, options);
[x_w_ref, ~, ex_w_ref] = fmincon(obj_weighted, x0_ref, ...
    [], [], [], [], lb, ub, @WangCon_fmincon, options);
V_w_ref   = WangObj_fmincon(x_w_ref);
eta_w_ref = WangEff(x_w_ref);

fprintf('=== 기준 초기점 결과 ===\n');
fprintf('단목적:  V* = %7.0f mm³  (exit=%d)\n', V_s_ref, ex_s_ref);
fprintf('가중합:  V* = %7.0f mm³,  η* = %.4f  (exit=%d)\n\n', ...
    V_w_ref, eta_w_ref, ex_w_ref);

%% ── 4. 랜덤 초기점 생성 (두 실험 공통) ──────────────────
N = 50;

% 결과 저장: [x(7) V eta exitflag]
res_s = zeros(N, 10);   % 단목적
res_w = zeros(N, 10);   % 가중합

fprintf('=== 랜덤 초기점 %d개 실험 중 ===\n', N);
fprintf('  %-4s  %10s  %10s  %10s\n', 'k', '단목적 V', '가중합 V', '가중합 η');
fprintf('  %s\n', repmat('-', 1, 42));

for k = 1:N
    x_init = lb + (ub - lb) .* rand(7, 1);

    % 단목적
    [xk_s, Vk_s, exk_s] = fmincon(obj_single, x_init, ...
        [], [], [], [], lb, ub, @WangCon_fmincon, options);
    res_s(k,:) = [xk_s', Vk_s, WangEff(xk_s), exk_s];

    % 가중합 (동일 초기점)
    [xk_w, ~, exk_w] = fmincon(obj_weighted, x_init, ...
        [], [], [], [], lb, ub, @WangCon_fmincon, options);
    Vk_w   = WangObj_fmincon(xk_w);
    etak_w = WangEff(xk_w);
    res_w(k,:) = [xk_w', Vk_w, etak_w, exk_w];

    fprintf('  [%2d]  %10.0f  %10.0f  %10.4f\n', k, Vk_s, Vk_w, etak_w);
end

%% ── 5. 수렴 해 분류 ──────────────────────────────────────
conv_s = res_s(res_s(:,9) > 0, :);
conv_w = res_w(res_w(:,9) > 0, :);

conv_s = sortrows(conv_s, 8);
conv_w = sortrows(conv_w, 8);

fprintf('\n=== 수렴 통계 ===\n');
fprintf('단목적 — 수렴: %d/%d, V 범위: %.0f ~ %.0f mm³\n', ...
    size(conv_s,1), N, conv_s(1,8), conv_s(end,8));
fprintf('가중합 — 수렴: %d/%d, V 범위: %.0f ~ %.0f mm³, η 범위: %.4f ~ %.4f\n', ...
    size(conv_w,1), N, conv_w(1,8), conv_w(end,8), ...
    min(conv_w(:,9)), max(conv_w(:,9)));

%% ── 6. 클러스터 분석 ─────────────────────────────────────
tol = 1000;   % 허용오차 mm³

function clusters = find_clusters(V_vals, tol)
    clusters = {};
    used = false(length(V_vals), 1);
    for i = 1:length(V_vals)
        if used(i), continue; end
        idx = abs(V_vals - V_vals(i)) < tol;
        clusters{end+1} = struct('V_mean', mean(V_vals(idx)), 'count', sum(idx));
        used(idx) = true;
    end
end

cl_s = find_clusters(conv_s(:,8), tol);
cl_w = find_clusters(conv_w(:,8), tol);

fprintf('\n=== 국소 최적해 클러스터 (허용오차 ±%d mm³) ===\n', tol);
fprintf('단목적 (%d개 클러스터):\n', length(cl_s));
for c = 1:length(cl_s)
    fprintf('  클러스터 %d: V ≈ %7.0f mm³  (%d개)\n', c, cl_s{c}.V_mean, cl_s{c}.count);
end
fprintf('가중합 (%d개 클러스터):\n', length(cl_w));
for c = 1:length(cl_w)
    fprintf('  클러스터 %d: V ≈ %7.0f mm³  (%d개)\n', c, cl_w{c}.V_mean, cl_w{c}.count);
end

%% ── 7. 시각화 ───────────────────────────────────────────
% 
% % (1) 체적 분포 히스토그램 비교
% figure('Name', '체적 분포 비교', 'NumberTitle', 'off');
% 
% subplot(1,2,1);
% histogram(conv_s(:,8), 12, 'FaceColor', [0.3 0.6 0.9], 'EdgeColor', 'w');
% hold on;
% xline(V_s_ref, 'r-', sprintf('기준 (%.0f)', V_s_ref), 'LineWidth', 2);
% xlabel('체적 V (mm³)'); ylabel('빈도');
% title('단목적 (체적만)');
% grid on; hold off;
% 
% subplot(1,2,2);
% histogram(conv_w(:,8), 12, 'FaceColor', [1.0 0.6 0.3], 'EdgeColor', 'w');
% hold on;
% xline(V_w_ref, 'r-', sprintf('기준 (%.0f)', V_w_ref), 'LineWidth', 2);
% xlabel('체적 V (mm³)'); ylabel('빈도');
% title('가중합 (w_V=0.5)');
% grid on; hold off;
% 
% sgtitle('다중 초기점 SQP — 체적 분포 비교');
% 
% % (2) 체적-효율 산점도 (두 경우 겹쳐서)
% figure('Name', '체적-효율 분포 비교', 'NumberTitle', 'off');
% hold on;
% scatter(conv_s(:,8), conv_s(:,9), 40, 'b', 'filled', ...
%     'DisplayName', '단목적 (랜덤)');
% scatter(conv_w(:,8), conv_w(:,9), 40, [1 0.5 0], 'filled', ...
%     'DisplayName', '가중합 (랜덤)');
% scatter(V_s_ref, WangEff(x_s_ref), 120, 'b', '^', 'filled', ...
%     'DisplayName', '단목적 (기준초기점)');
% scatter(V_w_ref, eta_w_ref, 120, [1 0.5 0], '^', 'filled', ...
%     'DisplayName', '가중합 (기준초기점)');
% scatter(V0, eta0, 120, 'k', 's', 'filled', 'DisplayName', '초기 설계');
% xlabel('체적 V (mm³)');
% ylabel('전달 효율 η');
% title('다중 초기점 SQP — 체적-효율 분포');
% legend('Location', 'best');
% grid on; hold off;
% 
% % (3) K1 분포 비교
% figure('Name', 'K1 분포 비교', 'NumberTitle', 'off');
% hold on;
% scatter(conv_s(:,5), conv_s(:,8), 40, 'b', 'filled', 'DisplayName', '단목적');
% scatter(conv_w(:,5), conv_w(:,8), 40, [1 0.5 0], 'filled', 'DisplayName', '가중합');
% xlabel('K_1 (단폭계수)');
% ylabel('체적 V (mm³)');
% title('최적해 분포: K_1 vs 체적');
% legend('Location', 'best');
% grid on; hold off;
