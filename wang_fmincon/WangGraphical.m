% WangGraphical.m
% =========================================================
% XW3형 사이클로이드 감속기 체적 최소화 — 도식해법
%
% 설계변수 2개 (나머지 4개는 단목적 최적값으로 고정):
%   x1 = drp  (핀 직경, mm)  : 7.0 ~ 10.4
%   x2 = K1   (단폭계수)     : 0.60 ~ 0.90
%
% 고정 설계변수: Dp=140, B=7, D=55, dsw=14, Dw=90
% 수업 연결: Lec 05 도식해법 — meshgrid + contour
% =========================================================

clc; clear; close all;

%% ── 0. 고정 파라미터 ──────────────────────────────────
Dp = 140; B = 7; D = 55; dsw = 14; Dw = 90;
zc = 43; zp = 44; zw = 8; Delta2 = 2;

%% ── 1. 격자 생성 ──────────────────────────────────────
drp_vec = linspace(7.0, 10.4, 600);
K1_vec  = linspace(0.60, 0.90, 600);
[drp, K1] = meshgrid(drp_vec, K1_vec);

%% ── 2. 목적함수 — 체적 V ─────────────────────────────
V = (pi/4) .* B .* ( ...
    (Dp - K1.*Dp/zp - drp).^2 ...
    - (dsw + 2*Delta2 + K1.*Dp/zp).^2 .* zw ...
    - D^2 ) ...
    + K1 .* (Dp/zp) .* zc .* B;

%% ── 3. 제약조건 함수 ─────────────────────────────────
% y1: 언더컷/첨예화 방지
g1 = drp./2 - Dp/2 .* sqrt(27.*(1-K1.^2).*(zp-1)./(zp+1)^3);

% y4: K2 하한 (K2 >= 1.0) → g4 = 1 - K2 <= 0
K2 = (Dp ./ drp) .* sin(pi/zp);
g4 = 1.0 - K2;

% y5: K2 상한 (K2 <= 1.6) → drp 최솟값 기준 6.24mm로 범위 밖이므로 생략

%% ── 4. Infeasible 영역 ───────────────────────────────
infeasible = (g1 > 0) | (g4 > 0);
feasible   = ~infeasible;

%% ── 5. 최적해 및 초기점 ──────────────────────────────
drp_opt = 9.990; K1_opt = 0.775;
V_opt = (pi/4)*B*((Dp-K1_opt*Dp/zp-drp_opt)^2 ...
    -(dsw+2*Delta2+K1_opt*Dp/zp)^2*zw-D^2) ...
    +K1_opt*(Dp/zp)*zc*B;

drp0 = 10.0; K10 = 0.6069;
V0 = (pi/4)*B*((Dp-K10*Dp/zp-drp0)^2 ...
    -(dsw+2*Delta2+K10*Dp/zp)^2*zw-D^2) ...
    +K10*(Dp/zp)*zc*B;

fprintf('최적해: drp* = %.3f mm,  K1* = %.4f,  V* = %.0f mm³\n', drp_opt, K1_opt, V_opt);
fprintf('초기점: drp0 = %.3f mm,  K10 = %.4f,  V0 = %.0f mm³\n', drp0, K10, V0);

%% ── 6. 도식해법 그래프 ───────────────────────────────
figure('Name','도식해법 — drp vs K1 체적 최소화', ...
       'NumberTitle','off','Position',[80 80 950 720]);
hold on;

%% (a) Infeasible 영역 음영
contourf(drp, K1, double(infeasible), [0.5 0.5], ...
    'FaceColor',[1 0.87 0.87], 'EdgeColor','none');

%% (b) 체적 등축선 — 전체 영역에 그려야 최적해가 경계 교점임을 확인 가능
V_levels = [52000, 53000, 54000, 55128, 56000, 57000, 59000, 62000, 66000];
[C_V, h_V] = contour(drp, K1, V, V_levels, ...
    'k--', 'LineWidth', 1.0);
clabel(C_V, h_V, [52000, 55128, 57000, 62000], ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3], 'LabelSpacing', 300);

%% (c) 제약조건 경계선
[~, hg1] = contour(drp, K1, g1, [0 0], 'k-',  'LineWidth', 2.5);
[~, hg4] = contour(drp, K1, g4, [0 0], 'b-',  'LineWidth', 2.5);

%% (d) 최적해
scatter(drp_opt, K1_opt, 160, 'r', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.5, ...
    'DisplayName', sprintf('최적해  V^* = %.0f mm³', V_opt));
text(drp_opt + 0.05, K1_opt + 0.008, ...
    sprintf('X^* = (%.3f,  %.3f)\nV^* = %.0f mm³', drp_opt, K1_opt, V_opt), ...
    'FontSize', 9, 'Color', 'r', 'FontWeight', 'bold', ...
    'BackgroundColor','w', 'EdgeColor','r', 'Margin', 3);

%% (e) 초기점
scatter(drp0, K10, 130, 'bs', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.5, ...
    'DisplayName', sprintf('초기점  V_0 = %.0f mm³', V0));
text(drp0 + 0.05, K10 + 0.008, ...
    sprintf('X_0 = (%.1f,  %.4f)\nV_0 = %.0f mm³', drp0, K10, V0), ...
    'FontSize', 9, 'Color', 'b', ...
    'BackgroundColor','w', 'EdgeColor','b', 'Margin', 3);

%% (f) 가용 영역 레이블
text(8.0, 0.835, 'Feasible Region', ...
    'FontSize', 13, 'FontWeight', 'bold', ...
    'Color', [0 0.5 0], 'HorizontalAlignment','center');

%% (g) Infeasible 레이블 (각 영역 중심에 배치)
% g1 위반 (우상단 빨간 영역)
text(9.0, 0.88, {'Infeasible'; '(y1: 언더컷 위반)'}, ...
    'FontSize', 8, 'Color', [0.7 0 0], 'HorizontalAlignment','center');
% g4 위반 (우측 — drp 너무 굵음)
text(10.25, 0.72, {'Infeasible'; '(y4: K_2 < 1.0)'}, ...
    'FontSize', 8, 'Color', [0 0 0.8], 'HorizontalAlignment','center', ...
    'Rotation', 90);

%% (h) 제약 경계선 레이블
% g1=0 경계 (곡선 위 적당한 위치)
text(8.3, 0.865, 'g_1 = 0', 'FontSize', 9, 'Color', 'k', ...
    'FontWeight','bold', 'BackgroundColor','w');
% g4=0 경계 (수직선 근처)
text(9.82, 0.83, 'g_4 = 0\n(K_2 = 1.0)', 'FontSize', 9, 'Color', 'b', ...
    'FontWeight','bold', 'BackgroundColor','w', 'HorizontalAlignment','right');
%% (i) 체적 감소 방향 화살표 (annotation 사용)
% 데이터 좌표 → normalized figure 좌표로 변환
ax = gca;
x_range = [7.0, 10.4]; y_range = [0.60, 0.90];
ax.XLim = x_range; ax.YLim = y_range;

% 화살표 시작점/끝점 (데이터 좌표)
ax1 = 8.3; ay1 = 0.770;   % 시작 (좌하단, 큰 체적)
ax2 = 9.0; ay2 = 0.785;   % 끝 (우상단, 작은 체적, 최적해 방향)

% 그래프 영역의 normalized 좌표로 변환
pos = ax.Position;
norm_x1 = pos(1) + (ax1 - x_range(1))/(x_range(2)-x_range(1)) * pos(3);
norm_y1 = pos(2) + (ay1 - y_range(1))/(y_range(2)-y_range(1)) * pos(4);
norm_x2 = pos(1) + (ax2 - x_range(1))/(x_range(2)-x_range(1)) * pos(3);
norm_y2 = pos(2) + (ay2 - y_range(1))/(y_range(2)-y_range(1)) * pos(4);

annotation('arrow', [norm_x1 norm_x2], [norm_y1 norm_y2], ...
    'Color', [0 0.6 0], 'LineWidth', 2.0, 'HeadWidth', 10, 'HeadLength', 10);
text(ax2 + 0.05, ay2 + 0.010, 'V 감소 방향', ...
    'FontSize', 9, 'Color', [0 0.6 0], 'FontWeight','bold');

%% 그래프 마무리
xlabel('d_{rp}  (핀 직경, mm)', 'FontSize', 12, 'FontWeight','bold');
ylabel('K_1  (단폭계수)', 'FontSize', 12, 'FontWeight','bold');
title({'XW3형 사이클로이드 감속기 체적 최소화 — 도식해법'; ...
    sprintf('(고정: D_p=%.0f, B=%.0f, D=%.0f, d_{sw}=%.0f mm)', Dp,B,D,dsw)}, ...
    'FontSize', 11);

legend({'Infeasible 영역', '체적 등축선 V (mm³)', ...
        'g_1=0  언더컷/첨예화 방지', ...
        'g_4=0  K_2 하한 (K_2 \geq 1.0)', ...
        sprintf('최적해  V^* = %.0f mm³', V_opt), ...
        sprintf('초기점  V_0 = %.0f mm³', V0)}, ...
    'Location','northwest', 'FontSize', 8.5);

xlim(x_range); ylim(y_range);
grid on; box on;
set(gca, 'FontSize', 10);
hold off;
