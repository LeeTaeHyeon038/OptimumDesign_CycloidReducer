% WangGraphical.m
% =========================================================
% XW3형 사이클로이드 감속기 체적 최소화 — 도식해법
%
% 설계변수 2개 (나머지 4개는 단목적 최적값으로 고정):
%   x1 = drp  (핀 직경, mm)  : 6.0 ~ 11.0  (여유 포함)
%   x2 = K1   (단폭계수)     : 0.55 ~ 0.95 (여유 포함)
%
% 고정 설계변수: Dp=140, B=7, D=55, dsw=14, Dw=90
% 수업 연결: Lec 05 도식해법 — meshgrid + contour
%
% 표시 제약:
%   y1  : 언더컷/첨예화 방지  (비선형 곡선)
%   y4  : K2 하한 (K2 >= 1.0) (수직선)
%   y5  : K2 상한 (K2 <= 1.6) (수직선, drp=6.24mm)
%   lb  : drp 하한 (drp >= 7.0) (수직선)
%   ub  : drp 상한 (drp <= 10.4) (수직선)
%   K1 lb: K1 하한 (K1 >= 0.60) (수평선)
%   K1 ub: K1 상한 (K1 <= 0.90) (수평선)
% =========================================================

clc; clear; close all;

%% ── 0. 고정 파라미터 ──────────────────────────────────
Dp = 140; B = 7; D = 55; dsw = 14; Dw = 90;
zc = 43; zp = 44; zw = 8; Delta2 = 2;

% 설계변수 상·하한
drp_lb = 7.0; drp_ub = 10.4;
K1_lb  = 0.60; K1_ub  = 0.90;

%% ── 1. 격자 생성 (축 범위를 넓게 설정) ──────────────
drp_vec = linspace(5.5, 11.5, 700);
K1_vec  = linspace(0.50, 0.95, 700);
[drp, K1] = meshgrid(drp_vec, K1_vec);

%% ── 2. 목적함수 — 체적 V ─────────────────────────────
V = (pi/4) .* B .* ( ...
    (Dp - K1.*Dp/zp - drp).^2 ...
    - (dsw + 2*Delta2 + K1.*Dp/zp).^2 .* zw ...
    - D^2 ) ...
    + K1 .* (Dp/zp) .* zc .* B;

%% ── 3. 제약조건 함수 ─────────────────────────────────
% y1: 언더컷/첨예화 방지 (K1 >= 0.60 범위에서 항상 첫 번째 수식)
g1 = drp./2 - Dp/2 .* sqrt(27.*(1-K1.^2).*(zp-1)./(zp+1)^3);

% K2 계산
K2 = (Dp ./ drp) .* sin(pi/zp);

% y4: K2 하한 (K2 >= 1.0)
g4 = 1.0 - K2;

% y5: K2 상한 (K2 <= 1.6)  → drp = Dp*sin(pi/zp)/1.6 = 6.242mm
g5 = K2 - 1.6;

%% ── 4. 가용 영역 마스크 ─────────────────────────────
% 모든 제약 + 설계변수 경계 고려
infeasible = (g1 > 0) | (g4 > 0) | (g5 > 0) ...
           | (drp < drp_lb) | (drp > drp_ub) ...
           | (K1  < K1_lb)  | (K1  > K1_ub);

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

% K2=1.6 경계값
drp_K2_16 = Dp * sin(pi/zp) / 1.6;
fprintf('y5 경계: drp = %.4f mm  (K2=1.6)\n', drp_K2_16);

%% ── 6. 도식해법 그래프 ───────────────────────────────
figure('Name','도식해법 — drp vs K1 체적 최소화', ...
       'NumberTitle','off','Position',[80 80 1000 750]);
hold on;

%% (a) Infeasible 영역 음영
contourf(drp, K1, double(infeasible), [0.5 0.5], ...
    'FaceColor',[1 0.87 0.87], 'EdgeColor','none');

%% (b) 체적 등축선 — 전체 영역
V_levels = [50000, 52000, 54000, 55128, 57000, 60000, 65000, 72000, 82000];
[C_V, h_V] = contour(drp, K1, V, V_levels, 'k--', 'LineWidth', 1.0);
clabel(C_V, h_V, [52000, 55128, 60000, 72000], ...
    'FontSize', 8, 'Color', [0.3 0.3 0.3], 'LabelSpacing', 300);

%% (c) 비선형 제약 경계선
[~, hg1] = contour(drp, K1, g1, [0 0], 'k-',  'LineWidth', 2.5);
[~, hg4] = contour(drp, K1, g4, [0 0], 'b-',  'LineWidth', 2.5);
[~, hg5] = contour(drp, K1, g5, [0 0], 'm-',  'LineWidth', 2.0);

%% (d) 설계변수 경계선 (상·하한)
x_range = [5.5, 11.5]; y_range = [0.50, 0.95];
% drp 하한
plot([drp_lb drp_lb], y_range, 'r-',  'LineWidth', 1.8);
% drp 상한
plot([drp_ub drp_ub], y_range, 'r-',  'LineWidth', 1.8);
% K1 하한
plot(x_range, [K1_lb K1_lb], 'r--', 'LineWidth', 1.8);
% K1 상한
plot(x_range, [K1_ub K1_ub], 'r--', 'LineWidth', 1.8);

%% (e) 최적해
scatter(drp_opt, K1_opt, 160, 'r', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.5, ...
    'DisplayName', sprintf('최적해  V^* = %.0f mm³', V_opt));
text(drp_opt + 0.08, K1_opt + 0.008, ...
    sprintf('X^* = (%.3f,  %.3f)\nV^* = %.0f mm³', drp_opt, K1_opt, V_opt), ...
    'FontSize', 9, 'Color', 'r', 'FontWeight','bold', ...
    'BackgroundColor','w', 'EdgeColor','r', 'Margin', 3);

%% (f) 초기점
scatter(drp0, K10, 130, 'bs', 'filled', ...
    'MarkerEdgeColor','k', 'LineWidth',1.5, ...
    'DisplayName', sprintf('초기점  V_0 = %.0f mm³', V0));
text(drp0 + 0.08, K10 + 0.008, ...
    sprintf('X_0 = (%.1f,  %.4f)\nV_0 = %.0f mm³', drp0, K10, V0), ...
    'FontSize', 9, 'Color', 'b', ...
    'BackgroundColor','w', 'EdgeColor','b', 'Margin', 3);

%% (g) 가용 영역 레이블
text(8.5, 0.78, 'Feasible Region', ...
    'FontSize', 13, 'FontWeight','bold', ...
    'Color', [0 0.5 0], 'HorizontalAlignment','center');

%% (h) Infeasible 레이블
text(8.5, 0.92, {'Infeasible'; '(y1: 언더컷 위반)'}, ...
    'FontSize', 8, 'Color', [0.6 0 0], 'HorizontalAlignment','center');
text(10.7, 0.72, {'Infeasible'; '(y4: K_2<1.0)'}, ...
    'FontSize', 8, 'Color', [0 0 0.8], 'HorizontalAlignment','center', ...
    'Rotation', 90);
text(5.85, 0.72, {'Infeasible'; '(y5: K_2>1.6)'}, ...
    'FontSize', 8, 'Color', [0.7 0 0.7], 'HorizontalAlignment','center', ...
    'Rotation', 90);
text(6.5, 0.54, {'Infeasible'; '(drp < 7.0)'}, ...
    'FontSize', 8, 'Color', [0.8 0.2 0], 'HorizontalAlignment','center');
text(11.0, 0.54, {'Infeasible'; '(drp > 10.4)'}, ...
    'FontSize', 8, 'Color', [0.8 0.2 0], 'HorizontalAlignment','center');

%% (i) 제약 경계선 레이블
text(7.8, 0.91, 'g_1=0', 'FontSize', 9, 'Color', 'k', ...
    'FontWeight','bold', 'BackgroundColor','w');
text(10.15, 0.86, 'g_4=0 (K_2=1.0)', 'FontSize', 8, 'Color', 'b', ...
    'FontWeight','bold', 'BackgroundColor','w', 'HorizontalAlignment','center', ...
    'Rotation', 90);
text(6.42, 0.86, 'g_5=0 (K_2=1.6)', 'FontSize', 8, 'Color', 'm', ...
    'BackgroundColor','w', 'HorizontalAlignment','center', 'Rotation', 90);
text(7.05, 0.88, 'd_{rp}=7.0', 'FontSize', 8, 'Color', 'r', ...
    'BackgroundColor','w', 'Rotation', 90);
text(10.45, 0.88, 'd_{rp}=10.4', 'FontSize', 8, 'Color', 'r', ...
    'BackgroundColor','w', 'Rotation', 90);
text(6.2, 0.605, 'K_1=0.60', 'FontSize', 8, 'Color', 'r', 'BackgroundColor','w');
text(6.2, 0.905, 'K_1=0.90', 'FontSize', 8, 'Color', 'r', 'BackgroundColor','w');

%% (j) 체적 감소 방향 화살표
ax = gca;
ax.XLim = x_range; ax.YLim = y_range;
pos = ax.Position;

ax1 = 7.8; ay1 = 0.730;
ax2 = 8.7; ay2 = 0.750;

norm_x1 = pos(1) + (ax1-x_range(1))/(x_range(2)-x_range(1))*pos(3);
norm_y1 = pos(2) + (ay1-y_range(1))/(y_range(2)-y_range(1))*pos(4);
norm_x2 = pos(1) + (ax2-x_range(1))/(x_range(2)-x_range(1))*pos(3);
norm_y2 = pos(2) + (ay2-y_range(1))/(y_range(2)-y_range(1))*pos(4);

annotation('arrow', [norm_x1 norm_x2], [norm_y1 norm_y2], ...
    'Color', [0 0.6 0], 'LineWidth', 2.0, 'HeadWidth', 10, 'HeadLength', 10);
text(ax2+0.08, ay2+0.010, 'V 감소 방향', ...
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
        'g_5=0  K_2 상한 (K_2 \leq 1.6)', ...
        '설계변수 상·하한', ...
        sprintf('최적해  V^* = %.0f mm³', V_opt), ...
        sprintf('초기점  V_0 = %.0f mm³', V0)}, ...
    'Location','northwest', 'FontSize', 8.5);

xlim(x_range); ylim(y_range);
grid on; box on;
set(gca, 'FontSize', 10);
hold off;
