% WangRun_fmincon.m  (v2 — Figure 1~3 시각화 추가)
% Wang et al.(2016) 수식 기반 fmincon(SQP) 체적 최소화
%
% 생성 Figure:
%   Figure 1: 수렴 이력 (체적 변화 및 초기→최적 비교)
%   Figure 2: 설계변수 비교 (초기값 vs 최적값)
%   Figure 3: KKT 라그랑지 승수 (활성 제약 및 경계 활성 변수)

clc;

%% 1. 초기값 및 변수 범위
x0 = [144,  10,   11,  53.5, 0.6069, 90,  12 ]';
lb = [140,   7,    7,  50,   0.60,   88,  11 ]';
ub = [155,  10.4, 12,  55,   0.90,  104,  14 ]';

names     = {'Dp','drp','B','D','K1','Dw','dsw'};
names_kor = {'D_p','d_{rp}','B','D','K_1','D_w','d_{sw}'};
units     = {'mm','mm','mm','mm','-','mm','mm'};

constraint_names = { ...
    'y1  언더컷/첨예 방지', 'y2  K1 하한', 'y3  K1 상한', ...
    'y4  K2 하한', 'y5  K2 상한', 'y6  접촉강도', ...
    'y7  핀 굽힘강도', 'y8  핀홀 접촉강도', 'y9  출력핀 굽힘강도', ...
    'y10 Dp 하한', 'y11 Dp 상한', 'y12 핀홀 조건1', ...
    'y13 핀홀 조건2', 'y14 B 하한', 'y15 B 상한', 'y16 베어링 수명'};

%% 2. 초기 설계 성능
V0   = WangObj_fmincon(x0);
eta0 = WangEff(x0);

fprintf('=== 초기 설계 ===\n');
fprintf('체적 V0 = %.0f mm³\n', V0);
fprintf('효율 η0 = %.4f\n\n', eta0);

%% 3. fmincon 옵션 (SQP)
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-6);

%% 4. fmincon 실행
fprintf('=== fmincon (SQP) 실행 중 ===\n');
[xopt, Vopt, exitflag, output, lambda] = fmincon( ...
    @WangObj_fmincon, x0, [], [], [], [], lb, ub, @WangCon_fmincon, options);

eta_opt = WangEff(xopt);

%% 5. 결과 출력
fprintf('\n=== 최적화 결과 ===\n');
fprintf('체적 V* = %.0f mm³  (%.1f%% 감소)\n', Vopt, 100*(V0-Vopt)/V0);
fprintf('효율 η* = %.4f  (%+.4f)\n\n', eta_opt, eta_opt-eta0);

fprintf('  %-20s  %8s  %8s\n', '변수', '초기값', '최적값');
fprintf('  %s\n', repmat('-', 1, 40));
for i = 1:7
    fprintf('  %-20s  %8.4f  %8.4f  %s\n', names_kor{i}, x0(i), xopt(i), units{i});
end

fprintf('\nexitflag: %d', exitflag);
switch exitflag
    case  1; fprintf(' → KKT 조건 만족\n');
    case  2; fprintf(' → 스텝 크기 수렴\n');
    case  0; fprintf(' → 최대 반복 초과\n');
    case -2; fprintf(' → 가용 영역 없음\n');
    otherwise; fprintf('\n');
end
fprintf('반복: %d회,  함수 평가: %d회\n', output.iterations, output.funcCount);

%% 6. KKT 라그랑지 승수
lam    = lambda.ineqnonlin;
lam_lb = lambda.lower;
lam_ub = lambda.upper;

fprintf('\n=== KKT 라그랑지 승수 ===\n');
for i = 1:16
    if abs(lam(i)) > 1e-4
        fprintf('  [활성] c(%2d) %-22s  λ = %.2f\n', i, constraint_names{i}, lam(i));
    end
end
for i = 1:7
    if abs(lam_lb(i)) > 1e-4
        fprintf('  [하한] %-6s  λ = %.4f\n', names{i}, lam_lb(i));
    end
    if abs(lam_ub(i)) > 1e-4
        fprintf('  [상한] %-6s  λ = %.4f\n', names{i}, lam_ub(i));
    end
end

%% 7. 제약 만족 확인
[c_opt, ~] = WangCon_fmincon(xopt);
if ~any(c_opt > options.ConstraintTolerance)
    fprintf('\n최적해: 모든 제약 만족\n');
end

%% 8. Wang 논문 결과와 비교
fprintf('\n=== Wang 논문 Table 6 결과와 비교 ===\n');
fprintf('  초기 설계:   %8.0f mm³\n', V0);
fprintf('  Wang GA:     %8d mm³  (-36.1%%)\n', 64693);
fprintf('  fmincon SQP: %8.0f mm³  (%.1f%%)\n', Vopt, 100*(V0-Vopt)/V0);

%% ══════════════════════════════════════════════════════════
%%  Figure 1: 단목적 SQP 최적화 결과
%% ══════════════════════════════════════════════════════════
figure('Name','Figure 1: 단목적 SQP 최적화 결과','NumberTitle','off');

subplot(1,2,1);
% categorical 대신 숫자 x축 사용
bar_vals = [V0, Vopt] / 1e3;
bar(1:2, bar_vals, 0.5, 'FaceColor',[0.4 0.65 0.9], 'EdgeColor','w');
set(gca, 'XTick', 1:2, 'XTickLabel', {'초기 설계','최적 설계'}, 'FontSize',10);
ylabel('체적 V (×10³ mm³)');
title(sprintf('(a) 체적 감소: %.0f → %.0f mm³\n(%.1f%% 감소, %d회 반복)', ...
    V0, Vopt, 100*(V0-Vopt)/V0, output.iterations));
ylim([0, V0/1e3 * 1.15]);
text(1, V0/1e3*1.02,   sprintf('%.0f', V0),   'HorizontalAlignment','center','FontSize',10);
text(2, Vopt/1e3*1.02, sprintf('%.0f', Vopt),  'HorizontalAlignment','center','FontSize',10);
grid on;

subplot(1,2,2);
% 세 방법 비교 (초기, Wang GA, SQP)
% categorical 대신 숫자 x축 사용 → 자동 정렬 문제 방지
vals2   = [V0, 64693, Vopt] / 1e3;   % 초기, Wang GA, SQP 순
colors2 = [0.7 0.7 0.7; 0.4 0.75 0.4; 0.4 0.65 0.9];
labels2 = {'초기 설계', 'Wang GA', 'SQP 최적'};

b2_ax = bar(1:3, vals2, 0.6);
b2_ax.FaceColor = 'flat';
b2_ax.CData     = colors2;
b2_ax.EdgeColor = 'w';
set(gca, 'XTick', 1:3, 'XTickLabel', labels2, 'FontSize',10);
ylabel('체적 V (×10³ mm³)');
title('(b) 방법별 결과 비교');
ylim([0, V0/1e3 * 1.15]);
for k = 1:3
    text(k, vals2(k)*1.02, sprintf('%.0f', vals2(k)*1e3), ...
        'HorizontalAlignment','center','FontSize',9);
end
grid on;

sgtitle('Figure 1: 단목적 SQP 최적화 결과');

%% ══════════════════════════════════════════════════════════
%%  Figure 2: 설계변수 비교
%% ══════════════════════════════════════════════════════════
figure('Name','Figure 2: 설계변수 비교','NumberTitle','off');

x_ax = 1:7;
bw   = 0.35;

hold on;
b1 = bar(x_ax - bw/2, x0',   bw, 'FaceColor',[0.5 0.7 1.0], 'EdgeColor','w');
b2 = bar(x_ax + bw/2, xopt', bw, 'FaceColor',[1.0 0.5 0.3], 'EdgeColor','w');
% 상·하한선 — 첫 번째 것만 변수로 받아서 범례에 사용
h_lb = plot([1-bw*1.2, 1+bw*1.2], [lb(1), lb(1)], 'r--', 'LineWidth',1.2);
h_ub = plot([1-bw*1.2, 1+bw*1.2], [ub(1), ub(1)], 'g--', 'LineWidth',1.2);
for i = 2:7
    plot([i-bw*1.2, i+bw*1.2], [lb(i), lb(i)], 'r--', 'LineWidth',1.2);
    plot([i-bw*1.2, i+bw*1.2], [ub(i), ub(i)], 'g--', 'LineWidth',1.2);
end
hold off;

set(gca, 'XTick', x_ax, 'XTickLabel', names_kor, 'FontSize',10);
ylabel('설계변수 값 (mm 또는 무차원)');
legend([b1, b2, h_lb, h_ub], {'초기값','최적값','하한 (LB)','상한 (UB)'}, 'Location','northeast');
title('Figure 2: 단목적 최적화 — 설계변수 초기값 vs 최적값');
text(7.3, lb(7), ' LB', 'Color','r', 'FontSize',8);
text(7.3, ub(7), ' UB', 'Color','g', 'FontSize',8);
grid on;

%% ══════════════════════════════════════════════════════════
%%  Figure 3: KKT 라그랑지 승수
%% ══════════════════════════════════════════════════════════
all_labels = {};
all_vals   = [];

for i = 1:16
    if abs(lam(i)) > 1e-4
        all_labels{end+1} = constraint_names{i};
        all_vals(end+1)   = lam(i);
    end
end
for i = 1:7
    if abs(lam_lb(i)) > 1e-4
        all_labels{end+1} = [names{i} ' (하한)'];
        all_vals(end+1)   = lam_lb(i);
    end
    if abs(lam_ub(i)) > 1e-4
        all_labels{end+1} = [names{i} ' (상한)'];
        all_vals(end+1)   = lam_ub(i);
    end
end

[all_vals_s, sidx] = sort(all_vals, 'descend');
all_labels_s = all_labels(sidx);

figure('Name','Figure 3: KKT 라그랑지 승수','NumberTitle','off');

barh(all_vals_s, 'FaceColor',[0.3 0.65 0.5], 'EdgeColor','w');
set(gca, 'YTick', 1:length(all_vals_s), ...
         'YTickLabel', all_labels_s, 'FontSize',10);
xlabel('라그랑지 승수 λ');
title('Figure 3: KKT 라그랑지 승수 — 활성 제약 및 경계 활성 변수');
grid on;
for i = 1:length(all_vals_s)
    text(all_vals_s(i)*1.01, i, sprintf(' %.0f', all_vals_s(i)), ...
        'FontSize',9, 'VerticalAlignment','middle');
end
