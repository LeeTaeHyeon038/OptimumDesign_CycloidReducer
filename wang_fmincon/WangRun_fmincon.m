% WangRun_fmincon.m
% Wang et al.(2016) 수식 기반 fmincon(SQP) 체적 최소화
%
% 목적: Król LM/SD와 다른 접근법(SQP + 명시적 제약)으로
%       수업의 KKT 조건·라그랑지 승수·후최적성 해석 실습
%
% 비교 기준: Wang 논문 Table 6 (단목적 체적 최소화)
%   - 초기 체적: 101,204 mm³
%   - 최적 체적: 64,693 mm³ (GA 기준, -36.1%)

clc;

%% 1. 초기값 및 변수 범위 (논문 Table 4)
% x = [Dp, drp, B, D, K1, Dw, dsw]
x0 = [144,  10,   11,  53.5, 0.6069, 90,  12 ]';
lb = [140, 7, 7, 50, 0.60, 88, 11]';
ub = [155,  10.4, 12,  55,   0.9,   104,  14 ]';

%% 2. 상수 (WangCon_fmincon과 동일하게 유지)
zc = 43; zp = 44; zw = 8;
Delta2 = 2; n = 1440; P = 0.75;
M  = 9550 * P / (n/zc);   % 출력 토크 (N·mm)

%% 3. 초기 설계 확인
V0 = WangObj_fmincon(x0);
[c0, ~] = WangCon_fmincon(x0);
fprintf('=== 초기 설계 (Wang Table 4) ===\n');
fprintf('체적 V0 = %.0f mm³  (논문값: 101,204)\n', V0);
fprintf('초기점 제약 위반 여부:\n');
any_viol = false;
for i = 1:16
    if c0(i) > 1e-6
        fprintf('  c(%d) = %.4f  위반\n', i, c0(i));
        any_viol = true;
    end
end
if ~any_viol
    fprintf('  모든 제약 만족 → 가용 시작점 확인\n');
end

%% 4. fmincon 옵션 (SQP)
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-6);

%% 5. fmincon 실행
fprintf('\n=== fmincon (SQP) 실행 중 ===\n');
[xopt, Vopt, exitflag, output, lambda] = fmincon( ...
    @WangObj_fmincon, x0, [], [], [], [], lb, ub, @WangCon_fmincon, options);

%% 6. 결과 출력
fprintf('\n=== 최적화 결과 ===\n');
fprintf('체적 V* = %.0f mm³  (%.1f%% 감소)\n', Vopt, 100*(V0-Vopt)/V0);
fprintf('\n설계변수 비교:\n');
names = {'Dp (핀중심원직경)', 'drp (핀직경)', 'B (기어폭)', ...
         'D (중심홀직경)', 'K1 (단폭계수)', 'Dw (출력핀중심원)', 'dsw (출력핀직경)'};
units = {'mm','mm','mm','mm','-','mm','mm'};
fprintf('  %-22s  %8s  %8s\n', '변수', '초기값', '최적값');
fprintf('  %s\n', repmat('-', 1, 42));
for i = 1:7
    fprintf('  %-22s  %8.4f  %8.4f  %s\n', names{i}, x0(i), xopt(i), units{i});
end

fprintf('\n종료 조건 (exitflag): %d', exitflag);
switch exitflag
    case  1; fprintf(' → KKT 조건 만족 (정상 수렴)\n');
    case  2; fprintf(' → 스텝 크기 허용오차 수렴\n');
    case  0; fprintf(' → 최대 반복 초과\n');
    case -2; fprintf(' → 가용 영역 없음\n');
    otherwise; fprintf('\n');
end
fprintf('반복 횟수: %d\n', output.iterations);
fprintf('함수 평가 횟수: %d\n', output.funcCount);

%% 7. 효율 계산 (참고값)
mu = 0.05; fw = 0.02; eta_gx = 0.995; eta_zx = 0.99;
Dp=xopt(1); drp=xopt(2); K1=xopt(5); Dw=xopt(6); dsw=xopt(7);
eta_nx = 1 - (Dp-drp)*4*mu / (K1*zc*Dp*pi);
eta_x  = eta_nx / (1 + zc*(1-eta_nx));
eta_sx = 1 - 4*fw*K1*dsw*Dp / (pi*Dw*(dsw+2*Delta2));
eta    = eta_x * eta_zx * eta_gx^2 * eta_sx;

Dp0=x0(1); drp0=x0(2); K10=x0(5); Dw0=x0(6); dsw0=x0(7);
eta_nx0 = 1 - (Dp0-drp0)*4*mu / (K10*zc*Dp0*pi);
eta_x0  = eta_nx0 / (1 + zc*(1-eta_nx0));
eta_sx0 = 1 - 4*fw*K10*dsw0*Dp0 / (pi*Dw0*(dsw0+2*Delta2));
eta0    = eta_x0 * eta_zx * eta_gx^2 * eta_sx0;

fprintf('\n효율 (참고):\n');
fprintf('  초기: η = %.4f\n', eta0);
fprintf('  최적: η = %.4f  (%.2f%%p 변화)\n', eta, (eta-eta0)*100);

%% 8. Wang 논문 결과와 비교
fprintf('\n=== Wang 논문 Table 6 결과와 비교 (단목적 체적 최소화) ===\n');
fprintf('  방법              체적(mm³)   초기 대비 감소\n');
fprintf('  초기 설계         %8.0f\n', V0);
fprintf('  Wang GA (표준)    %8d   -36.1%%\n', 64693);
fprintf('  fmincon SQP       %8.0f   %.1f%%\n', Vopt, 100*(V0-Vopt)/V0);

%% 9. KKT 라그랑지 승수 분석
fprintf('\n=== KKT 라그랑지 승수 분석 ===\n');
fprintf('활성 제약 (λ > 1e-4):\n');

constraint_names = { ...
    'y1  언더컷/첨예 방지', ...
    'y2  단폭계수 하한 K1>=0.65', ...
    'y3  단폭계수 상한 K1<=0.9', ...
    'y4  핀직경계수 K2 하한', ...
    'y5  핀직경계수 K2 상한', ...
    'y6  사이클로이드 접촉강도 <=1000MPa', ...
    'y7  핀기어 굽힘강도 <=150MPa', ...
    'y8  핀-핀홀 접촉강도 <=400MPa', ...
    'y9  핀 굽힘강도 <=200MPa', ...
    'y10 Dp 하한 >=140mm', ...
    'y11 Dp 상한 <=155mm', ...
    'y12 핀홀 직경 조건 1', ...
    'y13 핀홀 직경 조건 2', ...
    'y14 기어폭 하한 B>=0.05Dp', ...
    'y15 기어폭 상한 B<=0.1Dp', ...
    'y16 피벗 베어링 수명 >=5000h'};

lam = lambda.ineqnonlin;
active_any = false;
for i = 1:length(lam)
    if abs(lam(i)) > 1e-4
        fprintf('  [활성] c(%2d) %s\n         λ = %.6f\n', i, constraint_names{i}, lam(i));
        active_any = true;
    end
end
if ~active_any
    fprintf('  활성 제약 없음 (모든 λ ≈ 0)\n');
end

% lb/ub 라그랑지 승수
fprintf('\n경계 활성 변수 (λ_lb 또는 λ_ub > 1e-4):\n');
lam_lb = lambda.lower;
lam_ub = lambda.upper;
bound_any = false;
for i = 1:7
    if abs(lam_lb(i)) > 1e-4
        fprintf('  [하한 활성] %s >= %.4f  λ=%.6f\n', names{i}, lb(i), lam_lb(i));
        bound_any = true;
    end
    if abs(lam_ub(i)) > 1e-4
        fprintf('  [상한 활성] %s <= %.4f  λ=%.6f\n', names{i}, ub(i), lam_ub(i));
        bound_any = true;
    end
end
if ~bound_any
    fprintf('  경계 활성 변수 없음\n');
end

fprintf('\n해석:\n');
fprintf('  λ > 0  : 활성 제약 (현재 해를 제한하고 있음)\n');
fprintf('  λ 클수록: 그 제약을 완화하면 체적을 더 줄일 수 있음\n');
fprintf('  λ = 0  : 비활성 제약 (설계에 여유 있음)\n');

%% 10. 최적해 제약 만족 여부 확인
fprintf('\n=== 최적해 제약 위반 확인 ===\n');
[c_opt, ~] = WangCon_fmincon(xopt);
any_viol2 = false;
for i = 1:16
    if c_opt(i) > options.ConstraintTolerance
        fprintf('  c(%2d) = %.6f  위반\n', i, c_opt(i));
        any_viol2 = true;
    end
end
if ~any_viol2
    fprintf('  모든 제약 만족\n');
end
