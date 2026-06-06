% GearRun_fmincon.m
% fmincon(SQP)으로 사이클로이드 감속기 체적 최소화
% Król LM/SD 결과와 비교 + KKT 라그랑지 승수 분석

clc;

%% 1. 시작점 (논문 variant 1)
x0 = [2.8, 3.5, 0.7, 10, 5, 8, 30]';
% [e, q, m, h, Rs, Rh, Rw]

%% 2. 변수 경계
lb = [0.5,  0.5,  0.5,  1.0,  3.0,  3.1,  5.0 ]';
ub = [5.0,  8.0,  0.85, 20,   10,   20,   50  ]';

%% 3. fmincon 옵션
options = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxIterations',          500, ...
    'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance',    1e-6, ...
    'ConstraintTolerance',    1e-4);   % 약간 완화

%% 4. 고정 상수
zs=15; zk=16; zi=8; N=1000;
Mh=1550; nu1=0.3; nu2=0.3; Emod1=200000; Emod2=200000;

r0 = x0(1)*zk/x0(3);
V0 = GearVolume(r0, x0(1),zk,zs, x0(2),x0(3),N, x0(4),x0(5),x0(6),zi);
fprintf('=== 최적화 전 ===\n');
fprintf('체적 V0 = %.0f mm³\n', V0);
fprintf('설계변수: e=%.3f  q=%.3f  m=%.3f  h=%.3f  Rs=%.3f  Rh=%.3f  Rw=%.3f\n\n', x0);

%% 5. fmincon 실행
fprintf('=== fmincon (SQP) 실행 중 ===\n');
[xopt, Vopt, exitflag, output, lambda] = fmincon( ...
    @GearObj_fmincon, x0, [], [], [], [], lb, ub, @GearCon_fmincon, options);

%% 6. 결과 출력
fprintf('\n=== 최적화 결과 ===\n');
fprintf('체적 V* = %.0f mm³  (%.1f%% 감소)\n', Vopt, 100*(V0-Vopt)/V0);
fprintf('설계변수:\n');
labels_x = {'e','q','m','h','Rs','Rh','Rw'};
units_x  = {'mm','mm','-','mm','mm','mm','mm'};
for i=1:7
    fprintf('  %s = %.4f %s\n', labels_x{i}, xopt(i), units_x{i});
end
fprintf('종료 조건 (exitflag): %d', exitflag);
switch exitflag
    case  1; fprintf(' → KKT 조건 만족 (정상 수렴)\n');
    case  0; fprintf(' → 최대 반복 초과\n');
    case -2; fprintf(' → 가용 영역 없음 (제약 불만족)\n');
    otherwise; fprintf('\n');
end

%% 7. 제약 위반 확인
alp_opt = GearAlpha(xopt(1),zk,zs,xopt(3),xopt(2),0);
rho_opt = GearRho(xopt(1),zk,zs,xopt(3),xopt(2),alp_opt);
if rho_opt < 0
    lobe_alp = alp_opt; pit_alp = alp_opt + pi/zs;
else
    lobe_alp = alp_opt + pi/zs; pit_alp = alp_opt;
end
sig_pos = real(GearStressExPos(xopt(1),zk,zs,xopt(3),xopt(2),xopt(4),Mh,nu1,Emod1,nu2,Emod2,lobe_alp));
sig_neg = real(GearStressExNeg(xopt(1),zk,zs,xopt(3),xopt(2),xopt(4),Mh,nu1,Emod1,nu2,Emod2,pit_alp));
sig_in  = real(GearStressIn(zs,zi,xopt(4),Mh,xopt(7),xopt(5),xopt(6),nu1,Emod1,nu2,Emod2));

fprintf('\n=== 접촉응력 (한계 400 MPa) ===\n');
fprintf('  lobe (σ_EX+) = %.2f MPa  %s\n', sig_pos, okcheck(sig_pos,400));
fprintf('  pit  (σ_EX-) = %.2f MPa  %s\n', sig_neg, okcheck(sig_neg,400));
fprintf('  내부  (σ_IN) = %.2f MPa  %s\n', sig_in,  okcheck(sig_in, 400));

%% 8. KKT 라그랑지 승수 분석
fprintf('\n=== KKT 라그랑지 승수 ===\n');
constraint_labels = { ...
    'c(1)  편심률 하한 (언더컷 방지)', ...
    'c(2)  단폭계수 하한 m>=0.5', ...
    'c(3)  단폭계수 상한 m<=0.85', ...
    'c(4)  두께 하한 h>=0.2', ...
    'c(5)  내부슬리브 반경 하한 Rs>=3', ...
    'c(6)  홀반경 >= 슬리브반경 Rh>=Rs', ...
    'c(7)  외부슬리브 간섭 방지 q<=|rho|', ...
    'c(8)  Rw 하한 Rw>=2*Rs', ...
    'c(9)  Rw 상한 Rw<=r-Rs', ...
    'c(10) lobe 곡률반경 하한 >=9mm', ...
    'c(11) lobe 곡률반경 상한 <=100mm', ...
    'c(12) pit  곡률반경 하한 >=-100mm', ...
    'c(13) pit  곡률반경 상한 <=-2mm', ...
    'c(14) lobe 접촉응력 <=400MPa', ...
    'c(15) pit  접촉응력 <=400MPa', ...
    'c(16) 내부슬리브 접촉응력 <=400MPa'};

lam_ineq = lambda.ineqnonlin;
active_found = false;
for i = 1:length(lam_ineq)
    if abs(lam_ineq(i)) > 1e-3
        fprintf('  [활성] %s  λ=%.4f\n', constraint_labels{i}, lam_ineq(i));
        active_found = true;
    end
end
if ~active_found
    fprintf('  활성 제약 없음 (모든 λ ≈ 0)\n');
end
fprintf('\n  λ > 0 : 활성 제약 (현재 설계를 제한하는 제약)\n');
fprintf('  λ 클수록 : 그 제약을 완화하면 체적을 더 줄일 수 있음\n');

% lb/ub 라그랑지 승수도 확인
fprintf('\n=== 변수 경계 라그랑지 승수 ===\n');
lam_lb = lambda.lower;
lam_ub = lambda.upper;
for i=1:7
    if abs(lam_lb(i)) > 1e-3
        fprintf('  [하한 활성] %s >= %.2f  λ=%.4f\n', labels_x{i}, lb(i), lam_lb(i));
    end
    if abs(lam_ub(i)) > 1e-3
        fprintf('  [상한 활성] %s <= %.2f  λ=%.4f\n', labels_x{i}, ub(i), lam_ub(i));
    end
end

%% 9. Król LM 결과와 비교
fprintf('\n=== Król LM 결과와 비교 (variant 1) ===\n');
xLM = [2.3317, 4.1502, 0.6238, 8.2701, 3.5290, 12.4852, 34.784]';
rLM = xLM(1)*zk/xLM(3);
VLM = GearVolume(rLM,xLM(1),zk,zs,xLM(2),xLM(3),N,xLM(4),xLM(5),xLM(6),zi);
fprintf('  방법         체적(mm³)   시작점 대비 감소\n');
fprintf('  시작점       %8.0f\n', V0);
fprintf('  Król LM      %8.0f   %.1f%%\n', VLM, 100*(V0-VLM)/V0);
fprintf('  fmincon SQP  %8.0f   %.1f%%\n', Vopt, 100*(V0-Vopt)/V0);

%% 10. 기어 형상 (GearDraw는 인자 2개: xfin, xp)
figure;
GearDraw(xopt, x0);   % 파랑=최적화 후, 빨강=전
title('fmincon 최적화 전(빨강) / 후(파랑)');
axis equal; grid on;

%% 보조 함수
function s = okcheck(val, limit)
    if val <= limit; s = '✓ OK'; else; s = '✗ 초과'; end
end
