function [c, ceq] = WangCon_fmincon(x)
% Wang et al.(2016) fmincon용 제약함수
% 16개 부등호 제약 c(x) <= 0, 등호 제약 없음
% x = [Dp, drp, B, D, K1, Dw, dsw]
%
% 제약 목록:
%   c(1)        y1  : 언더컷/첨예 방지 (식 12)
%   c(2), c(3)  y2,y3: 단폭계수 K1 범위 (식 13, 14)
%   c(4), c(5)  y4,y5: 핀직경 계수 K2 범위 (식 16, 17)
%   c(6)        y6  : 사이클로이드 기어 접촉강도 (식 22)
%   c(7)        y7  : 핀 기어 굽힘강도 (식 25)
%   c(8)        y8  : 핀-핀홀 접촉강도 (식 27)
%   c(9)        y9  : 핀 굽힘강도 (식 29)
%   c(10),c(11) y10,y11: Dp 범위 (식 30, 31)
%   c(12),c(13) y12,y13: 핀홀 최대 직경 (식 32, 33)
%   c(14),c(15) y14,y15: 기어 폭 B 범위 (식 35, 36)
%   c(16)       y16 : 피벗 베어링 수명 (식 38)

    Dp  = x(1);
    drp = x(2);
    B   = x(3);
    D   = x(4);
    K1  = x(5);
    Dw  = x(6);
    dsw = x(7);

    %% 고정 상수 (XW3형, 전달비 43)
    zc     = 43;    % 사이클로이드 치형 수
    zp     = 44;    % 핀 기어 수
    zw     = 8;     % 출력 핀 수
    Delta2 = 2;     % 핀 슬리브 벽 두께 (mm)
    Delta  = 0.15;  % 핀홀-슬리브 간격 (mm)
    a      = 0;     % 중심 거리 (근사: 0으로 처리)
    Deltac = 2;     % 스페이서 링 두께 (mm, 일반값)
    Kw     = 1.4;   % 하중 분배 계수
    n      = 1440;  % 입력 회전수 (rpm)
    P      = 0.75;  % 입력 동력 (kW)
    M      = 9550 * P / (n / zc);   % 출력 토크 (N·mm): M = 9550*P*i/n * 1000
    L      = Dp / 2;  % 핀 지지 간격 (근사: Dp/2, mm)
    Ee     = 2.19e5;  % 등가 탄성계수 GCr15 강재 (MPa)

    % 허용 응력 (논문 Section 4)
    sigma_HP  = 1000;   % 접촉 허용응력 (MPa)
    sigma_FP  = 150;    % 굽힘 허용응력 핀기어 (MPa)
    sigma_ZHP = 400;    % 핀-핀홀 접촉 허용응력 (MPa, 일반 강재 기준)
    sigma_BP  = 200;    % 핀 굽힘 허용응력 (MPa)

    % 단폭계수 범위 (Table 1, zc=25~59)
    c1 = 0.60;  c2 = 0.9;   % K1 초기값 0.6069 포함하도록 완화

    % 핀직경 계수 범위 (Table 2, zp=36~60)
    d1 = 1.0;   d2 = 1.6;   % K2 하한/상한 순서 수정

    % Dp 범위 (Table 3, type 3: 140~155)
    e1 = 140;   e2 = 155;

    %% 최소 등가 곡률반경 rho_emin (접촉강도 제약용)
    % 식 (11) — K1 범위에 따라 두 경우
    threshold = (zc - 1) / (2*zc + 1);
    if K1 > threshold
        rho_emin = (Dp/2) * sqrt(27*(1-K1^2)*(zp-1) / (zp+1)^3);
    else
        rho_emin = Dp*(1-K1)^2 / (2*(zp*K1+1));
    end
    rho_emin = max(rho_emin, 1e-6);  % 수치 안정

    %% 핀홀 직경 dw (식 34)
    dw = dsw + 2*Delta2 + 2*a + Delta;

    %% 피벗 베어링 수명 계산 (식 37, 38)
    R  = 0.825 * M * zp / (K1 * Dp * zc);   % 공칭 반경 하중
    p  = 1.25 * R;                            % 동적 하중
    n1 = n * zp / (zp - 1);                  % 피벗 베어링 회전수
    Rp = R;
    C  = 9.8 * (6.4 * 0.1 * Rp)^1.85 * (2.065 + 0.91 * Rp^0.75);  % 정격 동하중
    if p > 0 && C > 0
        Lh = (1e6 / (60*n1)) * (C/p)^(10/3);
    else
        Lh = 0;
    end

    %% 제약 벡터 c(x) <= 0
    c = zeros(16, 1);

    % c(1): y1 — 언더컷/첨예 방지 (식 12)
    if K1 > threshold
        c(1) = drp/2 - (Dp/2)*sqrt(27*(1-K1^2)*(zp-1)/(zp+1)^3);
    else
        c(1) = drp/2 - Dp*(1-K1)^2 / (2*(zp*K1+1));
    end

    % c(2), c(3): y2,y3 — 단폭계수 범위 (식 13, 14)
    c(2) = c1 - K1;
    c(3) = K1 - c2;

    % c(4), c(5): y4,y5 — 핀직경 계수 K2 범위 (식 16, 17)
    K2 = (Dp/drp) * sin(pi/zp);
    c(4) = d1 - K2;
    c(5) = K2 - d2;

    % c(6): y6 — 사이클로이드 기어 접촉강도 (식 22)
    c(6) = 0.418 * sqrt(Ee * 4.4*M / (B * K1 * zc * Dp * rho_emin)) - sigma_HP;

    % c(7): y7 — 핀 기어 굽힘강도 (식 25)
    if Dp < 390
        c(7) = (1.41 * 4.4 * 9550*P*L) / (K1 * Dp * n * drp^2) - sigma_FP;
    else
        c(7) = (0.48 * 4.4 * 9550*P*L) / (K1 * Dp * n * drp^2) - sigma_FP;
    end

    % c(8): y8 — 핀-핀홀 접촉강도 (식 27)
    term8 = K1*M*Dp / (zw*Dw*B*(dsw/2+Delta2)^2*zp + 0.5*K1*Dp*(dsw/2+Delta2));
    if term8 >= 0
        c(8) = 300*sqrt(term8) - sigma_ZHP;
    else
        c(8) = 300*sqrt(abs(term8)) - sigma_ZHP;  % 방어
    end

    % c(9): y9 — 핀 굽힘강도 (식 29)
    c(9) = (4.4 * Kw * M * (1.5*B + Deltac)) / (0.1 * zw * Dw * dsw^3) - sigma_BP;

    % c(10), c(11): y10,y11 — Dp 범위 (식 30, 31)
    c(10) = e1 - Dp;
    c(11) = Dp - e2;

    % c(12), c(13): y12,y13 — 핀홀 최대 직경 (식 32, 33)
    c(12) = 0.06*Dp - Dw + dw + D;
    c(13) = 0.03*Dp - Dw*sin(pi/zw) + dw;

    % c(14), c(15): y14,y15 — 기어 폭 B 범위 (식 35, 36)
    c(14) = 0.05*Dp - B;
    c(15) = B - 0.1*Dp;

    % c(16): y16 — 피벗 베어링 수명 >= 5000h (식 38)
    c(16) = 5000 - Lh;

    ceq = [];
end
