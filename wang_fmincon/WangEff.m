function eta = WangEff(x)
% WangEff.m
% K-H-V형 사이클로이드 감속기 전달 효율 계산
%
% 효율 구성:
%   eta = eta_x * eta_zx * eta_gx^2 * eta_sx
%
%   eta_x  : 핀 기어-사이클로이드 기어 맞물림 효율
%   eta_zx : 피벗 암 베어링 효율 (고정 상수)
%   eta_gx : 롤링 베어링 효율    (고정 상수)
%   eta_sx : 출력부(핀-핀홀) 효율
%
% 입력:
%   x = [Dp, drp, B, D, K1, Dw, dsw]

    Dp  = x(1);
    drp = x(2);
    K1  = x(5);
    Dw  = x(6);
    dsw = x(7);

    %% 고정 상수 (WangCon_fmincon, WangObj_fmincon과 동일하게 유지)
    zc     = 43;
    mu     = 0.05;    % 핀-기어 마찰계수 (논문 범위 0.05~0.1, 우리 채택값)
    fw     = 0.02;    % 출력부 마찰계수  (논문 범위 0.008~0.08, 우리 채택값)
    eta_zx = 0.99;    % 피벗 암 베어링 효율
    eta_gx = 0.995;   % 롤링 베어링 효율
    Delta2 = 2;       % 핀 슬리브 벽 두께 (mm)

    %% 맞물림 효율 eta_x
    % eta_nx: 가상 순수 회전 기어의 효율
    % eta_x:  행성 기어 구조 환산 효율
    eta_nx = 1 - (Dp - drp) * 4 * mu / (K1 * zc * Dp * pi);
    eta_x  = eta_nx / (1 + zc * (1 - eta_nx));

    %% 출력부 효율 eta_sx
    eta_sx = 1 - 4 * fw * K1 * dsw * Dp / (pi * Dw * (dsw + 2 * Delta2));

    %% 전체 효율
    eta = eta_x * eta_zx * eta_gx^2 * eta_sx;

    % 비물리적 값 방어 (0~1 범위 강제)
    if ~isfinite(eta) || eta <= 0
        eta = 1e-6;
    elseif eta > 1
        eta = 1;
    end
end
