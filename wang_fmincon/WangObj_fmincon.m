function V = WangObj_fmincon(x)
% Wang et al.(2016) fmincon용 목적함수 — 체적 (식 10)
% x = [Dp, drp, B, D, K1, Dw, dsw]
%
% 감속기 사양 (XW3형, 전달비 43):
%   zc=43, zp=44, zw=8
%   Delta2=2 mm (핀 슬리브 벽 두께, 논문 미명시 → 일반값 사용)

    Dp  = x(1);   % 핀 중심원 직경
    drp = x(2);   % 핀 직경
    B   = x(3);   % 사이클로이드 기어 폭
    D   = x(4);   % 사이클로이드 기어 중심홀 직경
    K1  = x(5);   % 단폭계수
    Dw  = x(6);   % 출력 핀 중심원 직경
    dsw = x(7);   % 출력 핀 직경

    zp     = 44;   % 핀 기어 수 (= zc + 1)
    zc     = 43;   % 사이클로이드 기어 치형 수
    zw     = 8;    % 출력 핀 수
    Delta2 = 2;    % 핀 슬리브 벽 두께 (mm)

    % 식 (10): 체적
    V = (1/4) * pi * B * ( ...
          (Dp - K1*(Dp/zp) - drp)^2 ...
        - (dsw + 2*Delta2 + K1*(Dp/zp))^2 * zw ...
        - D^2 ) ...
        + K1*(Dp/zp)*zc*B;

    % 체적이 비물리적이면 큰 값 반환 (수치 안정)
    if ~isfinite(V) || V <= 0
        V = 1e9;
    end
end
