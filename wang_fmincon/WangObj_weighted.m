function f = WangObj_weighted(x, wV, wEta, V0, eta0)
% WangObj_weighted.m
% 가중합(Weighted Sum) 다목적 목적함수
%
% 목적함수:
%   f = wV * (V/V0) + wEta * ((1-eta)/(1-eta0))
%
% 두 항 모두 무차원화(초기값으로 나눔)하여 스케일 통일
%   - 체적항: 클수록 나쁨 → 최소화
%   - 효율항: (1-eta)는 효율 손실 → 클수록 나쁨 → 최소화
%
% 입력:
%   x     : 설계변수 [Dp, drp, B, D, K1, Dw, dsw]
%   wV    : 체적 가중치    (0 <= wV <= 1)
%   wEta  : 효율 가중치   (wEta = 1 - wV)
%   V0    : 초기 체적 (무차원화 기준)
%   eta0  : 초기 효율 (무차원화 기준)

    V   = WangObj_fmincon(x);        % 체적 계산 (기존 파일 재사용)
    eta = WangEff(x);                % 효율 계산

    % 효율 손실 무차원화: (1-eta)/(1-eta0)
    % eta0 = 1이면 분모가 0이 되므로 방어
    if abs(1 - eta0) < 1e-12
        eta_term = 0;
    else
        eta_term = (1 - eta) / (1 - eta0);
    end

    f = wV * (V / V0) + wEta * eta_term;
end
