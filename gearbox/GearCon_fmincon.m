function [c, ceq] = GearCon_fmincon(x)
% fmincon용 제약함수 — Król 페널티와 동일한 로직으로 구현
% c(x) <= 0 형태 (16개 부등호 제약)
% x = [e, q, m, h, Rs, Rh, Rw]

    e  = x(1); q  = x(2); m  = x(3); h  = x(4);
    Rs = x(5); Rh = x(6); Rw = x(7);

    zs=15; zk=16; zi=8; Mh=1550; N=1000;
    nu1=0.3; nu2=0.3; Emod1=200000; Emod2=200000;
    r = e * zk / m;

    % 한계값 (GearGetConstant와 동일)
    MAX_STRESS = 600;   % case 1
    PMIN_RHO   = 7;     % case 2: positive(pit) 곡률반경 하한
    PMAX_RHO   = 100;   % case 3: positive(pit) 곡률반경 상한
    NMIN_RHO   = -2;    % case 4: negative(lobe) 곡률반경 상한 (절댓값 하한)
    NMAX_RHO   = -100;  % case 5: negative(lobe) 곡률반경 하한 (절댓값 상한)
    M_MIN=0.5; M_MAX=0.85; H_MIN=0.2; RS_MIN=3;
    E_MIN = q*(zk+1)/(3*sqrt(3)*zk)*sqrt((zk+1)/(zk-1))*sqrt(m^2/(1-m^2));

    % 현재 체적
    V_now = GearVolume(r, e, zk, zs, q, m, N, h, Rs, Rh, zi);

    %% 곡률반경 극값 탐색 — Król 코드와 동일한 로직
    alp = GearAlpha(e, zk, zs, m, q, 0);
    rho_alp = GearRho(e, zk, zs, m, q, alp);

    % Rhon (negative, lobe): rho<0이면 alp, 아니면 alp+pi/zs
    if rho_alp < 0
        n_alpha = alp;           % lobe 위치
        p_alpha = alp + pi/zs;   % pit 위치
    else
        n_alpha = alp + pi/zs;   % lobe 위치
        p_alpha = alp;           % pit 위치
    end

    % Rhop (positive, pit): rho<0이면 alp+pi/zs, 아니면 alp
    % (위와 동일, p_alpha = pit 위치)

    rho_n = GearRho(e, zk, zs, m, q, n_alpha);  % 음수 (lobe)
    rho_p = GearRho(e, zk, zs, m, q, p_alpha);  % 양수 (pit)

    % 접촉응력 계산 (lobe=n_alpha, pit=p_alpha)
    sig_pos = real(GearStressExPos(e,zk,zs,m,q,h,Mh,nu1,Emod1,nu2,Emod2,n_alpha));
    sig_neg = real(GearStressExNeg(e,zk,zs,m,q,h,Mh,nu1,Emod1,nu2,Emod2,p_alpha));
    sig_in  = real(GearStressIn(zs,zi,h,Mh,Rw,Rs,Rh,nu1,Emod1,nu2,Emod2));

    c = zeros(17,1);

    %% 박스 제약 — 체적 차이로 환산 (Król 방식)
    % c(1): 편심률 하한 e >= E_MIN
    if e < E_MIN
        c(1) = GearVolume(E_MIN*zk/m, E_MIN,zk,zs,q,m,N,h,Rs,Rh,zi) - V_now;
    end

    % c(2): 단폭계수 하한 m >= 0.5
    if m < M_MIN
        c(2) = GearVolume(e*zk/M_MIN,e,zk,zs,q,M_MIN,N,h,Rs,Rh,zi) - V_now;
    end

    % c(3): 단폭계수 상한 m <= 0.85
    if m > M_MAX
        c(3) = V_now - GearVolume(e*zk/M_MAX,e,zk,zs,q,M_MAX,N,h,Rs,Rh,zi);
    end

    % c(4): 두께 하한 h >= 0.2
    if h < H_MIN
        c(4) = GearVolume(r,e,zk,zs,q,m,N,H_MIN,Rs,Rh,zi) - V_now;
    end

    % c(5): 내부슬리브 반경 하한 Rs >= 3
    if Rs < RS_MIN
        c(5) = GearVolume(r,e,zk,zs,q,m,N,h,RS_MIN,Rh,zi) - V_now;
    end

    % c(6): 홀반경 >= 슬리브반경 Rh >= Rs
    if Rh < Rs
        c(6) = GearVolume(r,e,zk,zs,q,m,N,h,Rs,Rs,zi) - V_now;
    end

    %% 기하학적 제약
    % c(7): 외부슬리브 간섭 방지 q <= |rho_n| (GearPenaltyqmin과 동일)
    rho_min_abs = abs(rho_n);
    if q > rho_min_abs
        c(7) = GearVolume(r,e,zk,zs,q,m,N,h,Rs,Rh,zi) - ...
               GearVolume(r,e,zk,zs,rho_min_abs,m,N,h,Rs,Rh,zi);
    end

    % c(8): Rw 하한 Rw >= 2*Rs (GearPenaltyRwmin과 동일)
    if Rw < 2*Rs
        sig_rw_low = real(GearStressIn(zs,zi,h,Mh,2*Rs,Rs,Rh,nu1,Emod1,nu2,Emod2));
        c(8) = sig_rw_low - sig_in;
    end

    % c(9): Rw 상한 Rw <= r-Rs (GearPenaltyRwmax와 동일)
    if Rw > (r - Rs)
        sig_rw_max = real(GearStressIn(zs,zi,h,Mh,r-Rs,Rs,Rh,nu1,Emod1,nu2,Emod2));
        c(9) = sig_in - sig_rw_max;
    end

    %% 곡률반경 범위 — Král 페널티와 동일한 조건
    % c(10): GearPenaltyRhopMin — pit 곡률반경 하한: rho_p >= PMIN_RHO(9)
    if rho_p < PMIN_RHO
        c(10) = PMIN_RHO - rho_p;
    end

    % c(11): GearPenaltyRhopMax — pit 곡률반경 상한: rho_p <= PMAX_RHO(100)
    if rho_p > PMAX_RHO
        c(11) = abs(rho_p) - PMAX_RHO;
    end

    % c(12): GearPenaltyRhonMax — lobe 곡률반경: rho_n >= NMAX_RHO(-100)
    if rho_n < NMAX_RHO
        c(12) = NMAX_RHO - rho_n;
    end

    % c(13): GearPenaltyRhonMin — lobe 곡률반경: rho_n <= NMIN_RHO(-2)
    % 즉 |rho_n| >= 2 (lobe 곡률반경 절댓값이 2mm 이상)
    if rho_n > NMIN_RHO
        c(13) = rho_n - NMIN_RHO;
    end

    %% 접촉응력 한계 — MAX_STRESS로 정규화
    c(14) = (sig_pos - MAX_STRESS) / MAX_STRESS;
    c(15) = (sig_neg - MAX_STRESS) / MAX_STRESS;
    c(16) = (sig_in  - MAX_STRESS) / MAX_STRESS;

    % c(17): Rh < r (홀이 피치원 안에 있어야 함)
    c(17) = Rh - (r - q);   % Rh <= r - q (슬리브 겹침 방지)

    ceq = [];
end
