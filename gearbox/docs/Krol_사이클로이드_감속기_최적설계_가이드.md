---
title: 사이클로이드 감속기 최적설계 — Król 논문 재현 가이드
tags:
  - 최적설계
  - 사이클로이드감속기
  - MATLAB
  - LevenbergMarquardt
  - 페널티함수법
created: 2026-06-05
참고논문: "Król et al. (2019), Optimization of the one stage cycloidal gearbox as a non-linear least squares problem"
---

# 사이클로이드 감속기 최적설계 — Król 논문 재현 가이드

> [!info] 이 문서의 목적
> Król et al.(2019)의 사이클로이드 감속기 단일목적(체적 최소화) 최적설계를 **정식화부터 MATLAB 구현·검증까지** 그대로 따라할 수 있도록 정리한 가이드. Zenodo 공개 코드(`DOI:10.5281/zenodo.2166718`)의 구조를 기반으로 함.
>
> 다목적(체적+효율) 확장이나 `fmincon` 버전은 이 문서 범위 밖이며, 추후 별도로 진행.

---

## 0. 큰 그림 — 무엇을 하는가

Król 논문의 핵심 아이디어는 다음과 같다.

사이클로이드 감속기의 **체적을 최소화**하되, 체적을 줄이면 접촉응력이 커지므로 **접촉응력 한계와 기하학적 제약을 페널티 함수로 묶어서** 하나의 비선형 최소제곱(non-linear least squares) 문제로 만든다. 그리고 이 문제를 **Levenberg-Marquardt(LM)** 와 **Steepest Descent(SD)** 두 경도 기반 알고리즘으로 푼다.

```
설계변수 x (7개)
   │
   ▼
잔차 벡터 r(x) = [체적, 접촉응력 3개, 페널티 17개]  ← 총 21개 성분
   │
   ▼
목적: minimize  f(x) = 0.5 · ‖r(x)‖²
   │
   ▼
LM 또는 SD 알고리즘으로 국소 최소점 탐색
   │
   ▼
최적 설계변수 x* + 기어 형상 시각화
```

> [!note] 왜 "최소제곱" 형태인가
> 일반적인 제약 최적화는 `min f(x) s.t. g(x)≤0` 형태지만, Król은 모든 항(목적함수+제약)을 **잔차 벡터의 성분**으로 넣고 그 제곱합을 최소화한다. 이렇게 하면 LM 알고리즘(원래 최소제곱 전용)을 그대로 쓸 수 있다. 제약조건은 **페널티 함수**로 잔차 벡터 안에 포함된다.

---

## 1. 정식화 (Problem Formulation)

### 1.1 설계변수 (Design Variables)

$$\mathbf{x} = [\,e,\ q,\ m,\ h,\ R_S,\ R_h,\ R_W\,]^T$$

| 기호 | 코드 변수 | 의미 | 단위 | 비고 |
|---|---|---|---|---|
| $e$ | `x(1)` | 편심률 (eccentricity) | mm | 입력축 편심 거리 |
| $q$ | `x(2)` | 외부 슬리브(롤러핀) 반경 | mm | external sleeve radius |
| $m$ | `x(3)` | 단폭계수 (short-width coefficient) | – | 치형 형상 결정 |
| $h$ | `x(4)` | 기어 두께 | mm | 모든 부품 동일 두께 가정 |
| $R_S$ | `x(5)` | 내부 슬리브 반경 | mm | internal sleeve radius |
| $R_h$ | `x(6)` | 내부 슬리브 홀 반경 | mm | 사이클로이드 기어의 구멍 |
| $R_W$ | `x(7)` | 내부 슬리브 홀 위치 반경 | mm | 구멍 배치 원 반경 |

> [!important] 고정 상수 (설계변수 아님)
> 다음 값들은 최적화 대상이 아니라 **고정 파라미터**다. 코드상 `GearObjective.m`과 `GearJac.m` 안에 하드코딩되어 있다.
>
> | 기호 | 코드 | 값 | 의미 |
> |---|---|---|---|
> | $z_s$ | `zs` | 15 | 사이클로이드 치형(lobe) 수 |
> | $z_k$ | `zk` | 16 | 외부 슬리브 수 ($=z_s+1$) |
> | $z_i$ | `zi` | 8 | 내부 슬리브 수 |
> | $M_h$ | `Mh` | 1550 | 입력 토크 (N·mm) |
> | $\nu_1,\nu_2$ | `nu1,nu2` | 0.3 | 푸아송 비 (기어/슬리브) |
> | $E_1,E_2$ | `Emod1,Emod2` | 200000 | 영률 (MPa) |
> | $N$ | `N` | 1000 | 체적 수치적분 분할 수 |

### 1.2 사이클로이드 치형 매개변수 방정식

설계변수로부터 치형 곡선을 그리는 핵심 식. 피치원 반경은 $r = \dfrac{e(z_s+1)}{m}$ 로 정의된다.

$$u(\alpha) = \frac{e(z_s+1)}{m}\cos\alpha + e\cos(z_k\alpha) - q\cos\!\left(\alpha + \arctan\frac{\sin(z_s\alpha)}{\frac{1}{m}+\cos(z_s\alpha)}\right)$$

$$v(\alpha) = \frac{e(z_s+1)}{m}\sin\alpha + e\sin(z_k\alpha) - q\sin\!\left(\alpha + \arctan\frac{\sin(z_s\alpha)}{\frac{1}{m}+\cos(z_s\alpha)}\right)$$

- $\alpha$ : 매개변수(각도), $0 \le \alpha \le 2\pi$
- $u, v$ : 직교좌표계에서 치형 곡선 위 점의 좌표

> [!tip] 코드 매핑
> 이 두 식은 `GearDraw.m`(형상 그리기)과 `GearVolume.m`(체적 적분) 안에 직접 구현되어 있다. 논문 식 (3), (4)에 해당.

### 1.3 목적함수 (Objective Function)

**진짜 최소화 대상은 체적** 하나다.

$$\min_{\mathbf{x}} \quad V_{GEAR}(\mathbf{x})$$

체적은 치형 곡선의 매개변수 방정식으로부터 수치적분해서 구한다.

$$V_{WHEEL} = h \int_0^{2\pi} |v(\alpha)| \cdot \left|\frac{du(\alpha)}{d\alpha}\right| d\alpha$$

$$V_{GEAR} = V_{WHEEL} + z_k V_{EX\_SLEEVE} + z_i V_{IN\_SLEEVE} - z_i V_{IN\_HOLE}$$

수치적분은 중점법(midpoint rule)으로 구현된다($dt = 2\pi/N$, $\alpha_i = i\cdot dt + dt/2$).

> [!note] 코드상 체적 식 (`GearVolume.m`)
> 마지막 줄을 보면:
> ```matlab
> vol = vol*h - pi*Rh*Rh*zi*h + pi*q*q*zk*h + pi*Rs*Rs*zi*h;
> ```
> - `vol*h` : 사이클로이드 휠 본체 ($V_{WHEEL}$)
> - `- pi*Rh²*zi*h` : 내부 슬리브 홀 부피 빼기
> - `+ pi*q²*zk*h` : 외부 슬리브 부피 더하기
> - `+ pi*Rs²*zi*h` : 내부 슬리브 부피 더하기

### 1.4 잔차 벡터 (Residual Vector) — 실제 최소화되는 것

LM 알고리즘은 다음 21×1 잔차 벡터의 제곱합 $0.5\|r(\mathbf{x})\|^2$ 를 최소화한다. (논문 식 1, 코드 `GearObjective.m`)

$$r(\mathbf{x}) = \big[\, \underbrace{V}_{1},\ \underbrace{\sigma_{EX}^{+},\ \sigma_{EX}^{-},\ \sigma_{IN}}_{2\text{–}4},\ \underbrace{P_{V},\ P_{E},\ P_{M}^{min},\ P_{M}^{max},\ P_{H},\ P_{Rs},\ P_{Rh},\ P_{q},\ P_{Rw}^{min},\ P_{Rw}^{max}}_{5\text{–}14},\ \underbrace{P_{\rho P}^{min},\ P_{\rho P}^{max},\ P_{\rho N}^{min},\ P_{\rho N}^{max}}_{15\text{–}18},\ \underbrace{P_{\sigma EX}^{+},\ P_{\sigma EX}^{-},\ P_{\sigma IN}}_{19\text{–}21} \,\big]^T$$

| # | 성분 | 의미 |
|---|---|---|
| 1 | $V$ | 체적 (진짜 목적함수) |
| 2 | $\sigma_{EX}^{+}$ | 외부 슬리브–lobe 접촉응력 (볼록부) |
| 3 | $\sigma_{EX}^{-}$ | 외부 슬리브–pit 접촉응력 (오목부) |
| 4 | $\sigma_{IN}$ | 내부 슬리브–기어 접촉응력 |
| 5–21 | $P_{\cdots}$ | **페널티 함수들** (제약조건, 아래 1.5절) |

> [!warning] 핵심 이해 포인트
> 응력 성분(2~4)도 잔차에 들어가 있으므로 LM은 "체적도 줄이고 응력도 줄이는" 방향으로 움직인다. 동시에 페널티 항(5~21)이 제약 위반 시 큰 값을 만들어 해를 가용 영역 안으로 끌어들인다.

### 1.5 제약조건 (Constraints) — 페널티 함수법

모든 제약은 **2차 외부 페널티(quadratic exterior penalty)** 형태로 잔차 벡터에 들어간다. 일반형:

$$P(\mathbf{x}) = \begin{cases} (\text{위반량})^2 & \text{제약 위반 시} \\ 0 & \text{제약 만족 시} \end{cases}$$

> [!note] Król의 영리한 트릭 — 체적 차이로 페널티 정의
> 박스 제약(예: $m \ge m_{min}$)을 단순히 $(m_{min}-m)^2$ 로 쓰지 않고, **그 위반이 체적에 미치는 영향**으로 환산한다. 예: `GearPenaltyMmin.m`은
> $$P = \big(V(m{=}m_{min}) - V(m{=}m_{\text{현재}})\big)^2$$
> 이렇게 하면 모든 페널티 항의 단위(scale)가 체적과 같아져서 잔차 벡터의 성분들이 균형을 이룬다. (코드에 주석 처리된 `%(M_MIN-m)^2`이 원래 단순형, 실제 사용은 체적 환산형)

#### 제약조건 목록

| 페널티                    | 코드 파일                    | 제약 내용                                    | 한계값                            |
| ---------------------- | ------------------------ | ---------------------------------------- | ------------------------------ |
| $P_V$                  | `GearPenaltyVmin`        | 체적 하한 $V \ge V_{min}$                    | `MIN_VOL=1000`                 |
| $P_E$                  | `GearPenaltyEmin`        | 편심률 하한 (언더컷 방지)                          | $E_{MIN}$ 수식 계산*               |
| $P_M^{min}$            | `GearPenaltyMmin`        | 단폭계수 하한 $m \ge 0.5$                      | `M_MIN=0.5`                    |
| $P_M^{max}$            | `GearPenaltyMmax`        | 단폭계수 상한 $m \le 0.85$                     | `M_MAX=0.85`                   |
| $P_H$                  | `GearPenaltyHmin`        | 두께 하한 $h \ge 0.2$                        | `H_MIN=0.2`                    |
| $P_{Rs}$               | `GearPenaltyRsmin`       | 내부 슬리브 반경 하한 $R_S \ge 3$                 | `RS_MIN=3`                     |
| $P_{Rh}$               | `GearPenaltyRhmin`       | 홀 반경 $\ge$ 슬리브 반경 $R_h \ge R_S$          | –                              |
| $P_q$                  | `GearPenaltyqmin`        | 외부 슬리브 반경 $q \le \|\rho_{min}\|$ (간섭 방지) | –                              |
| $P_{Rw}^{min}$         | `GearPenaltyRwmin`       | $R_W \ge 2R_S$                           | –                              |
| $P_{Rw}^{max}$         | `GearPenaltyRwmax`       | $R_W \le r - R_S$                        | –                              |
| $P_{\rho P}^{min,max}$ | `GearPenaltyRhopMin/Max` | 볼록부(lobe) 곡률반경 범위                        | `PMIN_RHO=9`, `PMAX_RHO=100`   |
| $P_{\rho N}^{min,max}$ | `GearPenaltyRhonMin/Max` | 오목부(pit) 곡률반경 범위                         | `NMIN_RHO=-2`, `NMAX_RHO=-100` |
| $P_{\sigma EX}^{+}$    | `GearPenaltyStressExPos` | lobe 접촉응력 $\le 400$ MPa                  | `MAX_CONTACT_STRESS=400`       |
| $P_{\sigma EX}^{-}$    | `GearPenaltyStressExNeg` | pit 접촉응력 $\le 400$ MPa                   | `400`                          |
| $P_{\sigma IN}$        | `GearPenaltyStressIn`    | 내부 슬리브 접촉응력 $\le 400$ MPa                | `400`                          |

\* $E_{MIN} = \dfrac{q(z_k+1)}{3\sqrt{3}\,z_k}\sqrt{\dfrac{z_k+1}{z_k-1}}\sqrt{\dfrac{m^2}{1-m^2}}$ (`GearPenaltyEmin.m`)

> [!important] 한계값은 모두 `GearGetConstant.m`에 모여 있다
> ```matlab
> case 1: 400   % MAX_CONTACT_STRESS (최대 접촉응력)
> case 2: 9     % PMIN_RHO (볼록 곡률반경 하한)
> case 3: 100   % PMAX_RHO (볼록 곡률반경 상한)
> case 4: -2    % NMIN_RHO (오목 곡률반경 하한)
> case 5: -100  % NMAX_RHO (오목 곡률반경 상한)
> case 6: 1000  % MIN_VOL (체적 하한)
> ```
> 다른 감속기 사양을 시도하려면 **여기 값만 바꾸면 된다.**

### 1.6 접촉응력 수식 (Hertz 이론)

논문 식 (20)~(25). 평면 접촉(2D Hertz) 기준.

$$\sigma_{EX}^{+} = 0.5642\sqrt{\frac{(F_{EX}/h)(\rho+q)}{\eta\,\rho\,q}} \qquad \text{(lobe, 볼록-볼록)}$$

$$\sigma_{EX}^{-} = 0.5642\sqrt{\frac{(F_{EX}/h)(|\rho|-q)}{\eta\,|\rho|\,q}} \qquad \text{(pit, 볼록-오목)}$$

$$\sigma_{IN} = 0.5642\sqrt{\frac{(F_{IN}/h)(R_h-R_S)}{\eta\,R_h\,R_S}} \qquad \text{(내부 슬리브)}$$

여기서 힘과 탄성계수는:

$$F_{EX} = \frac{4 M_c}{e\,z_s\,z_k}, \quad F_{IN} = \frac{2 M_c \cos\alpha}{\pi R_W}, \quad \eta = \frac{1-\nu_1^2}{E_1} + \frac{1-\nu_2^2}{E_2}$$

$M_c = M_h \cdot z_s$ (출력 토크 = 입력 토크 × 전달비). 코드상 `GearStressExPos/Neg.m`, `GearStressIn.m`.

### 1.7 곡률반경 $\rho(\alpha)$ 와 극값 탐색

접촉응력 계산에는 접촉점의 곡률반경 $\rho$가 필요하다.

$$\rho(\alpha) = \frac{\left(u'(\alpha)^2 + v'(\alpha)^2\right)^{3/2}}{|v''(\alpha)u'(\alpha) - u''(\alpha)v'(\alpha)|}$$

- 볼록부(convex): $\rho < 0$, 오목부(concave): $\rho > 0$ (코드 부호 규약)
- **극값 위치**(최대/최소 곡률 지점)는 뉴턴법으로 찾는다 → `GearAlpha.m`
  $$\alpha_{next} = \alpha_{prev} - \frac{\rho'(\alpha)}{\rho''(\alpha)}$$
- $\rho$ 자체는 `GearRho.m`(거대한 해석식, MATLAB 심볼릭으로 유도됨), 미분은 `GearRhoPrim`(1차), `GearRhoBis`(2차)에서 수치미분.

---

## 2. 최적화 알고리즘

### 2.1 Levenberg-Marquardt (LM)

$$\mathbf{x}_{i+1} = \mathbf{x}_i - \left(H(\mathbf{x}_i) + \lambda\cdot\text{diag}[H(\mathbf{x}_i)]\right)^{-1}\nabla r(\mathbf{x}_i)$$

- $H(\mathbf{x}) = J^T(\mathbf{x})J(\mathbf{x})$ : 근사 헤시안 (Gauss-Newton)
- $\nabla r(\mathbf{x}) = J^T(\mathbf{x})\,r(\mathbf{x})$ : 경도
- $\lambda$ : 감쇠계수. 개선되면 $\lambda \leftarrow \lambda/k$(Gauss-Newton에 가깝게), 악화되면 $\lambda \leftarrow \lambda\cdot k$(Steepest Descent에 가깝게). 코드상 `k=10`, 초기 `lambda=1000`.

코드: `GearSimpleLevMar.m`

### 2.2 Steepest Descent (SD)

$$\mathbf{x}_{i+1} = \mathbf{x}_i - \alpha\cdot\nabla r(\mathbf{x}_i)$$

- $\alpha$ : 고정 스텝 크기. 코드상 `lambda=1e-13` (매우 작음 → 매우 느림, 수렴은 보장)
- 코드: `GearSteepestDescent.m`

> [!tip] 두 방법 비교 (논문 결론)
> LM은 훨씬 빠르지만(수십~수백 iteration) 불리한 조건에서 멈출 수 있고, SD는 느리지만(수천~수만 iteration) 수렴이 안정적. 논문 Table 1에서 SD는 8000+ iteration, LM은 20~130 iteration.

---

## 3. 코드 파일 역할 정리

### 3.1 메인 실행

| 파일 | 역할 |
|---|---|
| **`GearRun.m`** | **메인 스크립트.** 시작점 설정 → LM 실행 → SD 실행 → 형상 그리기 → 응력·체적 비교 |

### 3.2 핵심 함수 (Core)

| 파일 | 역할 |
|---|---|
| `GearObjective.m` | 잔차 벡터 $r(\mathbf{x})$ (21×1) 반환 |
| `GearJac.m` | 야코비 행렬 $J(\mathbf{x})$ (21×7) 반환 (수치미분 조합) |
| `GearVolume.m` | 체적 수치적분 ($V_{GEAR}$) |
| `GearRho.m` | 곡률반경 $\rho(\alpha)$ 해석식 |
| `GearAlpha.m` | $\rho$ 극값 각도 뉴턴법 탐색 |

### 3.3 최적화 엔진

| 파일 | 역할 |
|---|---|
| `GearSimpleLevMar.m` | Levenberg-Marquardt 반복 |
| `GearSteepestDescent.m` | Steepest Descent 반복 |

### 3.4 응력 계산 + 미분

| 파일 | 역할 |
|---|---|
| `GearStressExPos/Neg.m` | 외부 슬리브 접촉응력 (lobe/pit) |
| `GearStressIn.m` | 내부 슬리브 접촉응력 |
| `GearStressEx d{e,q,m,h} Pos/Neg.m` | 외부 응력의 각 변수 편미분 (해석식) |
| `GearStressInd{h,Rs,Rh,Rw}.m` | 내부 응력의 각 변수 편미분 (해석식) |

### 3.5 체적 미분 (수치미분)

| 파일 | 역할 |
|---|---|
| `GearVolumed{e,q,m,h,Rs,Rh}.m` | 체적의 각 변수 편미분 (전진차분, $\Delta=0.001$) |

### 3.6 페널티 함수 (제약조건) — 각각 값/야코비 쌍

| 파일 | 역할 |
|---|---|
| `GearPenalty{Vmin,Emin,Mmin,Mmax,Hmin,Rsmin,Rhmin,qmin}.m` | 박스 제약 페널티 |
| `GearPenalty{RwminRwmax}.m` | $R_W$ 범위 페널티 |
| `GearPenaltyRho{p,n}{Min,Max}.m` | 곡률반경 범위 페널티 |
| `GearPenaltyStress{ExPos,ExNeg,In}.m` | 접촉응력 한계 페널티 |
| `...Jac.m` | 위 각 페널티의 야코비 성분 |
| `GearGetConstant.m` | 모든 한계 상수 중앙 관리 |

### 3.7 보조/검증용

| 파일 | 역할 |
|---|---|
| `GearDraw.m` | 최적화 전(빨강)/후(파랑) 기어 형상 plot |
| `GearRosenbrock.m`, `...Jac.m`, `GearTestRosenbrock.m` | LM 알고리즘 검증용 Rosenbrock 테스트 |
| `GearNewton.m` | (보조) 뉴턴법 |

> [!warning] 미사용/주의 파일
> `GearStressExdb*`, `GearVolumedb`, `GearPenaltybmin*`, `GearRhoBis`의 일부는 잔차/야코비에 직접 안 쓰이거나 보조용이다. 또 `GearLevMar`(주석 처리됨)는 코드에 없고 `GearSimpleLevMar`만 실제 사용된다. `GearNewton.m`은 `GearSymbolicsDiff`를 호출하는데 이 파일이 digest에 없으므로 **그대로 실행하면 에러** → `GearNewton`은 호출하지 말 것(메인 흐름에선 안 쓰임).

---

## 4. MATLAB에서 따라하기 (Step-by-Step)

### Step 1 — 파일 준비

> [!todo] 폴더 세팅
> 1. Zenodo(`DOI:10.5281/zenodo.2166718`)에서 zip을 받아 `gearbox/` 폴더의 `.m` 파일 전부를 한 폴더에 둔다.
> 2. MATLAB에서 그 폴더를 **현재 폴더(Current Folder)** 로 설정하거나 `addpath`로 등록.
> 3. (선택) 함수명·파일명이 정확히 일치하는지 확인. MATLAB은 파일명=함수명이어야 한다.

### Step 2 — LM 알고리즘 자체 검증 (Rosenbrock)

본 문제에 들어가기 전, LM 구현이 맞는지 표준 테스트 함수로 확인.

```matlab
% Rosenbrock 함수의 최소점은 (1,1)
GearTestRosenbrock
% 결과: xopt ≈ [1; 1], GearValueNext ≈ 0 이면 LM 정상 동작
```

> [!check] 기대 결과
> `xopt`이 `[1, 1]`에 가깝고 `GearValueNext`가 거의 0이면 LM 엔진이 정상.

### Step 3 — 시작점 설정 & 메인 실행

`GearRun.m`의 첫 줄이 시작점이다.

```matlab
xp = [3.0, 3.3, 0.56, 10, 5, 8, 40]';   % [e, q, m, h, Rs, Rh, Rw]
```

그대로 실행:

```matlab
GearRun
```

이 스크립트가 자동으로:
1. LM으로 최적화 → `xfin` (그림 1: 빨강=전, 파랑=후)
2. SD로 최적화 → `xfin2` (그림 2)
3. 최적화 전후 체적·접촉응력을 워크스페이스 변수로 저장

### Step 4 — 결과 확인

워크스페이스에서 다음 변수들을 비교한다.

| 변수 | 의미 |
|---|---|
| `VolumePrev` / `VolumeNext` | LM 최적화 전 / 후 체적 |
| `VolumeFin2` | SD 최적화 후 체적 |
| `StressExPosNext`, `StressExNegNext`, `StressInNext` | LM 후 세 접촉응력 |
| `xfin`, `xfin2` | LM / SD 최적 설계변수 |

> [!check] 논문 Table 1 재현 확인
> 시작점 `[2.8, 3.5, 0.7, 10, 5, 8, 30]`(논문 variant 1)으로 바꾸면 체적이 약 `1.17e5 → 6.07e4`(LM)로 줄어드는 결과를 재현할 수 있다. 응력은 400 MPa 한계 근처로 올라간다(체적↓ → 응력↑ 트레이드오프).

### Step 5 — 여러 시작점 실험 (논문 Table 2)

```matlab
% 논문의 4가지 시작점
starts = [2.8, 3.5, 0.7,  10, 5, 8, 30;
          1.5, 2.0, 0.65,  8, 3, 5, 20;
          2.9, 3.4, 0.7,  10, 5, 8, 30;
          3.0, 3.3, 0.7,  10, 5, 8, 30];

for i = 1:size(starts,1)
    xp = starts(i,:)';
    xfin = GearSimpleLevMar(3000, @GearObjective, @GearJac, xp);
    V0 = GearVolume(xp(1)*16/xp(3),  xp(1),16,15, xp(2),xp(3),1000, xp(4),xp(5),xp(6),8);
    V1 = GearVolume(xfin(1)*16/xfin(3), xfin(1),16,15, xfin(2),xfin(3),1000, xfin(4),xfin(5),xfin(6),8);
    fprintf('Variant %d:  V0=%.0f -> V1=%.0f  (%.1f%% 감소)\n', i, V0, V1, 100*(V0-V1)/V0);
end
```

> [!warning] 시작점이 가용 영역 밖이면
> 논문 variant 2처럼 시작점에서 이미 응력 한계(400 MPa)를 크게 초과하면, LM이 페널티 영역 안에서 멈춰 **체적이 오히려 늘어날 수 있다**. 이는 경도 기반 방법이 국소 최소점에 갇히는 현상으로, 정상적인 결과다(논문에서도 그렇게 설명).

---

## 5. 검증 체크리스트

> [!check] 구현이 제대로 됐는지 확인할 것들
> - [ ] `GearTestRosenbrock` → `[1,1]` 수렴 (LM 엔진 OK)
> - [ ] `GearVolume`에 시작점 넣어 양수의 합리적 체적값 나오는지
> - [ ] `GearRho`가 lobe에서 음수, pit에서 양수 부호 나오는지
> - [ ] `GearObjective(xp)`가 21×1 벡터 반환하는지 (`size` 확인)
> - [ ] `GearJac(xp)`가 21×7 행렬 반환하는지 (`size` 확인)
> - [ ] LM 결과 `xfin`이 모든 박스 제약 안에 있는지 (특히 $0.5\le m\le0.85$, $R_h\ge R_S$)
> - [ ] 최적화 후 세 접촉응력이 400 MPa 이하인지
> - [ ] `GearDraw` 그림에서 파란 곡선(최적화 후)이 빨간 곡선(전)보다 작은지

---

## 6. 자주 막히는 부분 (Troubleshooting)

> [!bug] 흔한 오류
> - **`GearSymbolicsDiff` 관련 에러**: `GearNewton.m`이 호출하지만 파일이 없다. 메인 흐름은 `GearAlpha`(자체 뉴턴법)를 쓰므로 `GearNewton`은 호출하지 말 것.
> - **`Undefined function` 에러**: `.m` 파일 일부가 폴더에 없거나 파일명 오타. digest의 76개 파일이 모두 있는지 확인.
> - **수렴 안 함 / NaN**: 시작점이 너무 극단적이거나 $m \to 1$ 근처($E_{MIN}$ 식에 $\sqrt{m^2/(1-m^2)}$ 발산). $m$을 0.5~0.85 사이로.
> - **응력이 비현실적으로 큼**: 곡률반경 $\rho$가 0 근처면 응력 발산. 곡률 페널티 한계(`PMIN_RHO` 등)가 막아주지만 시작점이 나쁘면 초기 발산 가능.
> - **`diag(diag(...))` 관련**: LM의 `jacTjac+lambda*diag(diag(jacTjac))`에서 $J^TJ$가 특이(singular)하면 역행렬 불안정. `lambda` 초기값을 키워볼 것.

---

## 7. 다음 단계 (이 문서 범위 밖, 추후 진행)

> [!note] 확장 아이디어 — 나중에
> 1. **`fmincon` 버전**: 잔차 벡터의 1번 성분(체적)만 목적함수로, 5~21번 페널티를 명시적 부등호 제약 `g(x)≤0`으로 분리. 수업의 SQP·KKT 해석과 직접 연결.
> 2. **다목적(체적+효율)**: Wang et al.(2016)의 효율 수식을 목적함수에 추가 → 가중합법으로 `fmincon` 반복 → Pareto front. (교재 18.4절 가중합법)
> 3. **라그랑지 승수 해석**: `fmincon` 출력 `lambda` 구조체에서 활성 제약의 승수를 뽑아 어떤 제약이 체적에 가장 큰 영향을 주는지 분석 (교재 후최적성 해석).

---

## 부록 A — 기호 빠른 참조

| 기호 | 의미 | 기호 | 의미 |
|---|---|---|---|
| $e$ | 편심률 | $\rho$ | 곡률반경 |
| $q$ | 외부 슬리브 반경 | $\eta$ | 탄성계수 |
| $m$ | 단폭계수 | $F_{EX}$ | 외부 슬리브 힘 |
| $h$ | 기어 두께 | $F_{IN}$ | 내부 슬리브 힘 |
| $R_S$ | 내부 슬리브 반경 | $M_h$ | 입력 토크 |
| $R_h$ | 홀 반경 | $M_c$ | 출력 토크 ($=M_h z_s$) |
| $R_W$ | 홀 위치 반경 | $r$ | 피치원 반경 ($=e(z_s{+}1)/m$) |
| $z_s$ | 치형 수 (15) | $z_k$ | 외부 슬리브 수 (16) |
| $z_i$ | 내부 슬리브 수 (8) | $\alpha$ | 매개변수 각도 |

## 부록 B — 참고 자료

- **주 논문**: Król, R., Wikło, M., Olejarczyk, K., Kołodziejczyk, K., Zieja, A. (2019). *Optimization of the one stage cycloidal gearbox as a non-linear least squares problem*. Advances in Mechanism and Machine Science, MMS 73, pp. 1039–1048.
- **코드**: Król, R. *Software for the one stage cycloidal gearbox optimization* (MATLAB scripts), `DOI:10.5281/zenodo.2166718`
- **교재**: Arora, J.S. *Introduction to Optimum Design*, 4th Ed. — Ch. 7(MATLAB), Ch. 11(경도법), 페널티 함수법, Ch. 18(다목적)
- **비교 논문**: Wang, J., Luo, S., Su, D. (2016). *Multi-objective optimal design of cycloid speed reducer based on genetic algorithm*. MMT 102, 135–148.
