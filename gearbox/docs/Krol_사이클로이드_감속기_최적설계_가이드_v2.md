---
title: 사이클로이드 감속기 최적설계 — Król 논문 재현 가이드 (v2)
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

## 0. 배경 — 왜 사이클로이드 감속기인가

### 0.1 하모닉 드라이브 vs 사이클로이드 드라이브

논문 서두에 "로봇처럼 작은 크기가 필요한 곳에는 하모닉 기어가 쓰이고, 일반적으로 다른 대안이 없다"는 문제 제기가 나온다. 이는 둘이 같다는 뜻이 아니라, **사이클로이드를 하모닉의 대안으로 만들겠다는 논문의 동기**이다.

| | 하모닉 드라이브 | 사이클로이드 드라이브 |
|---|---|---|
| 핵심 원리 | 유연한 금속 컵(flexspline)이 탄성 변형하며 맞물림 | 사이클로이드 디스크가 핀과 편심 회전하며 맞물림 |
| 크기 | 매우 작고 얇음 | 상대적으로 두껍고 무거움 |
| 충격/과부하 | 약함 (flexspline 피로 파손) | 강함 (다수 치형이 동시에 하중 분담) |
| 강성·수명 | 상대적으로 낮음 | 높음 |

사이클로이드의 체적을 최적화로 충분히 줄일 수 있다면, 하모닉의 단점(충격 취약)을 보완하면서도 소형화를 만족하는 진짜 대안이 된다. **이것이 체적 최소화 최적설계의 근본 동기다.**

### 0.2 싱글 디스크 vs 듀얼 디스크

**Król 논문은 싱글 디스크(1개) 모델**이다. 논문에 "single cycloidal wheel with 15 lobes"로 명시되어 있다.

Wang et al., Qi et al. 등 다른 논문들은 **듀얼 디스크(2개, 180° 위상차)** 를 기본으로 한다. 두 디스크가 편심 회전의 진동을 서로 상쇄하여 동적 균형을 맞추기 때문이다.

최적설계 수식 관점에서 둘의 차이:

| | 싱글 디스크 | 듀얼 디스크 |
|---|---|---|
| 체적 | 디스크 1개 | 디스크 2개 (약 2배) |
| 하중 분담 | 1개가 전체 토크 부담 | 2개가 절반씩 분담 → 접촉응력 감소 |
| 동적 균형 | 진동 발생 | 상쇄되어 안정 |

본 프로젝트는 Król의 싱글 디스크 모델을 그대로 따른다. 보고서에 "단일 디스크 모델 기준"임을 명시할 것.

### 0.3 큰 그림 — 무엇을 하는가

사이클로이드 감속기의 **체적을 최소화**하되, 체적을 줄이면 접촉응력이 커지므로 **접촉응력 한계와 기하학적 제약을 페널티 함수로 묶어서** 하나의 비선형 최소제곱(non-linear least squares) 문제로 만든다. 이 문제를 **Levenberg-Marquardt(LM)** 와 **Steepest Descent(SD)** 두 경도 기반 알고리즘으로 푼다.

```
설계변수 x (7개)
   │
   ▼
잔차 벡터 r(x) (21개 성분)
 = [체적(1), 접촉응력(3), 페널티(17)]
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
> 일반적인 제약 최적화는 `min f(x) s.t. g(x)≤0` 형태지만, Król은 최소화하고 싶은 것들(목적함수+제약)을 **잔차 벡터의 성분**으로 모두 넣고 그 제곱합을 최소화한다. 이렇게 하면 LM 알고리즘(원래 데이터 피팅용 최소제곱 전용)을 수정 없이 그대로 쓸 수 있다.

---

## 1. 정식화 (Problem Formulation)

### 1.1 설계변수 (Design Variables)

$$\mathbf{x} = [\,e,\ q,\ m,\ h,\ R_S,\ R_h,\ R_W\,]^T$$

| 기호 | 코드 변수 | 의미 | 단위 | 비고 |
|---|---|---|---|---|
| $e$ | `x(1)` | 편심률 (eccentricity) | mm | 입력축-디스크 중심 간 어긋난 거리 |
| $q$ | `x(2)` | 외부 슬리브 반경 | mm | external sleeve radius |
| $m$ | `x(3)` | 단폭계수 (short-width coefficient) | – | 치형 형상 결정, 범위 0.5~0.85 |
| $h$ | `x(4)` | 기어 두께 | mm | 모든 부품 동일 두께 가정 |
| $R_S$ | `x(5)` | 내부 슬리브 반경 | mm | internal sleeve radius |
| $R_h$ | `x(6)` | 내부 슬리브 홀 반경 | mm | 사이클로이드 기어의 구멍 반경 |
| $R_W$ | `x(7)` | 내부 슬리브 홀 위치 반경 | mm | 구멍 배치 원 반경 |

> [!note] 외부 슬리브 위치 반경은 설계변수가 아닌 이유
> 외부 슬리브 16개는 **피치원 반경** $r = e(z_s+1)/m$ 위에 균등 배치된다. 즉 피치원 반경 = 외부 슬리브 위치 반경이며, $e$와 $m$을 최적화하면 자동으로 결정된다. 별도 설계변수가 필요 없다. 반면 내부 슬리브는 출력축과 연결되는 부분으로 치형과 독립적으로 크기·위치를 설계자가 선택할 수 있어 설계변수로 들어간다.

> [!note] 논문마다 설계변수가 다른 이유
> - **Król**: 치형 기하학 파라미터 → "어떤 모양으로 만들까"
> - **Wang et al.**: 감속기 전체 치수(핀 직경, 기어 폭 등) → "이미 있는 감속기를 어떤 치수로 만들까"
>
> 두 논문의 목적 자체가 다르기 때문. 따라서 결과를 숫자로 직접 비교하는 것보다 방법론(경도법 vs GA, 단목적 vs 다목적)을 비교하는 것이 적절하다.

> [!important] 고정 상수 (설계변수 아님)
> 다음 값들은 최적화 대상이 아닌 **고정 파라미터**다. 코드상 `GearObjective.m`과 `GearJac.m` 안에 하드코딩되어 있다.
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

#### 매개변수 $\alpha$ 란

$\alpha$는 치형 곡선의 **매개변수(parameter)** 로, 물리적으로는 입력축이 회전하는 각도에 해당한다. $\alpha$가 $0$에서 $2\pi$까지 변하면 $(u(\alpha), v(\alpha))$가 치형 곡선을 **정확히 한 바퀴** 그린다.

```
α = 0        α = π/2      α = π       α = 3π/2     α = 2π
   ●              ●            ●            ●         ●(시작점 복귀)
(시작점)      (1/4 바퀴)   (반 바퀴)   (3/4 바퀴)
```

#### 치형 방정식

기어 맞물림 원리(envelope theory)를 적용해 유도된 **표준 확립 공식**이다. 사이클로이드 감속기 설계에서 공통으로 쓰인다. 피치원 반경은 $r = \dfrac{e(z_s+1)}{m}$ 로 정의된다.

$$u(\alpha) = \frac{e(z_s+1)}{m}\cos\alpha + e\cos(z_k\alpha) - q\cos\!\left(\alpha + \arctan\frac{\sin(z_s\alpha)}{\frac{1}{m}+\cos(z_s\alpha)}\right)$$

$$v(\alpha) = \frac{e(z_s+1)}{m}\sin\alpha + e\sin(z_k\alpha) - q\sin\!\left(\alpha + \arctan\frac{\sin(z_s\alpha)}{\frac{1}{m}+\cos(z_s\alpha)}\right)$$

> [!tip] 코드 매핑
> 이 두 식은 `GearDraw.m`(형상 그리기)과 `GearVolume.m`(체적 적분) 안에 직접 구현되어 있다. 논문 식 (3), (4)에 해당.

### 1.3 목적함수 (Objective Function)

**진짜 최소화 대상은 체적** 하나다.

$$\min_{\mathbf{x}} \quad V_{GEAR}(\mathbf{x})$$

#### 체적 계산 — Green 정리 적용

체적은 치형 곡선의 매개변수 방정식으로부터 수치적분으로 구한다. 핵심 아이디어는 **Green 정리**다.

Green 정리의 일반형:
$$\oint_C P\,dx + Q\,dy = \iint_D \left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right)dx\,dy$$

여기서 $P=0,\ Q=x$로 놓으면 $\partial Q/\partial x = 1,\ \partial P/\partial y = 0$이므로:
$$\oint_C x\,dy = \iint_D 1\,dx\,dy = A$$

즉 **내부 이중적분(면적)이 경계선 선적분으로 바뀐다.** 복잡한 사이클로이드 경계를 이중적분 없이, 경계선만 한 바퀴 돌면서 계산할 수 있다.

경계선이 $x=u(\alpha),\ y=v(\alpha)$로 표현되면 $dy = \frac{dv}{d\alpha}d\alpha$이므로:
$$A = \int_0^{2\pi} u(\alpha)\cdot\frac{dv(\alpha)}{d\alpha}\,d\alpha$$

Król 코드는 부호 처리상 $y\,dx$ 버전($P=-y, Q=0$)의 절댓값 형태를 쓴다:
$$V_{WHEEL} = h \int_0^{2\pi} |v(\alpha)| \cdot \left|\frac{du(\alpha)}{d\alpha}\right| d\alpha$$

#### 왜 수치적분인가

$u(\alpha)$의 $\arctan$ 안에 $\cos(z_s\alpha)$가 중첩된 삼각함수가 있어, $|v(\alpha)| \cdot |du/d\alpha|$의 원시함수(닫힌 형태)가 존재하지 않는다. $\int e^{-x^2}dx$처럼 **적분값은 존재하지만 수식으로 표현 불가능**한 경우다. 따라서 $\alpha$를 잘게 쪼개 숫자로 더하는 수치적분을 쓴다.

반면 슬리브와 홀은 원기둥이므로 $\pi r^2 h$로 해석적으로 바로 계산된다.

수치적분은 중점법(midpoint rule)으로 구현된다($dt = 2\pi/N$, $\alpha_i = i\cdot dt + dt/2$).

$$V_{GEAR} = V_{WHEEL} + z_k \cdot \pi q^2 h + z_i \cdot \pi R_S^2 h - z_i \cdot \pi R_h^2 h$$

> [!note] 코드상 체적 식 (`GearVolume.m`)
> 마지막 줄:
> ```matlab
> vol = vol*h - pi*Rh*Rh*zi*h + pi*q*q*zk*h + pi*Rs*Rs*zi*h;
> ```
> - `vol*h` : 사이클로이드 휠 본체 수치적분 결과 × 두께
> - `- pi*Rh²*zi*h` : 내부 슬리브 홀 8개 빼기 (재료 없음)
> - `+ pi*q²*zk*h` : 외부 슬리브 16개 더하기
> - `+ pi*Rs²*zi*h` : 내부 슬리브 8개 더하기

### 1.4 잔차 벡터 (Residual Vector) — 실제 최소화되는 것

#### 잔차 벡터의 정체

LM 알고리즘은 원래 데이터 피팅용으로, $\frac{1}{2}\|r\|^2$을 최소화한다. Król은 **최소화하고 싶은 것들의 목록**을 잔차 벡터 성분으로 집어넣어 LM을 최적설계에 활용한다.

> **잔차 벡터 = "동시에 작게 만들고 싶은 것들의 목록"**

모든 성분이 0에 가까울수록 좋은 설계다. 체적도 작고, 응력도 작고, 제약도 안 위반하는 설계가 $\|r\|^2$을 최소화하는 해가 된다.

#### 잔차 벡터 구성 (21×1)

LM 알고리즘은 다음 잔차 벡터의 제곱합 $0.5\|r(\mathbf{x})\|^2$를 최소화한다. (논문 식 1, 코드 `GearObjective.m`)

$$r(\mathbf{x}) = \big[\, \underbrace{V}_{1},\ \underbrace{\sigma_{EX}^{+},\ \sigma_{EX}^{-},\ \sigma_{IN}}_{2\text{–}4},\ \underbrace{P_{V},\ \cdots,\ P_{Rw}^{max}}_{5\text{–}14},\ \underbrace{P_{\rho P}^{min},\ \cdots,\ P_{\rho N}^{max}}_{15\text{–}18},\ \underbrace{P_{\sigma EX}^{+},\ P_{\sigma EX}^{-},\ P_{\sigma IN}}_{19\text{–}21} \,\big]^T$$

| # | 성분 | 코드 | 역할 |
|---|---|---|---|
| 1 | $V$ | `GearVolume` | 체적 값 그대로 (진짜 목적) |
| 2 | $\sigma_{EX}^{+}$ | `GearStressExPos` | 외부 lobe 응력 값 그대로 |
| 3 | $\sigma_{EX}^{-}$ | `GearStressExNeg` | 외부 pit 응력 값 그대로 |
| 4 | $\sigma_{IN}$ | `GearStressIn` | 내부 슬리브 응력 값 그대로 |
| 5–14 | $P_{\cdots}$ | `GearPenalty*` | 박스 제약 + $R_W$ 범위 페널티 |
| 15–18 | $P_{\rho\cdots}$ | `GearPenaltyRho*` | 곡률반경 범위 페널티 |
| 19 | $P_{\sigma EX}^{+}$ | `GearPenaltyStressExPos` | 외부 lobe 응력 한계 페널티 |
| 20 | $P_{\sigma EX}^{-}$ | `GearPenaltyStressExNeg` | 외부 pit 응력 한계 페널티 |
| 21 | $P_{\sigma IN}$ | `GearPenaltyStressIn` | 내부 슬리브 응력 한계 페널티 |

> [!warning] 응력이 두 번 들어가는 이유
> 세 접촉응력(외부 lobe/pit, 내부 슬리브) **모두** 잔차 벡터에 **두 번씩** 들어간다.
> - **2~4번 (값 그대로)**: "응력도 같이 줄이는 방향"으로 알고리즘을 유도
> - **19~21번 (페널티)**: "400 MPa을 넘으면 강하게 벌칙" → 경계 조건 역할
>
> 둘이 협력해서 응력을 400 MPa 근처로 수렴시킨다.

### 1.5 제약조건 (Constraints) — 페널티 함수법

모든 제약은 **2차 외부 페널티(quadratic exterior penalty)** 형태로 잔차 벡터에 들어간다. 일반형:

$$P(\mathbf{x}) = \begin{cases} (\text{위반량})^2 & \text{제약 위반 시} \\ 0 & \text{제약 만족 시} \end{cases}$$

> [!note] Król의 스케일 균형 트릭
> 박스 제약(예: $m \ge m_{min}$)을 단순히 $(m_{min}-m)^2$으로 쓰면 단위가 달라 잔차 성분 간 균형이 깨진다. Król은 **위반이 체적에 미치는 영향**으로 환산한다. 예: `GearPenaltyMmin.m`은
> $$P = \big(V(m{=}m_{min}) - V(m{=}m_{\text{현재}})\big)^2$$
> 이렇게 하면 모든 페널티 항의 단위(scale)가 체적과 같아져서 잔차 벡터 성분들이 균형을 이룬다.
> (코드에 주석 처리된 `%(M_MIN-m)^2`이 원래 단순형, 실제 사용은 체적 환산형)

#### 제약조건 목록 및 한계값 근거

| 페널티 | 코드 파일 | 제약 내용 | 한계값 | 근거 |
|---|---|---|---|---|
| $P_V$ | `GearPenaltyVmin` | $V \ge 1000$ | 1000 mm³ | 수치 안정성 |
| $P_E$ | `GearPenaltyEmin` | 편심률 하한 | $E_{MIN}$ 수식* | **이론값** — 곡률반경 $\rho_{min}>0$ 조건에서 유도 |
| $P_M^{min}$ | `GearPenaltyMmin` | $m \ge 0.5$ | 0.5 | 경험치 — 맞물림 성능 하한 |
| $P_M^{max}$ | `GearPenaltyMmax` | $m \le 0.85$ | 0.85 | $m\to1$에서 $E_{MIN}$ 발산 방지 |
| $P_H$ | `GearPenaltyHmin` | $h \ge 0.2$ | 0.2 mm | 제조 가능성 + 수치 안정성 |
| $P_{Rs}$ | `GearPenaltyRsmin` | $R_S \ge 3$ | 3 mm | 경험치 — 제조·강도 하한 |
| $P_{Rh}$ | `GearPenaltyRhmin` | $R_h \ge R_S$ | – | 기하학적 필수 조건 (슬리브가 홀 안에 들어가야 함) |
| $P_q$ | `GearPenaltyqmin` | $q \le \|\rho_{min}\|$ | – | 슬리브-치형 간섭 방지 |
| $P_{Rw}^{min}$ | `GearPenaltyRwmin` | $R_W \ge 2R_S$ | – | 슬리브 간 간섭 방지 |
| $P_{Rw}^{max}$ | `GearPenaltyRwmax` | $R_W \le r-R_S$ | – | 휠 외경 초과 방지 |
| $P_{\rho P}^{min,max}$ | `GearPenaltyRhopMin/Max` | lobe 곡률반경 범위 | 9~100 mm | 경험치 — 실용적 치형 형상 |
| $P_{\rho N}^{min,max}$ | `GearPenaltyRhonMin/Max` | pit 곡률반경 범위 | −100~−2 mm | 경험치 |
| $P_{\sigma EX}^{+}$ | `GearPenaltyStressExPos` | lobe 응력 $\le 400$ MPa | 400 MPa | 강재 허용 접촉응력 (재료 기준) |
| $P_{\sigma EX}^{-}$ | `GearPenaltyStressExNeg` | pit 응력 $\le 400$ MPa | 400 MPa | 동일 |
| $P_{\sigma IN}$ | `GearPenaltyStressIn` | 내부 응력 $\le 400$ MPa | 400 MPa | 동일 |

\* 편심률 하한: $E_{MIN} = \dfrac{q(z_k+1)}{3\sqrt{3}\,z_k}\sqrt{\dfrac{z_k+1}{z_k-1}}\sqrt{\dfrac{m^2}{1-m^2}}$

> [!note] 단폭계수 범위 근거 보충
> - **하한 0.5**: $m$이 작을수록 치형이 납작해져 맞물림 불량. 실무 경험치.
> - **상한 0.85**: $m\to1$이면 $\sqrt{m^2/(1-m^2)}$ 발산 → $E_{MIN}$ 무한대 → 현실적 설계 불가. Wang et al. 권장 범위(0.65~0.9)와 비교하면 Król이 약간 보수적.

> [!note] 접촉응력 한계 400 MPa 근거
> 코드가 `Emod=200000 MPa`, `nu=0.3`으로 **일반 강재**를 가정하므로 Hertz 접촉 허용응력 400~600 MPa 범위의 보수적 값을 택한 것이다. Wang et al.이 `σHP=1000 MPa`를 쓰는 것은 GCr15 베어링강(표면 경화)을 가정하기 때문으로, **재료와 열처리에 따라 달라지는 설계 입력값**이다.

> [!important] 한계값 중앙 관리 — `GearGetConstant.m`
> ```matlab
> case 1: 400   % MAX_CONTACT_STRESS
> case 2: 9     % PMIN_RHO (볼록 곡률반경 하한)
> case 3: 100   % PMAX_RHO (볼록 곡률반경 상한)
> case 4: -2    % NMIN_RHO (오목 곡률반경 하한)
> case 5: -100  % NMAX_RHO (오목 곡률반경 상한)
> case 6: 1000  % MIN_VOL
> ```
> 다른 사양을 시도하려면 여기 값만 바꾸면 된다.

### 1.6 접촉응력 수식 (Hertz 이론)

#### 이론 출처

1882년 Heinrich Hertz가 유도한 **탄성체 접촉 역학 공식**으로, 기어 설계의 표준 공식이다. 핵심 가정은 (1) 접촉 면적이 물체 크기에 비해 매우 작고, (2) 재료가 선형 탄성을 따른다는 것이다.

슬리브와 치형은 모두 원기둥 형태라 **선 접촉(2D Hertz)** 이 적용된다:

$$\sigma_{max} = 0.5642\sqrt{\frac{F/L}{\eta \cdot \rho^*}}$$

여기서 0.5642 = $\sqrt{2/\pi}$, $F/L$은 단위 길이당 하중, $\rho^*$는 등가 곡률반경이다.

#### 등가 곡률반경 $\rho^*$ — lobe와 pit가 다른 이유

접촉 형태에 따라 $\rho^*$가 달라진다.

**볼록-볼록 접촉 (lobe):**
$$\rho^* = \frac{\rho \cdot q}{\rho + q} \quad \Rightarrow \quad \sigma_{EX}^{+} = 0.5642\sqrt{\frac{(F_{EX}/h)(\rho+q)}{\eta\,\rho\,q}}$$

**볼록-오목 접촉 (pit):** 오목한 쪽이 볼록한 쪽을 감싸는 형태
$$\rho^* = \frac{|\rho| \cdot q}{|\rho| - q} \quad \Rightarrow \quad \sigma_{EX}^{-} = 0.5642\sqrt{\frac{(F_{EX}/h)(|\rho|-q)}{\eta\,|\rho|\,q}}$$

**홀-슬리브 접촉 (내부 슬리브):** 볼록-오목 동일 원리
$$\rho^* = \frac{R_h \cdot R_S}{R_h - R_S} \quad \Rightarrow \quad \sigma_{IN} = 0.5642\sqrt{\frac{(F_{IN}/h)(R_h-R_S)}{\eta\,R_h\,R_S}}$$

> [!note] $R_h > R_S$ 제약의 물리적 의미
> $R_h \le R_S$이면 $\rho^*$ 분모가 0 이하 → 응력 발산. 슬리브가 홀보다 커서 조립 불가능한 설계이기도 하다. `GearPenaltyRhmin`이 이를 강제하는 이유.

힘과 탄성계수:
$$F_{EX} = \frac{4 M_c}{e\,z_s\,z_k}, \quad F_{IN} = \frac{2 M_c \cos\alpha}{\pi R_W}, \quad \eta = \frac{1-\nu_1^2}{E_1} + \frac{1-\nu_2^2}{E_2}, \quad M_c = M_h \cdot z_s$$

코드: `GearStressExPos/Neg.m`, `GearStressIn.m`

### 1.7 곡률반경 $\rho(\alpha)$ 와 극값 탐색

#### 왜 극값을 찾나

접촉응력 식에서 $\rho$가 작을수록 $\sigma$가 커진다. 치형 곡선을 따라가면 $\rho$가 계속 바뀌는데, **lobe 끝에서 $|\rho|$ 최소 → 응력 최대**가 된다. 설계 검증은 항상 최악의 지점 기준이므로, 곡률반경이 극소가 되는 위치 $\alpha$를 찾아야 한다.

#### 곡률반경 공식

매개변수 곡선의 일반 곡률반경 공식:

$$\rho(\alpha) = \frac{\left(u'(\alpha)^2 + v'(\alpha)^2\right)^{3/2}}{|v''(\alpha)u'(\alpha) - u''(\alpha)v'(\alpha)|}$$

- 볼록부(lobe): $\rho < 0$, 오목부(pit): $\rho > 0$ (코드 부호 규약)
- $\rho(\alpha)$ 해석식은 `GearRho.m`에 MATLAB 심볼릭으로 유도된 거대한 식으로 구현됨

#### 왜 수치적으로 극값을 찾나

극값 위치를 해석적으로 구하려면 $\rho'(\alpha)=0$을 풀어야 한다. 그러나 `GearRho.m`의 식이 복잡한 삼각함수 합성이라 $\rho'(\alpha)=0$을 기호로 풀 수 없다. 체적 적분과 같은 이유다. 그래서 **뉴턴법**으로 수치적으로 찾는다.

#### 뉴턴법 (Newton's Method)

$f(x)=0$을 만족하는 $x$를 찾는 방법이다. 여기서는 $f=\rho'$로 놓아 $\rho'(\alpha)=0$인 $\alpha$를 탐색한다.

아이디어: 현재 점에서 접선을 그어 그 접선이 0이 되는 지점을 다음 추정값으로 삼는다.

$$\alpha_{next} = \alpha_{prev} - \frac{\rho'(\alpha_{prev})}{\rho''(\alpha_{prev})}$$

```
ρ'(α)
↑
|    ×  ← 현재점 (기울기 ≠ 0)
|   /
|  /  ← 접선
| /
+----------→ α
     ↑
    다음 추정값 (접선이 0 되는 곳)
```

이를 `|prev - next| < 1e-6`이 될 때까지 반복하면 빠르게 수렴한다. 뉴턴법은 **2차 수렴(quadratic convergence)** 으로 매 반복마다 유효 자릿수가 약 2배씩 늘어난다.

코드: `GearAlpha.m`. $\rho$ 미분은 `GearRhoPrim`(1차), `GearRhoBis`(2차)에서 수치미분.

---

## 2. 최적화 알고리즘

### 2.1 Levenberg-Marquardt (LM)

LM은 **Gauss-Newton법과 Steepest Descent를 적응적으로 결합**한 방법이다.

$$\mathbf{x}_{i+1} = \mathbf{x}_i - \left(H(\mathbf{x}_i) + \lambda\cdot\text{diag}[H(\mathbf{x}_i)]\right)^{-1}\nabla r(\mathbf{x}_i)$$

- $H(\mathbf{x}) = J^T(\mathbf{x})J(\mathbf{x})$ : 근사 헤시안 (Gauss-Newton)
- $\nabla r(\mathbf{x}) = J^T(\mathbf{x})\,r(\mathbf{x})$ : 경도(gradient)
- $\lambda$ : 감쇠계수 (damping factor)
  - $\lambda \to 0$ : Gauss-Newton에 가까워짐 (빠른 수렴)
  - $\lambda \to \infty$ : Steepest Descent에 가까워짐 (안정적)
  - 개선되면 $\lambda \leftarrow \lambda/k$, 악화되면 $\lambda \leftarrow \lambda\cdot k$로 자동 조정

코드상 `k=10`, 초기 `lambda=1000`. 코드: `GearSimpleLevMar.m`

### 2.2 Steepest Descent (SD)

경도의 반대 방향으로 고정 스텝만큼 이동한다.

$$\mathbf{x}_{i+1} = \mathbf{x}_i - \alpha\cdot\nabla r(\mathbf{x}_i)$$

- $\alpha$ : 고정 스텝 크기. 코드상 `lambda=1e-13` (매우 작음)
- 수렴은 보장되지만 매우 느림

코드: `GearSteepestDescent.m`

> [!tip] 두 방법 비교 (논문 결론)
> | | LM | SD |
> |---|---|---|
> | 속도 | 빠름 (20~130 iter) | 느림 (8000+ iter) |
> | 안정성 | 불리한 조건에서 멈출 수 있음 | 수렴 보장 |
> | 결과 | 국소 최소에 갇힐 수 있음 | 더 큰 체적 감소 달성 가능 |
>
> → 실행 결과에서 LM은 1% 감소, SD는 37% 감소를 달성한 것이 이 차이를 직접 보여준다.

### 2.3 야코비 행렬 (Jacobian Matrix)

두 알고리즘 모두 $J(\mathbf{x})$ (21×7 행렬)가 필요하다. 각 성분은 잔차의 각 설계변수에 대한 편미분이다.

$$J_{ij} = \frac{\partial r_i(\mathbf{x})}{\partial x_j}$$

- **체적 미분** (`GearVolumed*.m`): 수치미분 (전진차분, $\Delta=0.001$) — 해석식 불가능
- **응력 미분** (`GearStressEx/Ind*.m`): 해석식 — MATLAB 심볼릭으로 유도
- **페널티 미분** (`GearPenalty*Jac.m`): 체인룰로 위 미분들을 조합

코드: `GearJac.m`

---

## 3. 코드 파일 역할 정리

### 3.1 메인 실행

| 파일 | 역할 |
|---|---|
| **`GearRun.m`** | **메인 스크립트.** 시작점 → LM → SD → 형상 그리기 → 결과 비교 |

### 3.2 핵심 함수 (Core)

| 파일 | 역할 |
|---|---|
| `GearObjective.m` | 잔차 벡터 $r(\mathbf{x})$ (21×1) 반환 |
| `GearJac.m` | 야코비 행렬 $J(\mathbf{x})$ (21×7) 반환 |
| `GearVolume.m` | 체적 수치적분 ($V_{GEAR}$) |
| `GearRho.m` | 곡률반경 $\rho(\alpha)$ 해석식 |
| `GearAlpha.m` | $\rho$ 극값 각도 뉴턴법 탐색 |
| `GearGetConstant.m` | 모든 한계 상수 중앙 관리 |

### 3.3 최적화 엔진

| 파일 | 역할 |
|---|---|
| `GearSimpleLevMar.m` | Levenberg-Marquardt 반복 |
| `GearSteepestDescent.m` | Steepest Descent 반복 |

### 3.4 응력 계산 + 해석 미분

| 파일 | 역할 |
|---|---|
| `GearStressExPos/Neg.m` | 외부 슬리브 접촉응력 (lobe/pit) |
| `GearStressIn.m` | 내부 슬리브 접촉응력 |
| `GearStressExd{e,q,m,h}Pos/Neg.m` | 외부 응력의 각 변수 편미분 (해석식) |
| `GearStressInd{h,Rs,Rh,Rw}.m` | 내부 응력의 각 변수 편미분 (해석식) |

### 3.5 체적 미분 (수치미분)

| 파일 | 역할 |
|---|---|
| `GearVolumed{e,q,m,h,Rs,Rh}.m` | 체적의 각 변수 편미분 (전진차분 $\Delta=0.001$) |

### 3.6 페널티 함수 (제약조건)

| 파일 | 역할 |
|---|---|
| `GearPenalty{Vmin,Emin,Mmin,Mmax,Hmin,Rsmin,Rhmin,qmin}.m` | 박스 제약 페널티 |
| `GearPenalty{Rwmin,Rwmax}.m` | $R_W$ 범위 페널티 |
| `GearPenaltyRho{p,n}{Min,Max}.m` | 곡률반경 범위 페널티 |
| `GearPenaltyStress{ExPos,ExNeg,In}.m` | 접촉응력 한계 페널티 |
| `...Jac.m` | 위 각 페널티의 야코비 성분 |

### 3.7 보조/검증용

| 파일 | 역할 |
|---|---|
| `GearDraw.m` | 최적화 전(빨강)/후(파랑) 기어 형상 plot |
| `GearRosenbrock.m`, `GearRosenbrockJac.m`, `GearTestRosenbrock.m` | LM 알고리즘 검증용 Rosenbrock 테스트 |
| `GearRhoPrim.m`, `GearRhoBis.m` | 곡률반경 1·2차 수치미분 |

> [!warning] 주의 파일
> `GearNewton.m`은 `GearSymbolicsDiff`라는 존재하지 않는 함수를 호출한다. 메인 흐름에서는 `GearAlpha`(자체 뉴턴법)를 쓰므로 `GearNewton`은 호출하지 말 것.

---

## 4. MATLAB에서 따라하기 (Step-by-Step)

### Step 1 — 파일 준비

> [!todo] 폴더 세팅
> 1. Zenodo(`DOI:10.5281/zenodo.2166718`)에서 zip을 받아 76개 `.m` 파일을 한 폴더에 둔다.
> 2. MATLAB에서 그 폴더를 **Current Folder**로 설정한다. 이미 해당 폴더 안에 있으면 별도 `addpath` 불필요.
> 3. `which GearRun`으로 인식 여부 확인.

### Step 2 — LM 알고리즘 검증 (Rosenbrock)

```matlab
GearTestRosenbrock
% xopt ≈ [1; 1], GearValueNext ≈ 0 이면 LM 정상
```

### Step 3 — 메인 실행

```matlab
GearRun
% Figure 1: LM 결과 (빨강=전, 파랑=후)
% Figure 2: SD 결과
```

### Step 4 — 결과 확인

```matlab
VolumePrev      % 시작점 체적
VolumeNext      % LM 후 체적
VolumeFin2      % SD 후 체적
StressExPosNext % LM 후 lobe 응력 (≤400 확인)
StressExNegNext % LM 후 pit 응력
StressInNext    % LM 후 내부 슬리브 응력
xfin            % LM 최적해 [e,q,m,h,Rs,Rh,Rw]
xfin2           % SD 최적해
```

> [!check] 확인 포인트
> - 체적: `VolumeNext < VolumePrev`, `VolumeFin2 < VolumePrev`
> - 응력: 세 값 모두 ≤ 400 MPa (약간 초과는 페널티 방식의 한계, 정상)
> - 설계변수: `xfin(3)` ($m$)이 0.5~0.85, `xfin(6)` ($R_h$) > `xfin(5)` ($R_S$)

### Step 5 — 여러 시작점으로 논문 Table 2 재현

```matlab
starts = [2.8, 3.5, 0.7, 10, 5, 8, 30;
          1.5, 2.0, 0.65, 8, 3, 5, 20;
          2.9, 3.4, 0.7, 10, 5, 8, 30;
          3.0, 3.3, 0.7, 10, 5, 8, 30];

for i = 1:size(starts,1)
    xp = starts(i,:)';
    xfin = GearSimpleLevMar(3000, @GearObjective, @GearJac, xp);
    V0 = GearVolume(xp(1)*16/xp(3), xp(1),16,15, xp(2),xp(3),1000, xp(4),xp(5),xp(6),8);
    V1 = GearVolume(xfin(1)*16/xfin(3), xfin(1),16,15, xfin(2),xfin(3),1000, xfin(4),xfin(5),xfin(6),8);
    fprintf('Variant %d:  V0=%.0f -> V1=%.0f  (%.1f%%)\n', i, V0, V1, 100*(V0-V1)/V0);
end
```

> [!warning] 시작점이 가용 영역 밖이면
> variant 2처럼 시작점에서 이미 응력이 400 MPa를 크게 초과하면, LM이 페널티 영역 안에서 멈춰 **체적이 오히려 늘어날 수 있다.** 경도 기반 방법의 국소 최소점 문제이며, 논문에서도 그렇게 설명한다.

---

## 5. 검증 체크리스트

> [!check] 구현이 제대로 됐는지 확인할 것들
> - [ ] `GearTestRosenbrock` → `[1,1]` 수렴
> - [ ] `GearObjective(xp)` → 21×1 벡터 (`size` 확인)
> - [ ] `GearJac(xp)` → 21×7 행렬 (`size` 확인)
> - [ ] `GearVolume` 반환값이 양수이고 합리적인 크기
> - [ ] LM 결과 `xfin`: $0.5 \le m \le 0.85$, $R_h > R_S$
> - [ ] 최적화 후 세 접촉응력 ≤ 400 MPa
> - [ ] `GearDraw` 그림에서 파란 곡선이 빨간 곡선보다 작은지

---

## 6. Troubleshooting

> [!bug] 흔한 오류
> - **`GearSymbolicsDiff` 에러**: `GearNewton.m` 호출 금지. 메인 흐름에서 안 쓰임.
> - **`Undefined function`**: 76개 파일이 모두 현재 폴더에 있는지 확인.
> - **수렴 안 함 / NaN**: $m$이 0.85 초과하거나 0.5 미만. 시작점 $m$ 조정.
> - **응력 발산**: 곡률반경 $\rho \approx 0$이면 응력 폭발. 시작점이 나쁜 경우 초기 발산 가능.
> - **역행렬 불안정**: $J^TJ$가 특이(singular)이면 LM 업데이트 불안정. `lambda` 초기값을 키울 것.

---

## 7. 다음 단계 (추후 진행)

> [!note] 확장 방향
> 1. **`fmincon` 버전 (SQP)**: 체적만 목적함수로, 페널티 항들을 명시적 부등호 제약 `g(x)≤0`으로 분리 → KKT 조건·라그랑지 승수 해석 (교재 직접 연결)
> 2. **다목적 최적화 (체적+효율)**: Wang et al. 효율 수식 추가 → 가중합법으로 `fmincon` 반복 → Pareto front (교재 18.4절)
> 3. **후최적성 해석**: 활성 제약의 라그랑지 승수로 "어떤 제약이 체적 감소를 가장 제한하는가" 분석

---

## 부록 A — 기호 빠른 참조

| 기호 | 의미 | 기호 | 의미 |
|---|---|---|---|
| $e$ | 편심률 (입력축-디스크 중심 거리) | $\rho$ | 곡률반경 |
| $q$ | 외부 슬리브 반경 | $\eta$ | 탄성계수 조합 |
| $m$ | 단폭계수 (0.5~0.85) | $F_{EX}$ | 외부 슬리브 힘 |
| $h$ | 기어 두께 | $F_{IN}$ | 내부 슬리브 힘 |
| $R_S$ | 내부 슬리브 반경 | $M_h$ | 입력 토크 |
| $R_h$ | 홀 반경 ($> R_S$) | $M_c$ | 출력 토크 ($=M_h z_s$) |
| $R_W$ | 홀 위치 반경 | $r$ | 피치원 반경 = 외부 슬리브 위치 반경 |
| $z_s$ | 치형 수 (15) | $z_k$ | 외부 슬리브 수 (16 = $z_s+1$) |
| $z_i$ | 내부 슬리브 수 (8) | $\alpha$ | 치형 곡선 매개변수 (0~2π) |

## 부록 B — 참고 자료

- **주 논문**: Król et al. (2019). *Optimization of the one stage cycloidal gearbox as a non-linear least squares problem*. MMS 73, pp. 1039–1048.
- **코드**: Król. *MATLAB scripts for cycloidal gearbox optimization*. `DOI:10.5281/zenodo.2166718`
- **교재**: Arora. *Introduction to Optimum Design*, 4th Ed. — Ch. 11(경도법), 페널티함수법, Ch. 18(다목적)
- **비교 논문 1**: Wang et al. (2016). *Multi-objective optimal design of cycloid speed reducer based on GA*. MMT 102, 135–148.
- **비교 논문 2**: Qi et al. (2024). *Design principle and numerical analysis for cycloidal drive*. Alexandria Eng. J. 91, 403–418.
