---
title: Wang et al. (2016) 논문 정리 + fmincon(SQP) 구현 계획
tags:
  - 최적설계
  - 사이클로이드감속기
  - Wang2016
  - fmincon
  - SQP
  - KKT
created: 2026-06-06
참고논문: "Wang, Luo, Su (2016), Multi-objective optimal design of cycloid speed reducer based on genetic algorithm, MMT 102, 135–148"
---

# Wang et al. (2016) 논문 정리 + fmincon 구현 계획

> [!info] 이 문서의 목적
> Wang et al.(2016) 논문의 수식 정식화를 정리하고, 이를 바탕으로 **`fmincon`(SQP) + 단목적(체적 최소화)** 버전을 구현하는 계획을 담은 문서.
> Król 논문(LM/SD + 페널티 함수법)과 함께 수업 내용(경도법, SQP, KKT 조건, 후최적성 해석)을 모두 커버하는 프로젝트 구성을 목표로 함.

---

## 0. 지금까지의 맥락 — 왜 Wang 논문인가

### 프로젝트 진행 경과

1. **처음 계획**: Wang et al.(2016) GA 기반 다목적 최적설계 → 수업에서 배우지 않은 GA라서 보류
2. **Król et al.(2019) 채택**: LM/SD + 페널티 함수법 → 수업 내용(경도법, 페널티 함수)과 잘 맞음
3. **Król LM/SD 재현 성공**: 논문 Table 1, 2 결과 재현 완료
4. **fmincon(SQP) 추가 시도**: Král 수식을 `fmincon`으로 변환하려 했으나 계속 실패
   - 원인: Král의 체적 함수가 수치적분 기반이라 `ub=inf`에서 발산
   - 시작점이 가용 영역 밖 (곡률반경 제약 초기 위반)
   - 페널티 방식과 명시적 제약 방식의 구조적 불일치

5. **결론**: Král 수식으로 `fmincon`을 안정적으로 돌리기 어렵고, Wang의 **폐형식(closed-form) 수식**이 `fmincon`에 적합함

### 왜 수업 내용과 연결이 필요한가

수업(Arora 교재)의 핵심 주제:
- 경도 기반 최적화 (Steepest Descent, Conjugate Gradient) → **Król LM/SD로 커버**
- 페널티 함수법 → **Král 코드로 커버**
- SQP (`fmincon`) → **아직 미구현**
- KKT 조건, 라그랑지 승수, 후최적성 해석 → **fmincon 출력으로 분석 가능**

Wang 수식은 폐형식이라 `fmincon`과 잘 맞고, XW3형 감속기 기준 초기값/범위도 Table 4에 명시되어 있어 안정적 수렴이 기대됨.

---

## 1. Wang et al. 논문 개요

**연구 대상**: K-H-V형 사이클로이드 감속기 (XW3형, 전달비 $z_c = 43$)

**연구 목표**: 체적 최소화 + 효율 최대화를 동시에 달성하는 다목적 최적설계

**방법**: 유전 알고리즘 (GA) — 단목적 최적화와 비교 분석도 수행

**주요 결과** (다목적 GA 기준):
- 효율: 0.862 → 0.8947 (**+3.79%**)
- 체적: 101,204 → 68,301 mm³ (**-32.51%**)

---

## 2. 감속기 구조 (K-H-V형)

Wang 논문의 감속기는 Král과 다른 구조예요.

| 부품 | 수량 | 역할 |
|---|---|---|
| 사이클로이드 기어 (g) | 1개 | 핀 기어와 맞물리는 주요 기어 |
| 핀 기어 (pin wheel, b) | 44개 | 내접 스퍼기어 역할 |
| 핀 슬리브 | 8개 | 핀과 핀 홀 사이 슬라이딩 베어링 |
| 핀 (pin) | 8개 | 출력 전달 |
| 이중 편심 슬리브 | 1개 | 사이클로이드 기어의 편심 운동 생성 |

> [!note] Král과의 구조 차이
> - Král: 슬리브 수가 설계변수 (내부 8개, 외부 16개)
> - Wang: 핀 44개, 핀 슬리브 8개 고정
> - Král: 싱글 디스크 / Wang: 전달비 43 기준 XW3형 실제 감속기

---

## 3. 정식화 (Problem Formulation)

### 3.1 설계변수

$$\mathbf{X} = [D_p,\ d_{rp},\ B,\ D,\ K_1,\ D_w,\ d_{sw}]^T$$

| 기호 | 코드 | 의미 | 단위 | 초기값 | 범위 |
|---|---|---|---|---|---|
| $D_p$ | `x(1)` | 핀 중심원 직경 | mm | 144 | 140~155 |
| $d_{rp}$ | `x(2)` | 핀 직경 | mm | 10 | 7~10.4 |
| $B$ | `x(3)` | 사이클로이드 기어 폭 | mm | 11 | 7~12 |
| $D$ | `x(4)` | 사이클로이드 기어 중심홀 직경 | mm | 53.5 | 50~55 |
| $K_1$ | `x(5)` | 단폭계수 | – | 0.6069 | 0.65~0.9 |
| $D_w$ | `x(6)` | 핀 중심원 직경 (출력부) | mm | 90 | 88~104 |
| $d_{sw}$ | `x(7)` | 핀 직경 (출력부) | mm | 12 | 11~14 |

**고정 상수:**

| 기호 | 값 | 의미 |
|---|---|---|
| $z_c$ | 43 | 사이클로이드 기어 치형 수 (전달비) |
| $z_p$ | 44 | 핀 기어 수 ($= z_c + 1$) |
| $z_w$ | 8 | 출력 핀 수 |
| $\mu$ | 0.05~0.1 | 핀-기어 마찰계수 |
| $f_w$ | 0.008~0.08 | 출력부 마찰계수 |
| $\eta_{gx}$ | 0.995 | 롤링 베어링 효율 |
| $\eta_{zx}$ | 0.99 | 피벗 베어링 효율 |
| $\Delta_2$ | – | 핀 슬리브 벽 두께 |
| $n$ | 1440 | 입력 회전수 (rpm) |
| $P$ | 0.75 | 입력 동력 (kW) |
| $M$ | – | 출력 토크 (N·mm) |

---

### 3.2 목적함수

**목적함수 1 — 효율 최대화** (= $1-\eta$ 최소화, 식 9):

$$\min f_1 = 1 - \frac{1-\dfrac{(D_p-d_{rp})\cdot 4\mu}{K_1 z_c D_p \pi}}{1+\dfrac{(D_p-d_{rp})\cdot 4\mu}{K_1 D_p \pi}} \cdot \left(1 - \frac{4f_w K_1 d_{sw} D_p}{\pi D_w (d_{sw}+2\Delta_2)}\right) \cdot \eta_{zx} \eta_{gx}^2$$

**목적함수 2 — 체적 최소화** (식 10):

$$\min f_2 = \frac{1}{4}\pi B\left[\left(D_p - K_1\frac{D_p}{z_p} - d_{rp}\right)^2 - \left(d_{sw} + 2\Delta_2 + K_1\frac{D_p}{z_p}\right)^2 z_w - D^2\right] + K_1\frac{D_p}{z_p}z_c B$$

> [!note] 두 목적함수 모두 폐형식
> 수치적분 없이 설계변수의 사칙연산으로 바로 계산됩니다. 이것이 Král 수식과의 가장 큰 차이이고, `fmincon`이 안정적으로 동작할 수 있는 이유입니다.

> [!note] 효율식의 구성
> $$\eta = \eta_x \cdot \eta_{zx} \cdot \eta_{gx}^2 \cdot \eta_{sx}$$
> - $\eta_x$: 핀-사이클로이드 맞물림 효율 (식 4)
> - $\eta_{sx}$: 출력부 효율 (식 5) — 핀-핀홀 마찰
> - $\eta_{zx} = 0.99$, $\eta_{gx} = 0.995$: 베어링 효율 (상수)

---

### 3.3 제약조건 (16개)

#### 기하학적 제약

**y1 — 언더컷/첨예 방지** (식 12): 두 경우로 나뉨

$$\frac{d_{rp}}{2} - \frac{D_p}{2}\sqrt{\frac{27(1-K_1^2)(z_p-1)}{(z_p+1)^3}} < 0 \qquad \text{if } \frac{z_c-1}{2z_c+1} < K_1 < 1$$

$$\frac{d_{rp}}{2} - \frac{D_p(1-K_1)^2}{2(z_pK_1+1)} < 0 \qquad \text{if } K_1 \le \frac{z_c-1}{2z_c+1}$$

**y2, y3 — 단폭계수 범위** (식 13, 14):

$$c_1 - K_1 \le 0, \quad K_1 - c_2 \le 0 \qquad (c_1=0.65,\ c_2=0.9 \text{ for } z_c=25\text{~}59)$$

**y4, y5 — 핀직경 계수 $K_2$ 범위** (식 16, 17):

$$d_1 - \frac{D_p}{d_{rp}}\sin\frac{\pi}{z_p} \le 0, \quad \frac{D_p}{d_{rp}}\sin\frac{\pi}{z_p} - d_2 \le 0 \qquad (d_1=1.6,\ d_2=1.0 \text{ for } z_p=36\text{~}60)$$

#### 강도 제약

**y6 — 사이클로이드 기어 접촉강도** (식 22):

$$0.418\sqrt{\frac{E_e \cdot 4.4M}{B K_1 z_c D_p \rho_{emin}}} - \sigma_{HP} \le 0 \qquad (\sigma_{HP} = 1000 \text{ MPa})$$

**y7 — 핀 기어 굽힘강도** (식 25):

$$\frac{1.41 \times 4.4 \times 9550PL}{K_1 D_p n d_{rp}^2} - \sigma_{FP} \le 0 \qquad (D_p < 390, \quad \sigma_{FP} = 150 \text{ MPa})$$

**y8 — 핀-핀홀 접촉강도** (식 27):

$$300\sqrt{\frac{K_1 M D_p}{z_w D_w B \left(\frac{d_{sw}}{2}+\Delta_2\right)^2 z_p + \frac{1}{2}K_1 D_p\left(\frac{d_{sw}}{2}+\Delta_2\right)}} - \sigma_{ZHP} \le 0$$

**y9 — 핀 굽힘강도** (식 29):

$$\frac{4.4 K_w M (1.5B + \Delta_c)}{0.1 z_w D_w d_{sw}^3} - \sigma_{BP} \le 0 \qquad (\sigma_{BP} = 200 \text{ MPa},\ K_w = 1.4)$$

#### 치수 제약

**y10, y11 — 핀 중심원 직경 범위** (식 30, 31):

$$e_1 - D_p \le 0, \quad D_p - e_2 \le 0 \qquad (e_1=140,\ e_2=155)$$

**y12, y13 — 핀홀 최대 직경** (식 32, 33):

$$0.06D_p - D_w + d_w + D \le 0$$
$$0.03D_p - D_w\sin\frac{\pi}{z_w} + d_w \le 0$$

여기서 $d_w = d_{sw} + 2\Delta_2 + 2a + \Delta$ ($a$: 중심 거리, $\Delta=0.15$ mm: 핀홀-슬리브 간격)

**y14, y15 — 사이클로이드 기어 폭 범위** (식 35, 36):

$$0.05D_p - B \le 0, \quad B - 0.1D_p \le 0$$

**y16 — 피벗 베어링 수명** (식 38):

$$5000 - \frac{10^6}{60n_1}\left(\frac{C}{p}\right)^{10/3} \le 0$$

여기서 $C$, $p$는 베어링 정격하중과 동적하중으로, $M$, $D_p$, $K_1$, $z_c$, $z_p$, $n$의 함수.

---

## 4. 논문 결과 요약

Wang 논문은 세 가지 최적화를 비교합니다.

| 방법 | 체적 (mm³) | 효율 | 체적 감소 | 효율 증가 |
|---|---|---|---|---|
| 초기 설계 | 101,204 | 0.862 | – | – |
| 체적 최소화 (단목적) | 64,693 | 0.875 | **-36.1%** | +1.5% |
| 효율 최대화 (단목적) | 93,229 | 0.896 | -7.9% | **+3.9%** |
| **다목적 (체적+효율)** | **68,301** | **0.895** | **-32.5%** | **+3.8%** |

**결론**: 단목적 최적화는 하나의 목적만 최적화하면서 다른 성능을 희생하거나 무시하는 경향이 있음. 다목적 최적화가 두 목표를 균형 있게 달성하여 더 robust한 해를 제공함.

---

## 5. fmincon(SQP) 구현 계획

### 5.1 왜 단목적(체적 최소화)인가

수업 교재(Arora)의 `fmincon` 예제들은 모두 단일 목적함수 + 명시적 제약 구조입니다. 다목적은 교재 18장에서 **가중합법**으로 다루는데, 이는 단목적 `fmincon`을 가중치를 바꿔가며 반복하는 방식이에요. 따라서:

1. **1단계**: 체적만 최소화하는 단목적 `fmincon` 구현 → SQP, KKT 해석
2. **2단계 (선택)**: 가중합법으로 다목적 확장 → Pareto front

### 5.2 구현 구조

**`WangObj_fmincon.m`** — 목적함수 (체적):

```matlab
function V = WangObj_fmincon(x)
    % x = [Dp, drp, B, D, K1, Dw, dsw]
    zp = 44; zc = 43; zw = 8;
    Delta2 = 2;  % 핀 슬리브 벽 두께 (mm, 추정값)
    
    Dp=x(1); drp=x(2); B=x(3); D=x(4); K1=x(5); Dw=x(6); dsw=x(7);
    
    V = (1/4)*pi*B * ( (Dp - K1*Dp/zp - drp)^2 ...
                     - (dsw + 2*Delta2 + K1*Dp/zp)^2 * zw ...
                     - D^2 ) ...
        + K1*(Dp/zp)*zc*B;
end
```

**`WangCon_fmincon.m`** — 제약함수 (16개):

```matlab
function [c, ceq] = WangCon_fmincon(x)
    % 16개 부등호 제약 c(x) <= 0
    % 등호 제약 없음
    ...
    ceq = [];
end
```

**`WangRun_fmincon.m`** — 메인 실행:

```matlab
x0 = [144, 10, 11, 53.5, 0.6069, 90, 12]';  % 논문 Table 4 초기값
lb = [140, 7,   7,  50, 0.65, 88, 11]';       % 논문 Table 4 하한
ub = [155, 10.4, 12, 55, 0.9,  104, 14]';      % 논문 Table 4 상한

options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'OptimalityTolerance', 1e-6, 'ConstraintTolerance', 1e-6);

[xopt, Vopt, exitflag, ~, lambda] = fmincon(@WangObj_fmincon, x0, ...
    [], [], [], [], lb, ub, @WangCon_fmincon, options);
```

### 5.3 수업 내용 연결 포인트

**SQP 알고리즘**: `fmincon`이 내부적으로 매 반복마다 QP(이차계획법) 문제를 풀며 이동 방향을 결정. Arora 교재 12장.

**KKT 조건**: 최적점 $\mathbf{x}^*$에서 `lambda.ineqnonlin`이 라그랑지 승수 $\lambda_i$를 제공.
$$\nabla f(\mathbf{x}^*) + \sum_i \lambda_i \nabla g_i(\mathbf{x}^*) = 0, \quad \lambda_i \ge 0, \quad \lambda_i g_i(\mathbf{x}^*) = 0$$

**후최적성 해석 (Post-optimality analysis)**: 활성 제약($\lambda_i > 0$)의 $\lambda_i$ 크기로 "어떤 제약이 체적 감소를 가장 제한하는가" 분석. Arora 교재 7장.

**Wang GA 결과와 비교**: `fmincon` 결과(체적)를 Wang Table 6(단목적 체적 최소화)과 비교하면 SQP vs GA 방법론 비교가 됨.

### 5.4 Wang 수식이 fmincon에 적합한 이유

Král 수식과의 결정적 차이:

| | Král 체적 | Wang 체적 |
|---|---|---|
| 계산 방식 | 수치적분 ($\int$) | 폐형식 (사칙연산) |
| 도함수 | 수치미분 필요 | 해석적 도함수 가능 |
| 발산 위험 | $R_h$ 증가 시 음수 가능 | 물리적 한계 내에서 안정 |
| `fmincon` 적합성 | 낮음 | 높음 |

---

## 6. 전체 프로젝트 구성 (최종)

```
사이클로이드 감속기 최적설계 프로젝트
│
├── [1] Król et al. (2019) — 페널티 함수법 + 경도법
│     ├── LM (Levenberg-Marquardt) 구현 및 재현
│     ├── SD (Steepest Descent) 구현 및 재현
│     ├── 여러 시작점 실험 → 국소 최솟값 문제 분석
│     └── 수업 연결: 경도법, 페널티 함수법
│
└── [2] Wang et al. (2016) 수식 — fmincon(SQP)
      ├── 단목적(체적 최소화) fmincon 구현
      ├── KKT 조건 검증 및 라그랑지 승수 해석
      ├── Wang GA 결과(Table 6)와 비교
      └── 수업 연결: SQP, KKT 조건, 후최적성 해석
```

**비교 가능한 포인트:**

| 비교 항목 | Král + LM/SD | Wang + fmincon |
|---|---|---|
| 최적화 알고리즘 | 경도 기반 (LM, SD) | SQP |
| 제약 처리 | 페널티 함수법 | 명시적 부등호 제약 |
| 설계변수 | 치형 기하학 (7개) | 감속기 치수 (7개) |
| 수업 연결 | 경도법, 페널티 함수 | SQP, KKT, 후최적성 |
| 발산/수렴 특성 | 시작점 의존성 높음 | 안정적 수렴 기대 |

---

## 부록 — 기호 빠른 참조

| 기호 | 의미 |
|---|---|
| $D_p$ | 핀 중심원 직경 (pin wheel central circle) |
| $d_{rp}$ | 핀 직경 (pin wheel diameter) |
| $B$ | 사이클로이드 기어 폭 |
| $D$ | 사이클로이드 기어 중심홀 직경 |
| $K_1$ | 단폭계수 (short width coefficient) |
| $D_w$ | 핀 중심원 직경 (출력부, pin central circle) |
| $d_{sw}$ | 핀 직경 (출력부) |
| $z_c$ | 사이클로이드 기어 치형 수 (= 전달비 = 43) |
| $z_p$ | 핀 기어 수 ($= z_c + 1 = 44$) |
| $z_w$ | 출력 핀 수 (= 8) |
| $K_2$ | 핀직경 계수 |
| $\rho_{emin}$ | 최소 등가 곡률반경 |
| $\mu$ | 핀-기어 마찰계수 |
| $f_w$ | 출력부 마찰계수 |
| $\Delta_2$ | 핀 슬리브 벽 두께 |
| $\Delta_c$ | 스페이서 링 두께 |
| $K_w$ | 하중 분배 계수 (= 1.4) |

---

## 참고 자료

- **Wang et al. (2016)**: *Multi-objective optimal design of cycloid speed reducer based on genetic algorithm*, MMT 102, 135–148.
- **Král et al. (2019)**: *Optimization of the one stage cycloidal gearbox as a non-linear least squares problem*, MMS 73, pp. 1039–1048.
- **Arora (2016)**: *Introduction to Optimum Design*, 4th Ed. — Ch. 7(fmincon), Ch. 12(SQP), Ch. 3/4(KKT), Ch. 18(다목적)
