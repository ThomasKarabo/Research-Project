# CV-MTSE and BL-MF Portfolio Optimisation
## Exhaustive methodological and implementation guide for the UNISA South African equity research project

> **Purpose.** This document explains, step by step, how the Cross-Validation Multi-Target Shrinkage Estimator (CV-MTSE) and the Black–Litterman Multi-Factor (BL-MF) fusion model work, using the two uploaded papers as the primary methodological references and adapting them to the current South African research design.
>
> **Important distinction:** the papers provide the methodological foundation. The South African data, 2010–2019 training period, 2020–2025 out-of-sample period, Satrix 40 market benchmark, SARB 91-day Treasury Bill risk-free proxy, and two-factor MKT-RF/SMB specification are project-specific adaptations.

---

# 1. Research architecture

The current research pipeline is:

\[
\boxed{\text{EW} \rightarrow \text{MVO} \rightarrow \text{CV-MTSE} \rightarrow \text{BL-MF} \rightarrow \text{Integrated Model}}
\]

The roles are:

| Model | Main purpose |
|---|---|
| Equal Weight (EW) | Simple benchmark |
| Mean-Variance Optimisation (MVO) | Traditional optimisation benchmark |
| CV-MTSE | Improve covariance/risk estimation |
| BL-MF | Improve expected-return estimation by fusing Black–Litterman information with factor information |
| Integrated Model | Use BL-MF expected returns together with CV-MTSE covariance estimation |

EW and MVO are **benchmarks**, not components of the integrated model.

The research must preserve a strict chronological separation:

- **Training:** 2010–2019
- **Out-of-sample testing:** 2020–2025

No information from 2020–2025 may be used when estimating parameters for the beginning of the out-of-sample period.

---

# 2. Why these two methods solve different problems

Portfolio optimisation requires at least two major inputs:

\[
\mu = E(r)
\]

and

\[
\Sigma = \operatorname{Cov}(r).
\]

Here:

- \(\mu\) is the vector of expected excess returns.
- \(\Sigma\) is the covariance matrix of asset excess returns.

These inputs are difficult to estimate accurately.

## 2.1 CV-MTSE

CV-MTSE primarily addresses the **covariance estimation problem**.

Instead of trusting the noisy sample covariance matrix completely, it combines it with structurally simpler covariance targets.

Conceptually:

\[
\text{CV-MTSE}
\Rightarrow
\text{better estimate of }\Sigma.
\]

## 2.2 BL-MF

BL-MF primarily addresses the **expected-return/factor-estimation problem**.

It starts with Black–Litterman information and then treats the resulting BL expected-return estimate as an additional noisy observation that can be fused with a multi-factor model.

Conceptually:

\[
\text{BL-MF}
\Rightarrow
\text{better estimate of factor returns}
\Rightarrow
\text{better estimate of }\mu\text{ and }\Sigma.
\]

For the final integrated model:

\[
\boxed{
\text{BL-MF expected returns}
+
\text{CV-MTSE covariance}
\Rightarrow
\text{portfolio optimisation}
}
\]

---

# 3. CV-MTSE: mathematical foundation

## 3.1 Start with asset returns

Let there be \(N\) assets and \(T\) observations.

Define:

\[
r_t =
\begin{bmatrix}
r_{1,t}\\
r_{2,t}\\
\vdots\\
r_{N,t}
\end{bmatrix}.
\]

The sample mean vector is:

\[
\bar r =
\frac{1}{T}
\sum_{t=1}^{T}r_t.
\]

The sample covariance matrix is:

\[
\widehat{\Sigma}_{SCM}
=
\frac{1}{T-1}
\sum_{t=1}^{T}
(r_t-\bar r)(r_t-\bar r)'.
\]

The sample covariance matrix is attractive because it uses the data directly, but it can be unstable when the number of assets is large relative to the number of observations.

That instability is particularly important in portfolio optimisation because small errors in \(\widehat{\Sigma}\) can produce large changes in portfolio weights.

---

# 4. The idea of shrinkage

A shrinkage estimator deliberately moves the noisy sample covariance matrix toward a more stable target.

For a single target \(F\):

\[
\widehat{\Sigma}_{shrink}
=
(1-\delta)\widehat{\Sigma}_{SCM}
+
\delta F,
\]

where:

\[
0\leq\delta\leq1.
\]

Interpretation:

- \(\delta=0\): use only the sample covariance matrix.
- \(\delta=1\): use only the target.
- \(0<\delta<1\): combine the sample covariance with the target.

The key problem is deciding **how much shrinkage should be applied**.

---

# 5. Multi-target shrinkage

The CV-MTSE paper extends this idea to several targets.

With \(M\) target matrices:

\[
\boxed{
\widehat{\Sigma}(\delta)
=
\left(1-\sum_{i=1}^{M}\delta_i\right)
\widehat{\Sigma}_{SCM}
+
\sum_{i=1}^{M}\delta_iF_i
}
\]

where:

\[
0\leq\delta_i\leq1
\]

and

\[
\sum_{i=1}^{M}\delta_i\leq1.
\]

The coefficient on the sample covariance matrix is therefore:

\[
\delta_{SCM}
=
1-\sum_{i=1}^{M}\delta_i.
\]

The three weights in the project-specific two-target case are therefore:

\[
w_{SCM}=1-\delta_{SIM}-\delta_{IM},
\]

\[
w_{SIM}=\delta_{SIM},
\]

\[
w_{IM}=\delta_{IM}.
\]

They sum to one.

This convex-combination structure helps preserve symmetry and positive definiteness when the component matrices satisfy the required conditions.

---

# 6. The two targets used by CV-MTSE

The reference paper uses two targets:

1. Single Index Model (SIM)
2. Identity Matrix (IM)

Thus:

\[
M=2.
\]

The estimator becomes:

\[
\boxed{
\widehat{\Sigma}_{CV-MTSE}
=
(1-\delta_{SIM}-\delta_{IM})
\widehat{\Sigma}_{SCM}
+
\delta_{SIM}\widehat{\Sigma}_{SIM}
+
\delta_{IM}I
}
\]

where \(I\) is the identity matrix, subject to the relevant scaling/implementation convention discussed below.

---

# 7. Target 1: Single Index Model covariance matrix

The Single Index Model represents each asset's return as depending on a common market component plus an idiosyncratic component.

Conceptually:

\[
r_i
=
\alpha_i+\beta_i r_m+\epsilon_i.
\]

The covariance implied by the model is:

For \(i\neq j\):

\[
[\widehat{\Sigma}_{SIM}]_{ij}
=
\beta_i\beta_j\sigma_m^2.
\]

For \(i=j\):

\[
[\widehat{\Sigma}_{SIM}]_{ii}
=
\beta_i^2\sigma_m^2+\sigma_{\epsilon_i}^2.
\]

Here:

- \(\beta_i\) = asset \(i\)'s sensitivity to the market.
- \(\sigma_m^2\) = market-return variance.
- \(\sigma_{\epsilon_i}^2\) = residual/idiosyncratic variance of asset \(i\).

The important structural feature is that most cross-asset covariance comes through the common market factor.

---

# 8. Target 2: Identity matrix

The identity matrix is:

\[
I=
\begin{bmatrix}
1&0&\cdots&0\\
0&1&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
0&0&\cdots&1
\end{bmatrix}.
\]

It contains no off-diagonal covariance.

Its purpose is not to describe the entire market realistically. Its purpose is to provide a highly stable structural target that can reduce noisy covariance estimates.

### Important implementation issue

The paper describes the identity matrix as a target. In an implementation, the scale of the identity target matters because a literal unit-diagonal matrix may be on a different variance scale from the sample covariance matrix.

Therefore the implementation should reproduce the paper's stated convention exactly where possible and document any necessary scaling. Do **not** silently insert a different normalisation merely because it produces convenient numbers.

---

# 9. Why multiple targets are useful

The two targets represent different assumptions.

### SCM

Uses the observed data structure directly.

Advantage:
- flexible
- responsive

Weakness:
- noisy
- unstable in high dimensions

### SIM

Assumes covariance is largely explained by common market exposure.

Advantage:
- economically structured
- reduces estimation noise

Weakness:
- may oversimplify cross-sectional dependence

### IM

Assumes little useful covariance structure.

Advantage:
- extremely stable
- resistant to noisy off-diagonal estimates

Weakness:
- ignores correlations between assets

The multi-target estimator allows the data-driven procedure to decide how strongly to rely on each structure.

---

# 10. The shrinkage optimisation problem

Ideally, the shrinkage parameters would minimise the mean squared error:

\[
L(\delta)
=
E\left[
\left\|
\widehat{\Sigma}(\delta)-\Sigma
\right\|_F^2
\right].
\]

Here:

- \(\Sigma\) is the unknown true covariance matrix.
- \(\|\cdot\|_F\) is the Frobenius norm.

The ideal parameters would therefore be:

\[
\delta^*
=
\arg\min_{\delta}
L(\delta).
\]

But there is a fundamental problem:

\[
\Sigma
\]

is unknown.

Therefore the true MSE cannot be calculated directly.

This is exactly where cross-validation enters.

---

# 11. Cross-validation: the central idea

Cross-validation treats the shrinkage intensities as **hyperparameters**.

Instead of asking:

> Which shrinkage weights minimise error against the unknown true covariance matrix?

we ask:

> Which shrinkage weights produce a covariance estimate that performs best when evaluated against held-out data?

The paper therefore treats covariance estimation similarly to a machine-learning model-selection problem.

The data are divided into:

\[
D=D_{train}\cup D_{val}.
\]

The estimator is fitted using \(D_{train}\), and its covariance estimate is evaluated against information from \(D_{val}\).

The paper assumes the validation-set sample covariance is an approximation to the unknown true covariance for the purpose of selecting the hyperparameters.

---

# 12. Grid search

For two shrinkage parameters:

\[
\delta_{SIM},\delta_{IM}\in[0,1].
\]

The paper discretises each parameter as:

\[
\{0,0.05,0.10,\ldots,1.00\}.
\]

Thus there are 21 candidate values per parameter.

A naive full grid therefore contains:

\[
21^2=441
\]

candidate combinations for two targets.

Each candidate must satisfy:

\[
\delta_{SIM}+\delta_{IM}\leq1.
\]

The search therefore evaluates admissible pairs such as:

\[
(0,0),
(0.05,0),
(0,0.05),
(0.05,0.05),
\ldots
\]

but excludes combinations such as:

\[
(0.80,0.50)
\]

because their sum exceeds one.

---

# 13. One CV-MTSE iteration

Suppose the candidate is:

\[
\delta_{SIM}=0.30,\qquad
\delta_{IM}=0.20.
\]

Then:

\[
\delta_{SCM}=1-0.30-0.20=0.50.
\]

The candidate estimator is:

\[
\widehat{\Sigma}(0.30,0.20)
=
0.50\widehat{\Sigma}_{SCM}
+
0.30\widehat{\Sigma}_{SIM}
+
0.20I.
\]

The candidate covariance matrix is then evaluated against the validation data according to the paper's cross-validation loss construction.

The process is repeated for every admissible grid point.

The candidate with the smallest validation loss is selected.

---

# 14. Complete CV-MTSE algorithm

At a particular estimation date \(t\):

### Step 1 — Select the information set

Use only observations available before \(t\).

### Step 2 — Split the historical estimation sample

Create:

\[
D_{train}
\]

and

\[
D_{val}.
\]

For a time-series application, the split should preserve chronological ordering. A future validation block should not be allowed to leak into the training sample.

### Step 3 — Estimate the training SCM

Calculate:

\[
\widehat{\Sigma}_{SCM}^{train}.
\]

### Step 4 — Estimate the SIM target

Estimate market betas and residual variances using the training information.

Construct:

\[
\widehat{\Sigma}_{SIM}^{train}.
\]

### Step 5 — Construct the identity target

Construct the project's specified identity target.

### Step 6 — Generate the grid

Generate all admissible:

\[
(\delta_{SIM},\delta_{IM})
\]

pairs.

### Step 7 — Construct each candidate

For each candidate:

\[
\widehat{\Sigma}(\delta)
=
(1-\delta_{SIM}-\delta_{IM})
\widehat{\Sigma}_{SCM}^{train}
+
\delta_{SIM}\widehat{\Sigma}_{SIM}^{train}
+
\delta_{IM}I.
\]

### Step 8 — Evaluate on validation data

Calculate the validation loss for the candidate.

### Step 9 — Select the minimum-loss candidate

\[
\delta^*
=
\arg\min_{\delta}L_{CV}(\delta).
\]

### Step 10 — Construct the final covariance estimate

Using the selected shrinkage coefficients, produce the covariance matrix used by the portfolio optimiser.

---

# 15. Why CV-MTSE can adapt over time

If the procedure is repeated at every rebalancing date, the optimal parameters can change.

For example:

| Regime | SCM weight | SIM weight | IM weight |
|---|---:|---:|---:|
| Stable market | potentially lower | moderate | potentially higher |
| High-noise market | potentially lower | potentially higher | potentially higher |
| Rapid structural change | potentially higher | data-dependent | data-dependent |

These are conceptual interpretations, not guaranteed outcomes.

The important point is:

\[
\delta_t^*
\]

can vary with time.

Thus CV-MTSE is **adaptive** rather than permanently assigning one shrinkage coefficient.

---

# 16. CV-MTSE versus EW-MTSE

The reference paper also uses an equal-weighted multi-target estimator:

\[
\widehat{\Sigma}_{EW-MTSE}
=
\frac13\widehat{\Sigma}_{SCM}
+
\frac13\widehat{\Sigma}_{SIM}
+
\frac13I.
\]

This is a useful benchmark because its weights are fixed.

CV-MTSE instead estimates the weights from data.

The comparison therefore asks whether data-driven selection adds value over a simple one-third/one-third/one-third combination.

---

# 17. What CV-MTSE does and does not estimate

CV-MTSE is fundamentally a **covariance estimator**.

It does not, by itself, provide a sophisticated expected-return model.

Therefore, in the research pipeline:

\[
\boxed{
CV-MTSE\rightarrow\widehat{\Sigma}
}
\]

The expected-return vector must be supplied separately.

For the integrated model, the natural complementary source is BL-MF:

\[
\boxed{
BL-MF\rightarrow\widehat{\mu}
}
\]

and:

\[
\boxed{
Integrated:
(\widehat{\mu}_{BL-MF},\widehat{\Sigma}_{CV-MTSE})
}
\]

---

# 18. BL-MF: overall idea

The BL-MF paper begins with a problem in standard Black–Litterman.

Black–Litterman incorporates:

1. equilibrium expected returns;
2. investor views.

A multi-factor model provides:

1. systematic factor exposures;
2. factor return information;
3. residual/idiosyncratic risk decomposition.

The BL-MF model treats the Black–Litterman solution as an additional noisy observation of expected returns and fuses it with the multi-factor model.

The central idea is:

\[
\boxed{
\text{BL information}
+
\text{factor-model information}
\rightarrow
\text{fused factor-return estimate}
}
\]

---

# 19. Step 1 of BL-MF: the Black–Litterman model

Let:

\[
r\sim N(\mu,\Sigma)
\]

where:

- \(r\) = asset excess-return vector;
- \(\mu\) = expected excess-return vector;
- \(\Sigma\) = covariance matrix.

The BL prior assumes:

\[
\mu\sim N(\pi,C)
\]

with:

\[
C=\tau\Sigma.
\]

Here:

- \(\pi\) = equilibrium expected excess-return vector;
- \(C\) = uncertainty around the equilibrium returns;
- \(\tau>0\) = uncertainty scale.

---

# 20. Investor views

Investor views are represented as:

\[
P\mu=q+\eta,
\]

where:

- \(P\) = \(k\times N\) view matrix;
- \(q\) = vector of subjective views;
- \(\eta\) = view error.

The paper assumes:

\[
\eta\sim N(0,\Omega).
\]

Typically:

\[
\Omega
=
\operatorname{diag}
(\omega_1^2,\ldots,\omega_k^2).
\]

Thus each view has an explicit uncertainty level.

---

# 21. The Black–Litterman posterior

Combining the prior and views gives:

\[
\boxed{
\widehat{\mu}_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}
(C^{-1}\pi+P'\Omega^{-1}q)
}
\]

This is the Black–Litterman solution used by the BL-MF paper.

Its associated MSE matrix is:

\[
\boxed{
S_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}
}
\]

This second equation is crucial for BL-MF.

Why?

Because the BL-MF model does not merely use \(\widehat{\mu}_{BL}\).

It also uses:

\[
S_{BL}
\]

as a measure of how uncertain the BL estimate is.

---

# 22. Interpretation of the BL solution

The BL estimate is a precision-weighted combination of:

- the equilibrium prior;
- investor views.

The terms involving:

\[
C^{-1}
\]

represent the precision of the prior.

The terms involving:

\[
\Omega^{-1}
\]

represent the precision of the investor views.

Therefore:

- smaller view uncertainty \(\Omega\) → stronger influence of views;
- larger view uncertainty \(\Omega\) → weaker influence of views;
- smaller prior uncertainty \(C\) → stronger influence of equilibrium;
- larger prior uncertainty \(C\) → weaker influence of equilibrium.

This is the Bayesian logic underlying the posterior.

---

# 23. Investor views in the reference paper

The paper constructs views using momentum information.

It ranks stocks using recent six-month performance and constructs momentum portfolios.

The exact empirical construction is described in the paper's investor-view appendix.

For the South African research project, this must be treated as an **adaptation decision**, not as an automatic assumption.

If the project uses momentum-generated views, the view construction must use only information available at the relevant portfolio formation date.

No future return can enter \(P\), \(q\), or \(\Omega\).

---

# 24. Step 2 of BL-MF: the multi-factor model

The paper's factor model is:

\[
\boxed{
r=Xf+\epsilon
}
\]

where:

- \(r\) = \(N\times1\) asset excess-return vector;
- \(X\) = \(N\times M\) factor-loading matrix;
- \(f\) = \(M\times1\) factor-return vector;
- \(\epsilon\) = residual vector.

The model assumes:

\[
E(\epsilon)=0
\]

and:

\[
\operatorname{Cov}(\epsilon)=D.
\]

The residual covariance matrix is diagonal:

\[
D=
\operatorname{diag}
(\sigma_{\epsilon,1}^2,\ldots,\sigma_{\epsilon,N}^2).
\]

---

# 25. Project-specific two-factor model

The uploaded BL-MF paper uses a much larger factor set in its empirical application.

That is **not** the specification for this project.

The current research design explicitly uses two factors:

1. Market excess return, \(MKT-RF\)
2. SMB

Thus:

\[
f_t=
\begin{bmatrix}
MKT_t-R_{f,t}\\
SMB_t
\end{bmatrix}.
\]

The market factor is:

\[
\boxed{
MKT_t-R_{f,t}
=
R_{M,t}-R_{f,t}
}
\]

where:

- \(R_M\) = Satrix 40 market return;
- \(R_f\) = South African 91-day Treasury Bill rate.

SMB must be constructed from the South African equity universe using a documented market-capitalisation-based small-minus-big portfolio methodology.

No HML, MOM, RMW, CMA, or Carhart four-factor model should be added unless the research design is explicitly changed.

---

# 26. Factor-loading matrix

For two factors:

\[
X=
\begin{bmatrix}
\beta_{1,MKT}&\beta_{1,SMB}\\
\beta_{2,MKT}&\beta_{2,SMB}\\
\vdots&\vdots\\
\beta_{N,MKT}&\beta_{N,SMB}
\end{bmatrix}.
\]

For asset \(i\):

\[
r_{i,t}
=
\beta_{i,MKT}MKT_t
+
\beta_{i,SMB}SMB_t
+
\epsilon_{i,t}.
\]

The betas/loadings must be estimated using only information available at the estimation date.

---

# 27. Step 3 of BL-MF: reinterpret the BL estimate

This is the central innovation.

The BL-MF paper treats:

\[
\widehat{\mu}_{BL}
\]

as a noisy observation of the true expected excess-return vector:

\[
\boxed{
\widehat{\mu}_{BL}
=
\mu+v
}
\]

where:

\[
E(v)=0
\]

and:

\[
\operatorname{Cov}(v)=S_{BL}.
\]

Thus the BL estimate is treated as a measurement with known uncertainty.

This is what allows BL and the factor model to be placed into a common statistical framework.

---

# 28. Rewriting the BL measurement

The factor model implies:

\[
\mu=Xf
\]

at the level of expected returns.

The BL observation can therefore be written as:

\[
\widehat{\mu}_{BL}
=
Xf+(\mu+v-Xf).
\]

The term:

\[
\mu+v-Xf
\]

becomes an additional error term.

This allows the BL estimate to be written in the same structural form as the factor model.

---

# 29. Step 4: create the fused observation system

The original factor observation is:

\[
r=Xf+\epsilon.
\]

The BL observation is:

\[
\widehat{\mu}_{BL}
=
Xf+(\mu+v-Xf).
\]

Stack them:

\[
\boxed{
r_F=X_Ff+\epsilon_F
}
\]

where:

\[
r_F=
\begin{bmatrix}
\widehat{\mu}_{BL}\\
r
\end{bmatrix}
\]

and:

\[
X_F=
\begin{bmatrix}
X\\
X
\end{bmatrix}.
\]

The augmented error is:

\[
\epsilon_F=
\begin{bmatrix}
\mu+v-Xf\\
\epsilon
\end{bmatrix}.
\]

Thus the model has two information sources:

1. the BL estimate;
2. observed asset returns.

---

# 30. Covariance of the fused errors

The paper derives:

\[
\boxed{
\operatorname{Cov}(\epsilon_F)
=
\begin{bmatrix}
\Sigma-D+S_{BL}&0\\
0&D
\end{bmatrix}
}
\]

This is one of the most important equations in the entire BL-MF model.

The upper block represents uncertainty in the BL measurement relative to the factor model.

The lower block represents the residual noise in observed asset returns.

---

# 31. Why the error covariance matters

The model is not simply averaging the BL and factor estimates.

Instead, it performs a **weighted statistical fusion**.

A source with greater uncertainty receives less statistical weight.

This is visible in the inverse covariance matrices:

\[
(\Sigma-D+S_{BL})^{-1}
\]

and:

\[
D^{-1}.
\]

A covariance matrix with high uncertainty has lower precision.

A covariance matrix with lower uncertainty has higher precision.

Thus the model automatically determines how strongly each source contributes.

---

# 32. Step 5: estimate factor returns

The optimal weighted least-squares factor-return estimator is:

\[
\boxed{
\widehat f_{BL-MF}
=
\left[
X'(\Sigma-D+S_{BL})^{-1}X
+
X'D^{-1}X
\right]^{-1}
\left[
X'(\Sigma-D+S_{BL})^{-1}\widehat{\mu}_{BL}
+
X'D^{-1}r
\right]
}
\]

This is the central BL-MF estimation equation.

It combines:

- information from the BL expected-return estimate;
- information from observed asset returns;
- factor loadings;
- uncertainty in each information source.

---

# 33. Interpretation of the factor estimator

The estimator can be understood as a precision-weighted combination.

The BL contribution is:

\[
X'(\Sigma-D+S_{BL})^{-1}\widehat{\mu}_{BL}.
\]

The observed-return contribution is:

\[
X'D^{-1}r.
\]

The final estimate is scaled by the combined information matrix:

\[
X'(\Sigma-D+S_{BL})^{-1}X
+
X'D^{-1}X.
\]

Thus the model is not saying:

> BL is always better than the factor model.

Nor is it saying:

> The factor model is always better than BL.

Instead:

> Use both, and give each information source a weight determined by its estimated uncertainty.

---

# 34. The role of \(S_{BL}\)

The BL MSE matrix is:

\[
S_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}.
\]

If \(S_{BL}\) becomes large, the BL estimate is less precise.

The paper therefore shows that a larger \(S_{BL}\) reduces the influence of the BL observation.

This is statistically sensible:

\[
\text{larger error variance}
\Rightarrow
\text{lower precision}
\Rightarrow
\text{lower weight}.
\]

This is one of the most important conceptual advantages of BL-MF over a simple arithmetic blend.

---

# 35. The role of \(D\)

The matrix:

\[
D
\]

contains the idiosyncratic residual variances.

If residual noise is high:

\[
D^{-1}
\]

becomes smaller.

Therefore the observed-return information contributes less to the factor estimate.

Again:

\[
\text{more noise}
\Rightarrow
\text{less precision}
\Rightarrow
\text{less weight}.
\]

---

# 36. MSE of the factor estimator

The paper gives:

\[
\boxed{
\widehat V_{f,BL-MF}
=
\left[
X'(\Sigma-D+S_{BL})^{-1}X
+
X'D^{-1}X
\right]^{-1}
}
\]

This matrix represents the estimated mean-squared-error structure of the fused factor estimator.

The paper compares it with the MSE from the factor model alone:

\[
\widehat V_f
=
(X'D^{-1}X)^{-1}.
\]

The fusion estimator is therefore designed to benefit from the additional BL information.

---

# 37. Step 6: estimate expected asset returns

Once the factor returns have been estimated over time, calculate the sample mean factor return:

\[
\widehat{\mu}_f
=
E(\widehat f_{BL-MF}).
\]

With the two-factor project specification:

\[
\widehat{\mu}_f
=
\begin{bmatrix}
\widehat E(MKT-R_f)\\
\widehat E(SMB)
\end{bmatrix}.
\]

Then:

\[
\boxed{
\widehat{\mu}
=
X\widehat{\mu}_f
}
\]

This produces the asset expected excess-return vector.

---

# 38. Step 7: estimate covariance under BL-MF

The paper estimates the covariance as:

\[
\boxed{
\widehat{\Sigma}_{BL-MF}
=
X\widehat{\Sigma}_fX'
+
\widehat{\Pi}
}
\]

where:

- \(\widehat{\Sigma}_f\) = sample covariance of the estimated factor returns;
- \(\widehat{\Pi}\) = diagonal covariance matrix of residuals.

This is the standard factor-model decomposition:

\[
\text{total covariance}
=
\text{systematic covariance}
+
\text{idiosyncratic covariance}.
\]

For two factors:

\[
X\widehat{\Sigma}_fX'
\]

is an \(N\times N\) systematic covariance component.

---

# 39. Important distinction: BL-MF produces both \(\mu\) and \(\Sigma\)

Although BL-MF is especially important for expected returns, the reference model ultimately produces both:

\[
\widehat{\mu}_{BL-MF}
\]

and:

\[
\widehat{\Sigma}_{BL-MF}.
\]

For the final integrated research model, however, the project design deliberately gives CV-MTSE responsibility for the covariance estimate.

Thus:

### Stand-alone BL-MF

\[
(\widehat{\mu}_{BL-MF},
\widehat{\Sigma}_{BL-MF})
\]

### Integrated model

\[
\boxed{
(\widehat{\mu}_{BL-MF},
\widehat{\Sigma}_{CV-MTSE})
}
\]

This distinction must remain explicit in the implementation.

---

# 40. Full BL-MF algorithm

At each estimation/rebalancing date:

## Step 1 — Prepare returns

Construct asset excess returns:

\[
r_t=R_t-R_{f,t}.
\]

## Step 2 — Construct the two factors

\[
MKT_t=R_{M,t}-R_{f,t}
\]

and:

\[
SMB_t
=
R_{small,t}-R_{large,t}.
\]

## Step 3 — Estimate factor loadings

Estimate:

\[
X_t.
\]

## Step 4 — Construct the BL prior

Specify:

\[
\pi
\]

and:

\[
C=\tau\Sigma.
\]

## Step 5 — Specify investor views

Specify:

\[
P,\quad q,\quad\Omega.
\]

## Step 6 — Compute BL posterior

\[
\widehat{\mu}_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}
(C^{-1}\pi+P'\Omega^{-1}q).
\]

## Step 7 — Compute BL uncertainty

\[
S_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}.
\]

## Step 8 — Estimate residual covariance

Construct:

\[
D.
\]

## Step 9 — Build fused observation system

\[
r_F=
\begin{bmatrix}
\widehat{\mu}_{BL}\\
r
\end{bmatrix},
\qquad
X_F=
\begin{bmatrix}
X\\
X
\end{bmatrix}.
\]

## Step 10 — Estimate factor returns

Use the weighted least-squares equation.

## Step 11 — Estimate factor mean and covariance

Calculate:

\[
\widehat{\mu}_f
\]

and:

\[
\widehat{\Sigma}_f.
\]

## Step 12 — Estimate asset expected returns

\[
\widehat{\mu}=X\widehat{\mu}_f.
\]

## Step 13 — Estimate BL-MF covariance

\[
\widehat{\Sigma}_{BL-MF}
=
X\widehat{\Sigma}_fX'
+
\widehat{\Pi}.
\]

## Step 14 — Pass the required outputs to portfolio optimisation

---

# 41. CV-MTSE and BL-MF side by side

| Feature | CV-MTSE | BL-MF |
|---|---|---|
| Primary problem | Covariance estimation | Expected-return/factor estimation |
| Core idea | Shrink noisy SCM toward structured targets | Fuse BL and factor information |
| Main data structures | SCM, SIM, IM | BL posterior, factor model, residual covariance |
| Data-driven parameter | \(\delta_{SIM},\delta_{IM}\) | Precision determined by \(S_{BL}\) and \(D\) |
| Selection method | Grid-search cross-validation | Weighted least squares |
| Main output | \(\widehat{\Sigma}_{CV-MTSE}\) | \(\widehat{\mu}_{BL-MF}\), and \(\widehat{\Sigma}_{BL-MF}\) |
| Project role | Risk/covariance model | Expected-return model |
| Integration | Supplies covariance | Supplies expected returns |

---

# 42. The integrated model

The final proposed model is not simply:

\[
CV-MTSE+BL-MF
\]

as an informal average.

It should be defined explicitly at the portfolio-input level.

The integrated inputs are:

\[
\boxed{
\widehat{\mu}_{INT}
=
\widehat{\mu}_{BL-MF}
}
\]

and:

\[
\boxed{
\widehat{\Sigma}_{INT}
=
\widehat{\Sigma}_{CV-MTSE}.
}
\]

The portfolio is then optimised using these two estimates.

Therefore:

\[
\boxed{
w^*
=
\arg\max_w
f
\left(
w;
\widehat{\mu}_{BL-MF},
\widehat{\Sigma}_{CV-MTSE}
\right)
}
\]

subject to the chosen portfolio constraints.

---

# 43. Portfolio optimisation

A mean-variance framework can be written generally as:

\[
\max_w
\left[
w'\widehat{\mu}
-
\frac{\gamma}{2}
w'\widehat{\Sigma}w
\right].
\]

Subject to:

\[
\mathbf{1}'w=1.
\]

Additional constraints must be specified before implementation, for example:

\[
w_i\geq0
\]

if short selling is prohibited.

The same optimisation objective and constraints should be applied across models wherever possible so that the comparison isolates the effect of the input estimation methods.

---

# 44. What changes across the models?

The optimisation machinery should remain as consistent as practical.

What changes is primarily the estimated input.

### EW

\[
w_i=\frac1N.
\]

### MVO

Uses the project's baseline estimates:

\[
\widehat{\mu}_{MVO},
\qquad
\widehat{\Sigma}_{MVO}.
\]

### CV-MTSE

Uses:

\[
\widehat{\mu}_{baseline},
\qquad
\widehat{\Sigma}_{CV-MTSE}.
\]

### BL-MF

Uses:

\[
\widehat{\mu}_{BL-MF},
\qquad
\widehat{\Sigma}_{BL-MF}.
\]

### Integrated

Uses:

\[
\widehat{\mu}_{BL-MF},
\qquad
\widehat{\Sigma}_{CV-MTSE}.
\]

This creates a clean research comparison.

---

# 45. Chronological implementation

The most important practical rule is:

> Every quantity used to construct a portfolio at time \(t\) must be computable using information available no later than \(t\).

For example, if the portfolio is formed at the end of December 2019:

- covariance estimation can use historical observations up to that date;
- BL views can use information available up to that date;
- SMB must be constructed using information available up to that date;
- the risk-free rate must be aligned correctly;
- no 2020 return can be used.

---

# 46. Training and testing design

The current project has:

\[
\boxed{\text{Training}=2010\text{--}2019}
\]

and:

\[
\boxed{\text{Out-of-sample}=2020\text{--}2025}.
\]

The training period is used for:

- model calibration;
- covariance estimation;
- factor estimation;
- expected-return estimation;
- BL parameters;
- shrinkage calibration;
- any other parameters required by the models.

The out-of-sample period is used for evaluating realised portfolio performance.

---

# 47. Rolling versus fixed estimation

There are two concepts that must not be confused.

### Fixed training model

Estimate everything once using 2010–2019 and hold the parameters fixed.

### Rolling/re-estimated model

At each rebalancing date, use only the historical data available up to that date and re-estimate the model.

The reference papers use rolling/rebalancing logic in their empirical portfolio construction.

For this project, the exact rolling-window/rebalancing convention should be declared explicitly in the methodology rather than being left implicit.

A defensible implementation is:

\[
\text{information available at }t
\rightarrow
\text{estimate}
\rightarrow
\text{portfolio at }t
\rightarrow
\text{hold}
\rightarrow
\text{rebalance}.
\]

---

# 48. Cross-validation must itself be inside the information set

A common implementation error is to perform CV-MTSE using the entire 2010–2025 sample.

That would leak future information.

The correct logic is:

\[
\text{Historical information available at }t
\rightarrow
\text{CV split}
\rightarrow
\text{choose }\delta_t^*
\rightarrow
\text{estimate covariance}
\rightarrow
\text{optimise portfolio}.
\]

The 2020–2025 observations cannot be used to choose a parameter for a portfolio formed before those observations occurred.

---

# 49. Factor construction must also be chronological

For MKT-RF:

\[
MKT_t-R_{f,t}.
\]

The risk-free rate must be transformed into the same return frequency as the asset and market returns.

For SMB, the portfolio sorts and weights must be based on market-capitalisation information available at the sorting date.

The implementation must document:

- the size breakpoint;
- portfolio weighting;
- rebalancing frequency;
- treatment of missing market capitalisation;
- treatment of firms entering/leaving the universe;
- whether the factor return is value-weighted or equal-weighted.

These are methodological choices, not coding details.

---

# 50. South African adaptation

The reference CV-MTSE paper uses Vietnamese equities.

The reference BL-MF paper uses Chinese equities and a substantially larger factor set.

The research project adapts them to South Africa.

Current specification:

### Assets

Satrix 40 / JSE Top 40 universe.

### Asset and market data

Yahoo Finance, subject to documented data-availability limitations.

### Risk-free rate

South African Reserve Bank 91-day Treasury Bill rate.

### Factors

\[
MKT-RF
\]

and:

\[
SMB.
\]

### Training

2010–2019.

### Out-of-sample

2020–2025.

### Integrated model

\[
\widehat{\mu}_{BL-MF}
\]

combined with:

\[
\widehat{\Sigma}_{CV-MTSE}.
\]

---

# 51. Important methodological distinction: source method versus adaptation

The following are directly grounded in the uploaded papers:

- MTSE combines SCM and multiple targets.
- CV-MTSE uses SIM and IM as targets.
- shrinkage parameters are selected by grid-search cross-validation.
- candidate shrinkage values are discretised in increments of 0.05.
- BL-MF treats the BL solution as a noisy observation.
- \(S_{BL}\) is used as its measurement-error covariance.
- BL and factor information are combined through weighted least squares.
- the factor model is used to estimate expected returns and covariance.

The following are adaptations for this research project:

- South African equities.
- Satrix 40.
- SARB 91-day Treasury Bill.
- two factors only.
- MKT-RF and SMB.
- 2010–2019 training period.
- 2020–2025 out-of-sample period.
- use of CV-MTSE covariance with BL-MF expected returns in the integrated model.

These should be clearly labelled as adaptations in the dissertation.

---

# 52. Practical Python architecture

The implementation should be modular.

A useful structure is:

```text
project/
│
├── data/
│
├── notebooks/
│   └── portfolio_research.ipynb
│
├── src/
│   ├── data_pipeline.py
│   ├── factors.py
│   ├── cv_mtse.py
│   ├── black_litterman.py
│   ├── bl_mf.py
│   ├── optimisation.py
│   ├── backtest.py
│   └── evaluation.py
│
└── results/
```

The notebook should orchestrate the research rather than contain one enormous block of code.

---

# 53. Suggested CV-MTSE functions

A transparent implementation can contain functions such as:

```python
def sample_covariance(returns):
    ...

def estimate_sim_target(asset_returns, market_returns):
    ...

def identity_target(n_assets, scale=None):
    ...

def generate_shrinkage_grid(step=0.05):
    ...

def mtse_covariance(sample_cov, sim_cov, identity_cov,
                    delta_sim, delta_identity):
    ...

def validation_loss(candidate_cov, validation_returns):
    ...

def cv_mtse(train_returns, validation_returns, market_returns):
    ...
```

The exact loss function should follow the paper's methodology rather than being replaced by an arbitrary portfolio-performance metric.

---

# 54. Suggested BL-MF functions

```python
def black_litterman_posterior(pi, C, P, q, Omega):
    ...

def black_litterman_mse(C, P, Omega):
    ...

def estimate_factor_loadings(asset_excess_returns, factor_returns):
    ...

def estimate_residual_covariance(residuals):
    ...

def bl_mf_factor_estimator(mu_bl, S_bl, X,
                            asset_returns, Sigma, D):
    ...

def bl_mf_expected_returns(X, factor_returns):
    ...

def bl_mf_covariance(X, factor_covariance, residual_covariance):
    ...
```

The functions should make the mathematical equations recognisable in the code.

---

# 55. Numerical stability

The models involve matrix inverses.

In practice, avoid naïvely calculating:

```python
np.linalg.inv(A)
```

where a linear solve is possible.

Prefer numerically stable approaches such as:

```python
np.linalg.solve(A, b)
```

or Cholesky-based solutions where appropriate.

Before optimisation, verify:

- symmetry;
- finite values;
- positive eigenvalues where positive definiteness is required;
- sensible condition numbers.

This is especially important for covariance matrices.

---

# 56. Essential diagnostic checks for CV-MTSE

For every estimation date, record:

\[
\delta_{SIM,t},
\qquad
\delta_{IM,t},
\qquad
\delta_{SCM,t}.
\]

Check:

\[
0\leq\delta_i\leq1
\]

and:

\[
\delta_{SIM,t}
+
\delta_{IM,t}
+
\delta_{SCM,t}
=1.
\]

Also record:

- validation loss;
- eigenvalues;
- condition number;
- covariance matrix symmetry;
- number of candidate grid points evaluated.

The time series of shrinkage weights is itself useful empirical evidence.

---

# 57. Essential diagnostic checks for BL-MF

At every estimation date verify:

### BL prior

\[
C=C'
\]

and positive definiteness where required.

### View uncertainty

\[
\Omega
\]

must be positive definite.

### BL MSE

\[
S_{BL}
=
S_{BL}'.
\]

### Residual covariance

\[
D
\]

should be diagonal and contain non-negative residual variances.

### Fused covariance

\[
\Sigma-D+S_{BL}
\]

must be suitable for the required inverse/linear solve.

### Factor loadings

Check dimensions:

\[
X\in\mathbb R^{N\times2}.
\]

### Factor covariance

Check:

\[
\widehat{\Sigma}_f
\in\mathbb R^{2\times2}.
\]

---

# 58. Common implementation errors

## Error 1 — Using future data in CV

Incorrect:

```text
2010–2025 → calculate covariance → cross-validation
```

for an early out-of-sample portfolio.

Correct:

```text
information available at t
→ training/validation split
→ choose shrinkage
→ estimate covariance
```

## Error 2 — Treating CV-MTSE as an expected-return model

It is primarily a covariance estimator.

## Error 3 — Averaging BL and factor returns

The BL-MF model is not a simple:

\[
0.5\widehat{\mu}_{BL}
+
0.5\widehat{\mu}_{MF}.
\]

It uses uncertainty-weighted statistical fusion.

## Error 4 — Ignoring \(S_{BL}\)

\(S_{BL}\) is central to the BL-MF fusion because it quantifies BL estimation uncertainty.

## Error 5 — Adding unnecessary factors

The current project uses:

\[
MKT-RF,\quad SMB.
\]

Do not silently add HML, MOM, RMW, or CMA.

## Error 6 — Mixing annual and periodic quantities

Returns, covariance, risk-free rates, and factor returns must be expressed consistently.

## Error 7 — Treating the identity matrix casually

Document its scale/normalisation.

## Error 8 — Optimising model parameters for maximum historical Sharpe ratio

Parameter selection must follow the methodology, not whichever configuration gives the best historical result.

---

# 59. Conceptual example of CV-MTSE

Suppose:

\[
\delta_{SIM}=0.25
\]

and:

\[
\delta_{IM}=0.35.
\]

Then:

\[
\delta_{SCM}=0.40.
\]

The final covariance is:

\[
\widehat{\Sigma}
=
0.40\widehat{\Sigma}_{SCM}
+
0.25\widehat{\Sigma}_{SIM}
+
0.35I.
\]

Interpretation:

- 40% weight on the observed covariance structure;
- 25% weight on market-factor structure;
- 35% weight on the stabilising identity structure.

These numbers are illustrative only. The actual values must come from the cross-validation procedure.

---

# 60. Conceptual example of BL-MF

Suppose the BL model produces:

\[
\widehat{\mu}_{BL}
\]

with relatively low:

\[
S_{BL}.
\]

Then the BL observation has high precision.

Suppose the factor model has residual covariance:

\[
D.
\]

If \(D\) is also small, both information sources are precise and both materially influence the factor estimate.

If \(S_{BL}\) becomes large, the BL information is down-weighted.

If \(D\) becomes large, the observed-return information is down-weighted.

Thus:

\[
\boxed{
\text{BL-MF is adaptive to information quality}
}
\]

rather than merely adaptive to return magnitude.

---

# 61. Conceptual example of the integrated model

Suppose BL-MF estimates:

\[
\widehat{\mu}_{BL-MF}
=
\begin{bmatrix}
0.08\\
0.05\\
0.03
\end{bmatrix}
\]

for three assets.

Suppose CV-MTSE estimates:

\[
\widehat{\Sigma}_{CV-MTSE}.
\]

The integrated optimiser receives exactly these two objects:

```text
expected_returns = mu_bl_mf
covariance = sigma_cv_mtse
```

and solves the same portfolio optimisation problem used for the comparison models.

The integrated model therefore tests a specific hypothesis:

> Better expected-return information from BL-MF and more stable covariance information from CV-MTSE may jointly produce a better portfolio than either estimation method alone.

The empirical results must determine whether this is true.

---

# 62. Research hypotheses that follow naturally

The methodology can support hypotheses such as:

### H1 — Covariance estimation

CV-MTSE produces a more stable covariance estimate than SCM.

### H2 — Expected-return estimation

BL-MF provides useful expected-return information by combining investor views and factor information.

### H3 — Integrated optimisation

The integrated BL-MF + CV-MTSE model improves out-of-sample portfolio performance relative to the individual model specifications.

These are empirical hypotheses, not assumptions that the models must outperform.

---

# 63. Backtesting sequence

At every rebalancing date:

```text
1. Identify assets available at t
        ↓
2. Gather historical information available at t
        ↓
3. Calculate asset returns
        ↓
4. Calculate market excess return
        ↓
5. Construct SMB
        ↓
6. Estimate baseline covariance
        ↓
7. Run CV-MTSE
        ↓
8. Construct BL prior and investor views
        ↓
9. Calculate BL posterior and S_BL
        ↓
10. Estimate two-factor loadings
        ↓
11. Run BL-MF fusion
        ↓
12. Obtain BL-MF expected returns
        ↓
13. Obtain CV-MTSE covariance
        ↓
14. Optimise integrated portfolio
        ↓
15. Hold portfolio
        ↓
16. Record realised return
        ↓
17. Rebalance at the next date
```

The same chronology should be used for the stand-alone model comparisons.

---

# 64. Performance evaluation

Do not rely on Sharpe ratio alone.

Useful measures include:

- annualised return;
- annualised volatility;
- Sharpe ratio;
- Sortino ratio;
- maximum drawdown;
- Calmar ratio;
- downside risk;
- turnover;
- transaction costs, if included;
- cumulative wealth.

All models must be evaluated on the same out-of-sample dates and under the same return/risk conventions.

---

# 65. Turnover

Portfolio turnover is important because a theoretically attractive portfolio can be impractical if it changes dramatically at every rebalance.

A simple turnover measure is:

\[
Turnover_t
=
\sum_i
|w_{i,t}^{new}-w_{i,t}^{old}|.
\]

If transaction costs are included, they should be applied consistently across all models.

---

# 66. Survivorship bias

The project uses a JSE/Satrix 40 universe and Yahoo Finance data.

Care must be taken because a dataset constructed only from securities that survive until the end of the sample can introduce survivorship bias.

If historical constituents are unavailable, the limitation should be disclosed.

Do not claim that the backtest is free from survivorship bias unless the constituent history supports that claim.

---

# 67. Why the methods complement each other

The methods attack different sources of estimation error.

CV-MTSE asks:

> How can the covariance matrix be made more stable?

BL-MF asks:

> How can information about expected returns be combined more intelligently?

Therefore the integrated model has a coherent economic/statistical rationale:

\[
\boxed{
\text{Risk estimation improvement}
+
\text{return estimation improvement}
}
\]

rather than an arbitrary combination of unrelated models.

---

# 68. Final mathematical summary

## CV-MTSE

\[
\boxed{
\widehat{\Sigma}_{CV-MTSE}
=
(1-\delta_{SIM}-\delta_{IM})
\widehat{\Sigma}_{SCM}
+
\delta_{SIM}\widehat{\Sigma}_{SIM}
+
\delta_{IM}I
}
\]

with:

\[
(\delta_{SIM}^*,\delta_{IM}^*)
=
\arg\min_{\delta}
L_{CV}(\delta).
\]

---

## Black–Litterman

\[
\boxed{
\widehat{\mu}_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}
(C^{-1}\pi+P'\Omega^{-1}q)
}
\]

and:

\[
\boxed{
S_{BL}
=
(P'\Omega^{-1}P+C^{-1})^{-1}.
}
\]

---

## Multi-factor model

\[
\boxed{
r=Xf+\epsilon.
}
\]

For this project:

\[
\boxed{
f=
\begin{bmatrix}
MKT-RF\\
SMB
\end{bmatrix}.
}
\]

---

## BL-MF fusion

\[
\boxed{
\widehat f_{BL-MF}
=
\left[
X'(\Sigma-D+S_{BL})^{-1}X
+
X'D^{-1}X
\right]^{-1}
\left[
X'(\Sigma-D+S_{BL})^{-1}\widehat{\mu}_{BL}
+
X'D^{-1}r
\right]
}
\]

then:

\[
\boxed{
\widehat{\mu}_{BL-MF}
=
X\widehat{\mu}_{f,BL-MF}.
}
\]

---

## Integrated model

\[
\boxed{
\widehat{\mu}_{INT}
=
\widehat{\mu}_{BL-MF}
}
\]

\[
\boxed{
\widehat{\Sigma}_{INT}
=
\widehat{\Sigma}_{CV-MTSE}
}
\]

and:

\[
\boxed{
w^*
=
\arg\max_w
f
\left(
w;
\widehat{\mu}_{BL-MF},
\widehat{\Sigma}_{CV-MTSE}
\right).
}
\]

---

# 69. The most important conceptual takeaway

The easiest way to remember the complete methodology is:

### CV-MTSE

\[
\boxed{
\text{Noisy SCM}
\rightarrow
\text{SCM + SIM + IM}
\rightarrow
\text{cross-validation chooses the weights}
\rightarrow
\widehat{\Sigma}_{CV-MTSE}
}
\]

### BL-MF

\[
\boxed{
\text{Equilibrium}
+
\text{Investor views}
\rightarrow
\widehat{\mu}_{BL}
+
S_{BL}
}
\]

then:

\[
\boxed{
\widehat{\mu}_{BL}
+
\text{observed asset returns}
+
\text{factor loadings}
+
\text{uncertainty matrices}
\rightarrow
\widehat f_{BL-MF}
\rightarrow
\widehat{\mu}_{BL-MF}
}
\]

### Integrated model

\[
\boxed{
\widehat{\mu}_{BL-MF}
+
\widehat{\Sigma}_{CV-MTSE}
\rightarrow
\text{portfolio weights}
\rightarrow
\text{out-of-sample evaluation}
}
\]

---

# 70. Implementation checklist

Before declaring the implementation complete, verify:

## Data

- [ ] Asset universe documented.
- [ ] Satrix 40 market series documented.
- [ ] SARB 91-day Treasury Bill documented.
- [ ] Return frequency documented.
- [ ] Missing observations handled explicitly.
- [ ] Corporate actions handled appropriately.

## CV-MTSE

- [ ] SCM implemented.
- [ ] SIM target implemented.
- [ ] Identity target implemented.
- [ ] Shrinkage constraints enforced.
- [ ] Grid uses the specified 0.05 increments.
- [ ] Cross-validation is chronological.
- [ ] Validation data are not used to fit training estimates.
- [ ] Selected shrinkage coefficients are recorded.
- [ ] Final covariance matrix is positive definite/suitable for optimisation.

## BL-MF

- [ ] BL prior specified.
- [ ] \(P\) specified.
- [ ] \(q\) specified.
- [ ] \(\Omega\) specified.
- [ ] \(C=\tau\Sigma\) documented.
- [ ] \(\widehat{\mu}_{BL}\) calculated.
- [ ] \(S_{BL}\) calculated.
- [ ] MKT-RF constructed.
- [ ] SMB construction documented.
- [ ] Two-factor loadings estimated.
- [ ] Residual covariance \(D\) estimated.
- [ ] BL-MF weighted least-squares estimator implemented.
- [ ] Factor means/covariances estimated.
- [ ] BL-MF expected returns calculated.

## Integrated model

- [ ] BL-MF supplies expected returns.
- [ ] CV-MTSE supplies covariance.
- [ ] Same optimiser and constraints used for comparable models.
- [ ] No future information enters either estimator.

## Backtest

- [ ] 2010–2019 training information respected.
- [ ] 2020–2025 remains out-of-sample.
- [ ] Rebalancing dates are defined.
- [ ] Portfolio returns are calculated consistently.
- [ ] Turnover is recorded.
- [ ] Performance metrics are calculated consistently.
- [ ] Limitations are disclosed.

---

# 71. Methodological caution for the dissertation

The uploaded papers should not be presented as if they were originally designed for this exact South African research problem.

A defensible dissertation should state that:

1. CV-MTSE is adapted from the Vietnamese-market application to South African equities.
2. BL-MF is adapted from the Chinese-market application.
3. The BL-MF factor dimension is reduced to two factors for this project.
4. The factor definitions are specifically constructed for the South African market.
5. The integrated model is the project's proposed combination of BL-MF expected-return estimation and CV-MTSE covariance estimation.
6. The empirical question is whether the adaptation improves out-of-sample performance, not whether improvement is assumed in advance.

---

# 72. Primary methodological sources

### CV-MTSE

Tran, M., Nguyen, N. M., & Tran, T. A. (2025). *Enhancing portfolio optimization in emerging markets: A cross-validation multi-target shrinkage approach*. Results in Control and Optimization, 21, 100611.

### BL-MF

Yuan, J., Jin, L., & Lan, F. (2025). *A BL-MF fusion model for portfolio optimization: Incorporating the Black–Litterman solution into multi-factor model*. Finance Research Letters, 80, 107464.

The uploaded papers are the primary methodological references for the implementation. Project-specific choices must be identified as adaptations rather than attributed to the papers.
