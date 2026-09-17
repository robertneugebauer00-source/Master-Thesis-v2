# Single-site chemical co-occurrence networks — complete methods breakdown

**LIONESS, BONOBO and lioness_glasso as implemented in `R/05_lioness.R` (pipeline steps S6–S8)**
Robert Neugebauer, MSc thesis — sgh-chemical-networks. 16.09.2026.

---

## 1. Data and preprocessing

The input to every estimator is the same matrix the mainline glasso sees:

- $X$: $n$ sites × $p$ compounds, built per spatial unit (basin, country group, river section).
- Node set: the **global pool** — compounds detected at ≥ 10 % of all reference sites in the analysed basins; within a unit, a compound additionally needs ≥ 2 detections (zero-variance guard).
- Non-detects (MEC = 0) are replaced by a compound-specific **MDL/2 floor** (fallback 1 ng/L); detected values are used exactly as measured.
- $X_{\log} = \log_{10}$ of the floored matrix. All estimators below work on $X_{\log}$.

Chemical meaning: a correlation on $X_{\log}$ compares *relative* concentration patterns across sites, not absolute levels — two chemicals co-vary if they rise and fall together along the river, which is the signature of a **shared source** (same WWTP plume, same agricultural season).

## 2. The aggregate network (the reference every site is compared against)

**Correlation scale.** $S = \mathrm{cor}(X_{\log})$ — the $p \times p$ sample correlation matrix over the unit's sites.

**Graphical lasso (glasso).** The mainline estimates a sparse precision matrix $\Theta$ by solving

$$\hat{\Theta} = \arg\min_{\Theta \succ 0}\; \mathrm{tr}(S\Theta) - \log\det\Theta + \lambda \sum_{i \neq j} |\Theta_{ij}|$$

and converts it to **partial correlations**:

$$r_{ij} = -\frac{\hat{\Theta}_{ij}}{\sqrt{\hat{\Theta}_{ii}\,\hat{\Theta}_{jj}}}, \qquad r_{ii} = 0 .$$

A non-zero $r_{ij}$ means chemicals $i$ and $j$ co-vary *conditionally on all other measured compounds* — a residual co-source link that the rest of the mixture cannot explain.

**StARS $\lambda$ selection.** The penalty is chosen by Stochastic Approach to Regularization Selection: draw $B = 100$ subsamples of size $b = \lfloor 10\sqrt{n} \rfloor$ (for $n > 144$; else $\lfloor 0.8n \rfloor$), refit the glasso on a grid of 30 log-spaced penalties from $\lambda_{\max} = \max_{i \neq j}|S_{ij}|$ down to $0.1\,\lambda_{\max}$, and measure edge instability $\bar{\xi}_{ij}(1-\bar{\xi}_{ij})$ (mean occurrence frequency $\bar{\xi}$ across subsamples), averaged over edges and scaled by $4/(p(p-1))$. The selected $\lambda^*$ is the sparsest grid point whose total instability stays below 0.05 (one step back from the first crossing). The support of $\hat\Theta$ at $\lambda^*$ is the **skeleton** — the fixed edge set all per-site scores are read on.

## 3. LIONESS — linear leave-one-out extrapolation

**Reference:** Kuijjer et al. (2019), *iScience* 11:226–243.

**Idea.** You cannot compute a correlation from one site. Instead, ask: *how does the aggregate estimate change because this site exists?* Let $G$ be any aggregate network statistic computed from all $n$ sites, and $G^{(-q)}$ the same statistic computed without site $q$. The LIONESS single-sample network for site $q$ is

$$\boxed{\;e_q = n \cdot G - (n-1) \cdot G^{(-q)}\;}$$

applied element-wise to every edge. In the pipeline, $G$ is the sample correlation matrix $S$ (`lioness_cor`).

**Why this form.** It is the jackknife/pseudo-value construction: if the aggregate were exactly the mean of the per-site networks, $G = \frac{1}{n}\sum_q e_q$, then solving for $e_q$ gives precisely this formula. So LIONESS *defines* the per-site network as the pseudo-value that makes the decomposition exact. For the sample correlation the identity holds up to an $O(1/n)$ jackknife bias, i.e. $\frac{1}{n}\sum_q e_q \approx S$.

**Properties.**

| property | consequence |
|---|---|
| linear extrapolation | entries are **not bounded** to $[-1, 1]$ (values like $+1.08$ occur) |
| $e_q$ need not be positive semidefinite | not a valid correlation matrix |
| leave-one-out differences are amplified by factor $n$ | **leverage inflation**: sites with extreme multivariate profiles (very polluted, or chemically unusual) get extreme weights on many edges |
| exactness | $\frac{1}{n}\sum_q e_q \approx G$; no tuning parameters |

**Chemical reading.** $e_{q,ij} > G_{ij}$: at site $q$ the pair $(i,j)$ co-occurs *more* strongly than the basin average — the site reinforces that source relationship. $e_{q,ij} < 0$ on a positive aggregate edge: this site's profile actively *works against* the usual co-occurrence (e.g. $i$ high, $j$ absent). Because of the $n$-fold amplification, a single extreme cocktail (e.g. EUS_029, sum-TU = 276) is extrapolated into extreme edge weights — this is the mechanism behind LIONESS's spurious *positive* strength–stress coupling in S7.

## 4. BONOBO — Bayesian single-sample correlation networks

**Reference:** Saha, Fanfani et al. (2024), *Genome Research*, doi:10.1101/gr.279117.124; implementation calibrated to netZooPy `compute_bonobo`.

**Idea.** Model each site's covariance as a Bayesian posterior mean between what the *other* sites say (prior) and what only this site says (likelihood):

$$\boxed{\;\Sigma_q = \delta_q \, dx_q\, dx_q^{\top} + (1-\delta_q)\, S^{(-q)}\;}$$

- $S^{(-q)} = \mathrm{cov}(X_{\log}^{(-q)})$: the leave-one-out sample covariance (ddof = 1) — the **prior**, the bulk correlation structure without site $q$.
- $dx_q = x_q - \bar{x}$: site $q$'s row minus the **full-data** column means; $dx_q dx_q^{\top}$ is the site's own outer-product deviation — the **likelihood**, what only this site contributes.
- $\delta_q \in (0,1)$: the posterior weight, tuned from the data (netZooPy form):

$$\delta_q = \frac{1}{3 + 2\,\dfrac{\mathrm{mean}\big(\sqrt{\mathrm{diag}(S^{(-q)})}\big)}{\mathrm{var}\big(\mathrm{diag}(S^{(-q)})\big)}}$$

with $\mathrm{var}$ the **population** variance (ddof = 0) over the $p$ diagonal entries; flat fallback $\delta_q = 1/3$ if the variance is 0. Heterogeneous compound variances (large var) → smaller $\delta_q$ → more shrinkage to the prior.

**Scaling to a correlation.**

$$R_q = D_q^{-1/2}\,\Sigma_q\,D_q^{-1/2}, \qquad D_q = \mathrm{diag}(\Sigma_q)$$

(zero/NA diagonal entries replaced by 1). Since $\Sigma_q$ is a convex combination of two PSD matrices it is PSD, so **$R_q$ is a valid correlation matrix: PSD, diagonal 1, all entries in $[-1, 1]$**.

**Optional sparsification** (computed in S6 with $\alpha = 0.05$; the S6/S7 summaries use the *dense* nets). The approximate sampling variance of the entries under the posterior, with $g = p$ and $d = g + 1/\delta_q$:

$$a_1 = \frac{d-g+1}{(d-g)(d-g-3)}, \qquad a_2 = \frac{d-g-1}{(d-g)(d-g-3)}$$
$$\mathrm{sd}_{jk} = \sqrt{a_1\,\Sigma_{q,jk}^2 + a_2\,\Sigma_{q,jj}\,\Sigma_{q,kk}}, \qquad z_{jk} = \frac{\Sigma_{q,jk}}{\mathrm{sd}_{jk}}$$

keep entries with $|z_{jk}| > \Phi^{-1}(1-\alpha/2)$ (two-sided normal test; diagonal kept at 1).

**Chemical reading.** $R_{q,ij}$ reads like an ordinary correlation: near $+1$ = this site strongly supports the co-occurrence of $i$ and $j$; near $0$ or negative = the site does not express that source relationship, or expresses it inversely. The shrinkage prior guarantees that an extreme site is pulled toward the leave-one-out bulk instead of being extrapolated — the exact opposite of LIONESS's leverage inflation. This is why BONOBO is the defensible estimator for any score that assumes a valid correlation matrix.

## 5. lioness_glasso — the partial-correlation backbone (S8)

The S8 robustness variant applies the LIONESS extrapolation to **glasso partial correlations** instead of raw correlations:

$$e_q = n \cdot A_{\mathrm{agg}} - (n-1) \cdot A^{(-q)}$$

where $A_{\mathrm{agg}}$ and $A^{(-q)}$ are glasso partial-correlation matrices fitted at the **fixed** aggregate StARS $\lambda^*$ (no per-site re-selection), and the result is **restricted to the aggregate skeleton** ($e_{q,ij} = 0$ wherever $A_{\mathrm{agg},ij} = 0$). It changes the association measure (conditional instead of marginal co-occurrence) while keeping edges comparable across sites. Cost: $n$ glasso solves per unit. Numerical note: leave-one-out fits can diverge under coordinate descent when $p > n$ with near-perfectly correlated compound pairs (chemical near-duplicates, e.g. co-emitted tracers); the delivered S8 run solves those with a validated ADMM fallback (99.98 % edge agreement with the CD solver where both converge).

## 6. From per-site networks to site scores

All downstream tests use per-site summaries **on the aggregate StARS skeleton** (same edges for every site — the fixed-node-set idea, one level down). For site $q$ with per-site network $W_q$ and skeleton edge set $\mathcal{E}$:

$$\text{strength}_q = \sum_{(i,j) \in \mathcal{E}} |W_{q,ij}|, \qquad \overline{w}_q = \frac{1}{|\mathcal{E}|} \sum_{(i,j) \in \mathcal{E}} W_{q,ij}$$

- **strength**: how strongly the site's whole mixture participates in the basin's co-occurrence structure. High = the site reinforces many edges (chemically "typical"); low = the site's profile is discordant with the usual source structure (chemically "atypical/idiosyncratic").
- These scores feed the site-level stress-gradient tests (S7: Spearman ρ vs sum-TU per unit, pooled by Stouffer's method with a within-unit permutation null; external proxy UDF; trivial-score benchmark; biomarker F-test) and the S8 backbone sensitivity.

## 7. Worked example (real numbers, Elbe, site EUS_434, edge Carbamazepine–Sucralose)

Both are wastewater markers; their aggregate correlation is strong.

| quantity | value |
|---|---|
| $n$, $p$ | 153 sites, 182 compounds |
| full-data correlation $S_{ij}$ | $+0.7904$ |
| leave-one-out correlation $S^{(-q)}_{ij}$ | $+0.7886$ |
| **LIONESS** $e_{q,ij} = 153(0.7904) - 152(0.7886)$ | $\mathbf{+1.0784}$ — outside $[-1,1]$ |
| $\mathrm{var}(\mathrm{diag}\,S^{(-q)})$ (ddof 0), $\mathrm{mean}\sqrt{\mathrm{diag}}$ | $0.7556$, $0.8946$ |
| $\delta_q = 1/(3 + 2 \cdot 0.8946/0.7556)$ | $0.1863$ |
| $dx_{q,i}$, $dx_{q,j}$ (log10 units) | $+1.309$, $+1.123$ |
| $\Sigma_{q,ij} = 0.1863(1.309 \cdot 1.123) + 0.8137(0.8542)$ | $+0.9687$ |
| **BONOBO** $R_{q,ij}$ (scaled) | $\mathbf{+0.8387}$ — a valid correlation |

Same site, same edge: both methods say EUS_434 reinforces the Carbamazepine–Sucralose wastewater edge more than the basin average (aggregate $+0.79$), but LIONESS extrapolates to an impossible correlation ($+1.08$) while BONOBO returns a bounded, interpretable value ($+0.84$).

## 8. Summary

| | LIONESS | BONOBO | lioness_glasso |
|---|---|---|---|
| equation | $nG - (n-1)G^{(-q)}$ | $\delta_q\, dx_q dx_q^\top + (1-\delta_q) S^{(-q)}$ | $nA_{\mathrm{agg}} - (n-1)A^{(-q)}$ |
| input | correlations | covariances | glasso partial correlations at fixed $\lambda^*$ |
| valid correlation matrix? | no (not PSD, unbounded) | **yes** (PSD, $[-1,1]$) | no (unbounded) |
| behaviour at extreme sites | **inflates** (leverage × $n$) | **shrinks** toward LOO bulk | inflates; on skeleton ≈ detection-richness score |
| tuning | none | $\delta_q$ from data | $\lambda^*$ from aggregate StARS |
| S7/S8 verdict | positive coupling, vanishes vs UDF | **negative coupling, survives UDF + trivial-score controls** | positive coupling, collapses vs trivial scores |

**References.** Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K (2019) Estimating sample-specific regulatory networks. *iScience* 11:226–243. · Saha A, Fanfani V, et al. (2024) BONOBO: Bayesian single-sample networks. *Genome Research* doi:10.1101/gr.279117.124 (netZooPy implementation). · Liu H, Roeder K, Wasserman L (2010) StARS stability selection. · Friedman J, Hastie T, Tibshirani R (2008) Sparse inverse covariance estimation with the graphical lasso. *Biostatistics* 9:432–441.
