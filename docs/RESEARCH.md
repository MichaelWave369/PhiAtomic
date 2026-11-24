# Research protocol

## Hypotheses and falsification

H0: recursive features provide no reliable out-of-sample improvement over standard composition and shell features. H1: they provide reproducible incremental information. Binding energy per nucleon is the first target, not half-life or stability. The literal social-media claim is already inconsistent with bounded positive-real Mandelbrot behavior: Fe-56 has c=26/56 > 1/4 and escapes at iteration 5 from w0=0. Ni-62 has the larger binding energy per nucleon; neither fact establishes a fractal mechanism.

## Mappings

Define A=Z+N and f=Z/A. Iterate w(t+1)=w(t)^2+c from w0=0, with escape radius 2 and cap 64. Compare c=f, c=f+i(N-Z)/A, and c=2f-1+i(N-Z)/A. The latter two are heuristics. They are deterministic functions of f and do not add independent composition information. All three collide for equal proton fraction; finite iteration non-escape is not a mathematical proof of boundedness. No physical time is represented.

## Baseline and benchmark

Liquid-drop features: volume, surface, Coulomb, asymmetry, pairing, expressed per nucleon; coefficients are fitted rather than fixed to a published parameterization. Shell proxy: smooth proximity to known magic numbers, not a microscopic shell calculation. Recursive features: capped escape iteration, log magnitude, escape indicator for each mapping. Fit ridge regression with fixed alpha 0.1 and train-fold-only feature scaling. Split by proton number with seeded group assignment. Report pooled and fold MAE. The demo is small, selected, and includes hydrogen, which creates unstable extrapolation when held out; the results are not evidence for H1.

## Next validation gates

1. Import a versioned evaluated ground-state mass table with provenance and uncertainty. Separate extrapolations from observations. Audit units and duplicates.
2. Add raw ratio polynomial/spline and matched-capacity nonlinear composition baselines. Since the recursion is only a transformation of f, demonstrate improvement against flexible nonrecursive transformations, not merely a linear baseline.
3. Freeze mappings, model capacity, splits, metrics and success thresholds before final testing. Nested grouped cross-validation for tuning; independent evaluation or contiguous holdout region. Compare paired errors and uncertainty by group. Include negative controls and multiple-testing correction.
4. Explore shell effects and topology on residuals after baseline subtraction. Only claim incremental predictive performance if it survives the controls. Mechanistic interpretation requires independent physics and validation.
5. A future Altermath/Nested Bubble Gear adapter could represent states and relationships (Z,N,R) as typed graph data. R denotes model interactions, not a new particle. Keep the adapter independent of the benchmark; require falsifiable predictions.

## Sources

- AME2020 II, Wang et al., Chinese Physics C 45 030003 (2021), https://doi.org/10.1088/1674-1137/abddaf
- IAEA LiveChart API notebook: https://iaea-nds.github.io/lc_api_notebook/
- Argonne mass data resources: https://www.anl.gov/phy/atomic-mass-data-resources
