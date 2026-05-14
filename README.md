# Microhabitat use by *Atelognathus reverberii* — Kass et al.

[![Python](https://img.shields.io/badge/python-3.9%2B-blue)](https://www.python.org/)

Reproducibility repository for:

> Kass, N.A., Kass, C.A., Tettamanti, G., Passamonti, L., Zarini, O.,
> Ipas, M.C., Kacoliris, F.P., and Williams, J.D.
> **Microhabitat use by *Atelognathus reverberii* at Laguna Azul,
> Patagonia: selection for the rocky shoreline–water ecotone.**
> *Ichthyology & Herpetology*, sent for review.

All figures, tables, and statistics reported in the paper are produced
by running `scripts/microhabitat_analysis.py` on the four datasets in
`data/`. Results are bit-for-bit identical to the published values
given the fixed random seeds documented in the manuscript Methods.

---

## Repository structure

```
atelognathus-reverberii-microhabitat/
│
├── data/
│   ├── quadrants_rocky_zone.csv         ← rocky shoreline (n = 80 quadrants)
│   ├── quadrants_lacustrine_zone.csv    ← exposed lacustrine sediment (n = 80)
│   ├── quadrants_steppe_zone.csv        ← intermediate Patagonian steppe (n = 80)
│   ├── independent_encounters_rocky_zone.csv  ← incidental encounters (n = 24, 40 frogs)
│   └── README_data.md                   ← full column descriptions
│
├── scripts/
│   └── microhabitat_analysis.py         ← all analyses, figures, tables
│
├── figures/                             ← publication figures (300 DPI PNG)
│   ├── fig1_zone_patterns.png
│   └── fig2_within_zone.png
│
├── outputs/                             ← regenerated when the script runs
│
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Installation

### Requirements

- Python 3.9 or later
- pip

### Step-by-step setup

**1. Clone the repository**

```bash
git clone https://github.com/nicolaskass/atelognathus-reverberii-microhabitat
cd atelognathus-reverberii-microhabitat
```

**2. Create a virtual environment** (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate        # Linux / macOS
# .venv\Scripts\activate         # Windows
```

**3. Install dependencies**

```bash
pip install -r requirements.txt
```

| Package | Version | Purpose |
|---------|---------|---------|
| `pandas` | ≥2.0 | Data loading and tabular operations |
| `numpy` | ≥1.24 | Numerical computations, bootstrap, permutation |
| `matplotlib` | ≥3.7 | All publication figures |
| `scipy` | ≥1.11 | Mann–Whitney U, Pearson chi-square, Fisher's exact, Spearman, Beta distribution (Jeffreys bound) |

---

## Running the analysis

### Full pipeline

```bash
python scripts/microhabitat_analysis.py \
    --datadir data/ \
    --outdir outputs/
```

This produces in `outputs/`:

- `fig1_zone_patterns.png` — Figure 1 of the paper (3 panels: detection rate, Manly's α, substrate cover)
- `fig2_within_zone.png` — Figure 2 of the paper (3 panels: detection by distance, used vs available substrate, IE distance)
- `within_zone_mw_results.csv` — Table 2 of the paper (Mann–Whitney + Holm)
- `within_zone_transect_summary.csv` — per-transect detection positions (n = 5 positive transects)
- `within_zone_collinearity.csv` — full Spearman matrix among all 5 within-zone predictors
- `reportable_quantities.json` — single source of every numerical value reported in the paper

And prints a full summary to stdout including all statistics reported in the paper.

### Parameters explained

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--datadir` | `data/` | Directory containing the four CSV datasets. |
| `--outdir` | `outputs/` | Directory where figures, tables, and JSON are saved. |

---

## Analysis details

### Two-level test for between-zone heterogeneity

Quadrants are nested within transects (10 contiguous quadrants per
transect, 8 transects per zone), so treating individual quadrants as
independent units would be pseudoreplicated. We therefore conduct
between-zone tests at two levels:

- **Quadrant level** (descriptive pooled summary): Pearson χ² on the
  2 × 3 detected/not-detected by zone contingency table.
  Result: χ² = 10.85, p = 0.004.
- **Transect level** (confirmatory inferential test): each transect
  scored as positive (≥ 1 detection) or negative; Pearson χ² on the
  transect-level 2 × 3 table; verified with a permutation test
  (n_perm = 50,000, seed = 42).
  Result: χ² = 7.66, p = 0.022 (asymptotic), p = 0.033 (permutation).

### Habitat electivity (Manly's α)

Manly's standardised selection-ratio (Manly et al. 1993, eq. 4.10):

$$\alpha_i = \frac{o_i / n_i}{\sum_j (o_j / n_j)}$$

Random-use threshold for K = 3 zones is 1/K = 0.333. Bootstrap 95% CIs
use 10,000 resamples (seed = 42). For zones with zero detections the
bootstrap CI degenerates to a point; we replace it with the upper 97.5%
Jeffreys credible bound on the underlying detection probability under a
Beta(0.5, 80.5) posterior.

Result for the lacustrine sediment zone: α = 0.000, upper Jeffreys bound ≈ 0.198.

### Within-zone use vs availability

Substrate cover was recorded for **all 80 rocky-zone quadrants**
regardless of detection outcome, enabling a formal use-vs-availability
analysis with full availability information rather than a
presence-only proxy.

- **Quadrant-level exploratory tests:** Mann–Whitney U comparisons of
  each of five microhabitat variables (distance, vegetation, gravel,
  rock, soil) between the 8 occupied and 72 unoccupied quadrants;
  Holm's sequential Bonferroni correction across the five tests.
  Result: distance to shore (p = 0.002, surviving Holm) and
  vegetation cover (p = 0.012, surviving Holm) are the two
  predictors that pass.

- **Transect-level test for distance:** within-transect permutation of
  the presence labels (preserving the marginal number of detections
  per transect), recomputing the median distance among detected
  quadrants under each permutation; one-sided p = proportion of
  permutations with median ≤ observed.
  Result: observed median = 1.0 m, null median = 4.87 m
  (95% prediction interval [2.0, 7.5]); permutation p = 0.0015.
  Leave-one-positive-transect-out worst-case: p ≤ 0.018 — robust to
  any single transect being dropped.

- **Predictor collinearity:** Spearman ρ among all five predictors
  flags strong correlation between distance and vegetation cover
  (ρ = 0.60, p < 0.001), supporting interpretation as a single
  composite microhabitat axis ("exposed near-shore basalt within
  ~2 m of water") rather than two independent dimensions.

### Pooled-adjacent transect sensitivity (spatial autocorrelation)

Eight transects per zone were deployed sequentially along one continuous
lagoon perimeter, raising the possibility of inter-transect spatial
autocorrelation. We pooled consecutive transect pairs into n=4
block-transects per zone and re-ran the between-zone test.
Result: χ² = 4.80, p = 0.091, permutation p = 0.21. The rocky-vs-
lacustrine contrast (3/4 vs 0/4) survives the aggregation; the rocky-
vs-intermediate contrast does not. The manuscript therefore frames
the between-zone result as **consistent with a strong rocky/lacustrine
contrast in detection rate**, not as a confirmed three-zone heterogeneity
test.

### Substrate-conditioned detectability sensitivity (3 models)

Under what relative detection probability between zones would the
observed lacustrine zero be compatible with equal underlying use? We
report three parameterisations:

| Model | Assumption | E[lacustrine \| equal use] | P(observe 0) |
|-------|------------|----------------------------|--------------|
| Linear | p ∝ rock cover | 0.47 events | 0.62 |
| Step | p = 0 if rock < 5% | 0.00 events | 1.00 |
| Saturating | Hill (K = rock_rocky / 4) | 1.91 events | 0.15 |

All three converge on the same conclusion: the lacustrine zero is
plausible under equal underlying use given the rock-turning protocol's
substrate-conditioned detectability. The qualitative claim that
surface emergence is concentrated in the rocky zone is robust to this
parameterisation; the quantitative magnitude is not.

### Manly's α from quadrant counts vs event counts

To verify that the event-based numerator is not inflated by within-
quadrant aggregation (the rocky zone has ~6.5 events per occupied
quadrant; the intermediate zone ~1.5), we recomputed Manly's α from
**occupied-quadrant counts** (8/2/0) instead of detection-event counts
(52/3/0). The two formulations give identical α values
(rocky=0.800, intermediate=0.200, lacustrine=0.000) because the
standardised selection ratio is determined by between-zone rates.

### Aquatic-quadrant exclusion sensitivity

Excluding the 3 aquatic quadrants (distance=0, structurally false-
negative under turbid water) from the rocky-zone denominator yields a
detection rate of 8/77 = 10.4% vs the reported 8/80 = 10.0% — no
qualitative change in any test.

### Independent encounters (IE) vs quadrant distance distribution

Mann-Whitney one-sided test of the hypothesis that IE distances exceed
quadrant-detected distances: median IE = 7.5 m, median quadrant = 1 m,
p = 0.004. The IE distribution extends to 25 m, well beyond the
9-m terrestrial reach of the systematic transects.

### Rock-turnability gradient and cavity-density dissociation

Initially hypothesised that rocks at the shoreline edge are more
"turnable" (smaller, looser) than at the rocky-steppe interior,
biasing protocol detection toward the shore. The data refute this on
both proportion and absolute count: the proportion of small+medium
movable rocks per quadrant position is 0.44 at d=0 m, 0.52 at d=1 m,
0.71 at d=2 m (increasing inward), and the **absolute number of
movable rocks per quadrant** is 5.8 at d=1 m versus 12.5 at d=3 m and
12.6 at d=5 m. The systematic 1-m detection peak therefore cannot be
attributed to either a movable-rock-fraction gradient or a movable-
rock-count gradient — at the distances where 5 of 8 detections occur,
there are fewer turnable shelters available than at distances where
zero detections occurred. This positively rules out the simplest
cavity-density artefact.

### Within-zone Mann-Whitney with vs. without aquatic quadrants

The 3 aquatic quadrants (distance = 0, positioned 1 m into the lagoon)
are structural false-negatives (protocol cannot detect submerged frogs
through turbid water). Re-running the within-zone MW for distance with
those 3 excluded yields **p = 0.001** versus the reported p = 0.002
with them retained. Inclusion of the aquatic quadrants is therefore
conservative; the contrast strengthens when they are excluded.

### Statistical power

A Monte-Carlo simulation (n_sim = 20,000, normal-shift model
calibrated to the target rank-biserial r, seed = 42) estimates power
for the n₁ = 8 vs n₂ = 72 design at α = 0.05 (uncorrected, no Holm):

| Effect size (rank-biserial r) | Power |
|-------------------------------|-------|
| 0.30 (small) | 0.27 |
| 0.40 (moderate) | 0.47 |
| 0.50 (medium-large) | 0.69 |

These are uncorrected upper bounds; Holm-corrected power for
non-leading tests in the family is correspondingly lower. Variables
that did not survive Holm correction are therefore underpowered
rather than null; this is reported as a Limitation of the manuscript.

---

## Reproducing exact paper values

All random operations use fixed seeds (42 unless otherwise stated).
Running `microhabitat_analysis.py` on the provided datasets reproduces
every value in the paper exactly:

| Statistic | Paper | Script output |
|-----------|-------|---------------|
| Rocky-zone detections | 8/80 | ✓ |
| Intermediate-zone detections | 2/80 | ✓ |
| Lacustrine-zone detections | 0/80 | ✓ |
| α (rocky) | 0.800 | ✓ |
| Quadrant χ² | 10.85 | ✓ |
| Transect χ² | 7.66 | ✓ |
| Permutation p (transect) | 0.033 | ✓ |
| Mann–Whitney p (distance) | 0.002 | ✓ |
| Mann–Whitney p (vegetation) | 0.012 | ✓ |
| Within-transect permutation p (distance) | 0.0015 | ✓ |
| Spearman ρ (distance vs vegetation) | 0.60 | ✓ |
| Power for r = 0.40 | 0.47 | ✓ |

---

## Citing this work

If you use this code or dataset, please cite:

```bibtex
@article{kass2026microhabitat,
  author  = {Kass, Nicol{\'a}s Ariel and Kass, Camila Alejandra and
             Tettamanti, Germ{\'a}n and Passamonti, Leandro and
             Zarini, Ornela and Ipas, M{\'o}nica Cecilia and
             Kacoliris, Federico Pablo and Williams, Jorge Daniel},
  title   = {Microhabitat use by \textit{Atelognathus reverberii} at
             Laguna Azul, {Patagonia}: selection for the rocky
             shoreline--water ecotone},
  journal = {Ichthyology \& Herpetology},
  year    = {2026},
  note    = {Sent}
}
```

---

## Contact

Nicolás Ariel Kass · nicolaskass@gmail.com
Sección Herpetología, División Zoología Vertebrados
Facultad de Ciencias Naturales y Museo, UNLP · Paseo del Bosque s/n,
La Plata, Argentina
