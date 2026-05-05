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

- **Transect-level confirmatory test for distance:** within-transect
  permutation of the presence labels (preserving the marginal number
  of detections per transect), recomputing the median distance among
  detected quadrants under each permutation; one-sided p = proportion
  of permutations with median ≤ observed.
  Result: observed median = 1.0 m, null median = 4.87 m
  (95% prediction interval [2.0, 7.5]); permutation p = 0.0015.
  This confirms the distance result is not an artefact of within-transect
  spatial autocorrelation.

- **Predictor collinearity:** Spearman ρ among all five predictors
  flags strong correlation between distance and vegetation cover
  (ρ = 0.60, p < 0.001), supporting interpretation as a single
  composite microhabitat axis ("exposed near-shore basalt within
  ~2 m of water") rather than two independent dimensions of selection.

### Statistical power

A Monte-Carlo simulation (n_sim = 20,000, normal-shift model
calibrated to the target rank-biserial r, seed = 42) estimates power
for the n₁ = 8 vs n₂ = 72 design at α = 0.05:

| Effect size (rank-biserial r) | Power |
|-------------------------------|-------|
| 0.30 (small) | 0.27 |
| 0.40 (moderate) | 0.47 |
| 0.50 (medium-large) | 0.69 |

Variables that did not survive Holm correction are therefore underpowered
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
