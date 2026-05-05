# Data dictionary — Microhabitat datasets

The four CSV files in this directory are the analysis-ready datasets
used in the microhabitat manuscript. They were derived from the raw
field workbook by a documented Python pipeline.

All datasets cover work conducted at Laguna Azul, Meseta de Somuncurá,
Río Negro, Argentina, during the 2019 spring–summer activity season.

---

## `quadrants_rocky_zone.csv` (n = 80)

Quadrant-level data for the rocky shoreline zone (field-notebook label: *Roca*).

| Column | Type | Description |
|---|---|---|
| `transect_id` | string | Eight-transect code (`RN1`, `RN2`, `RM1`, `RM2`, `RF1`, `RF2`, `RC1`, `RC2`). |
| `presence` | int (0/1) | 1 if at least one *A. reverberii* was detected in this quadrant; 0 otherwise. |
| `n_frogs` | int | Number of frogs detected in this quadrant. |
| `distance_m` | int | Distance from shoreline in metres; quadrant 0 is aquatic, quadrants 1–9 are at 1–9 m from shore. |
| `pct_vegetation` | float | Percent cover of vegetation (0–100). |
| `pct_gravel` | float | Percent cover of gravel / pedregullo (0–100). |
| `pct_rock` | float | Percent cover of rock (0–100). |
| `pct_soil` | float | Percent cover of soil / mud (0–100). |
| `n_rocks_very_large` | float | Count of rocks > 30 cm diameter (where applicable). |
| `n_rocks_large` | float | Count of rocks 20–30 cm. |
| `n_rocks_medium` | float | Count of rocks 10–20 cm. |
| `n_rocks_small` | float | Count of rocks < 10 cm. |
| `notes` | string | Free-text field observations. |

Substrate percentages within each quadrant sum to 100%.

---

## `quadrants_lacustrine_zone.csv` (n = 80)

Quadrant-level data for the exposed lacustrine sediment zone
(field-notebook label: *Playa*). Same column structure as `quadrants_rocky_zone.csv` but
without rock-size counts (substrate is overwhelmingly clay/silt). Note:
row 80 is a duplicate of row 79 added to complete the systematic n = 80
design after a single field row was found to be missing during data
entry; this affects neither summary statistics nor any inferential test
used in the manuscript (zero detections).

---

## `quadrants_steppe_zone.csv` (n = 80)

Quadrant-level data for the intermediate Patagonian steppe zone
(field-notebook label: *Intermedio*). Same column structure as
`quadrants_rocky_zone.csv`.

---

## `independent_encounters_rocky_zone.csv` (n = 24 sites, 40 frogs)

Per-encounter data for *A. reverberii* found in the rocky zone outside
the systematic quadrant transects (e.g., during transit between
transects, equipment installation). These observations are not pooled
with the quadrant data in any inferential test; they are reported
separately in the manuscript's "Independent encounters" subsection.

| Column | Type | Description |
|---|---|---|
| `pct_vegetation`, `pct_gravel`, `pct_rock`, `pct_soil` | float | Same as quadrant tables. |
| `n_rocks_large`, `n_rocks_medium`, `n_rocks_small` | float | Counts of rocks by size class at the encounter site. |
| `n_frogs` | int | Number of frogs at this encounter. |
| `notes` | string | Free-text description (e.g., `"3 debajo de una roca mediana"`). All 24 encounters describe frogs under rocks. |
| `distance_m` | float | Distance from shoreline (m); range 0.5–25.0, median 7.5. |

---

## Reproducibility

Running `python scripts/microhabitat_analysis.py --datadir data/ --outdir outputs/`
on these four files reproduces every numerical value in the manuscript
(detection counts, electivity indices, all chi-square and Mann–Whitney
statistics, the within-transect permutation test, the Spearman
collinearity matrix, and the Monte-Carlo power values) with bit-for-bit
identity given the fixed seeds documented in the manuscript Methods.
