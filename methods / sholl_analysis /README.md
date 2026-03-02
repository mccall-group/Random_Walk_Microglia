# Sholl-style multi-radius trace-back (BLUE) — README

This folder documents the **multi-radius trace-back** workflow (often used like a Sholl-style validation) built on a 2D probability map `prob2d` and a fixed soma coordinate `(x0, y0)`.  
It generates **BLUE random-walk paths** from boundary seeds at multiple radii and provides **cross-round consistency filtering** based on seed angles.

---

## 1. Inputs and expected objects

### Required inputs
- `prob2d`: 2D matrix (H × W). Higher values indicate more likely/process-supporting pixels.
- `x0, y0`: soma center coordinates in the same pixel coordinate system as `prob2d`.

### Optional upstream (how `prob2d` is typically obtained)
- `prob2d <- so@detected_signal_stack[,, z]` from MGtoolbox after `detected_territory()`.

---

## 2. Core functions and what each block does

### Block A — `rw_path_2d_multi_to_center()` (BLUE trace-back)

**Purpose:** Trace a path from a boundary seed `(x_start, y_start)` toward soma `(x0, y0)` using a biased random-walk.

**Per-step logic:**
1. Sample `n_cand` candidate moves `(dx, dy)` within `[-step_size, step_size]`, excluding `(0,0)`.
2. Compute candidate locations `(xx, yy) = (x[t]+dx, y[t]+dy)`.
3. Filter candidates to:
   - stay inside the image bounds,
   - stay inside the circle of radius `radius` around soma,
   - enforce monotone-ish progress toward soma: `d_new^2 <= d_cur^2 + 3`.
4. Compute weights for remaining candidates:
   - `w = (prob2d[yy,xx] + eps)^alpha`.
5. Sample next step from candidates with probability proportional to `w`.
6. Stop early if the path reaches within ~2 pixels of soma.

**Key parameters:**
- `alpha` (bias strength): larger → stronger preference for high-`prob2d` pixels.
- `step_size`: local move radius per step.
- `n_cand`: number of proposals per step (higher → more robust but slower).
- `radius`: defines the active circle (“round”) for this trace-back run.

---

### Block B — helpers for masking and coverage

#### `circle_mask(H, W, x0, y0, radius)`
Creates a boolean mask for pixels inside the circle centered at soma.

#### `build_coverage(paths_list, H, W)`
Aggregates all visited pixels across a list of paths into an integer coverage map.

Used in the repo mainly for:
- diagnostics,
- later modules (forbid/residual logic),
but coverage is **not required** for the BLUE-only Sholl-style rounds.

---

### Block C — boundary seeds + utilities

#### `get_boundary_points(prob2d, x0, y0, radius, n_theta)`
Samples points along the circle perimeter by sweeping `theta ∈ [0, 2π)`:
- Generates integer boundary coordinates `(xs, ys)`.
- Looks up `prob2d[ys, xs]` as `val`.
- Returns unique boundary points.

#### Seed selection by separation
From boundary candidates sorted by `val` decreasing, select up to `K_boundary` seeds such that pairwise distances are at least `min_sep_boundary`.

This produces `seeds` with columns `x, y, val`.

---

## 3. Single-round BLUE trace-back (radius = r)

### What this round produces
- `seeds`: chosen boundary seeds on circle of radius `r`
- `paths_back`: one BLUE path per seed
- `coverage`: optional coverage map of all BLUE paths

### Code pattern
1. Build boundary candidates on radius `r`.
2. Pick `K_boundary` seeds with minimum separation.
3. For each seed, run `rw_path_2d_multi_to_center(...)` with `radius = r`.
4. Optionally compute coverage via `build_coverage()`.

---

## 4. Two-round (or multi-round) Sholl-style trace-back

### Manual two-round example
- Round 1: `radius = 120`
- Round 2: `radius = 70`

Each round is run independently with its own boundary seeds and paths, then overlaid.

### Archived objects in code
- `blue_120`: list containing `radius_lim`, `(x0,y0)`, `H,W`, `seeds`, `paths_back`, `coverage1`
- `blue_90`: (named in code) but `radius_lim = 70` in the example block

> Note: object names may not match radius value (e.g., `blue_90` storing radius 70).  
> Recommended: rename to `blue_r70`, `blue_r120`, etc.

---

## 5. Automated multi-round trace-back (`auto_blue_3rounds`)

### What it does
Generates a list of BLUE trace-back rounds with shrinking radii:
- radii: `radius_max * shrink^(0:(n_round-1))` (rounded, lower bounded)

For each round k:
1. pick boundary seeds on radius `r_k`
2. trace BLUE paths back to soma inside `r_k`
3. archive round result as a list element

### Output
`blue_auto`: list of length `n_round`, each element has:
- `radius_lim`, `seeds`, `paths_back`, `coverage1`, plus metadata

---

## 6. Angle-based support filtering (consistency check)

### Helper: seed angle
`angle_of_seed(x,y,x0,y0) = atan2(y-y0, x-x0)`

### Helper: circular distance
`circ_diff(a,b) = min(|a-b|, 2π-|a-b|)`

### Filter logic (round1-by-support)
For each seed in round1, keep it only if its angle is “supported” by seeds in selected later rounds within tolerance `angle_tol_deg`:
- Support event in a round r: there exists at least one seed angle in round r within ±tol of the round1 seed angle.
- The implementation can require:
  - support from **all** specified rounds (strict), or
  - support from **any** specified round (lenient), depending on the chosen boolean rule.

### What it modifies
- Only the target round’s `seeds` and `paths_back` are subsetted.
- Coverage for that round can be recomputed.

---

## 7. Plotting / visualization blocks

### `plot_blue_only(prob2d, blue_obj, main)`
Plots:
- background `prob2d`
- circle boundary of `blue_obj$radius_lim`
- BLUE paths
- soma marker and seeds

### Overlay plots
Overlay multiple rounds:
- draw outer circle once
- draw all paths (sometimes in red)
- optionally show inner circles for debugging

---

## 8. Minimal “how to run” checklist

1. Obtain `prob2d` and set `(x0,y0)`.
2. Choose `radius_max` and number of rounds.
3. Run:
   - manual rounds (`blue_120`, `blue_r70`) **or**
   - `blue_auto <- auto_blue_3rounds(...)`
4. (Optional) apply angle support filter to round1:
   - `blue_auto <- filter_round1_by_support(...)`
5. Plot:
   - per-round summaries with `plot_blue_only()`
   - full overlay with `plot_overlay_all_rounds()` or custom overlay code

---

## 9. Parameter notes (practical)

- `alpha_rw`: too small → paths ignore signal; too large → paths get stuck on local maxima.
- `n_cand_rw`: higher → more stable sampling but slower.
- `step_rw`: larger steps can jump across gaps but may overshoot fine structures.
- `min_sep_boundary`: controls seed diversity around the circle; too small → redundant seeds.
- `shrink`: controls how quickly radii decrease across rounds; typical 0.5–0.7.

---

## 10. Known implementation caveats

- Ensure `prob2d` indexing uses `[y, x]` consistently (as in the code).
- Object naming should reflect actual radius values to avoid confusion.
- If `prob2d` has many near-zero values, the density of valid paths can drop; consider smoothing or adjusting `alpha`.

---
