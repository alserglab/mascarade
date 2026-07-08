# mascarade label-placement integration — implementation plan

Canonical plan for packaging the boundary-seed label placer. **Reread this file and the
integration design document (`label_placement_integration_design/report.html`) at the start
of every stage** to prevent drift. Each stage ends by committing and re-reading both.

## Reference
- Design doc: `/workspace/label_placement_integration_design/report.html`
- Prototype (source of truth for the algorithm): `/workspace/label_iters/polish_cpp.cpp`,
  `/workspace/label_iters/multi_setup.R` + `common_setup.R`, `tmp/test_onepass.R`.
- Package: `R/`, (new) `src/`, `tests/testthat/`, `DESCRIPTION`.

## Context that shapes the plan (from recon)
1. **Package is currently pure R (no `src/`).** Integration makes mascarade a **compiled
   package** (Rcpp + LinkingTo BH). BH is installed; Rcpp must be added.
2. **Poles are not computed in the package.** `generateMask()` returns only mask-border
   points. Poles (leader targets) are new — source = `polylabelr::poi` (new dependency).
3. **Labels today are positioned by `geom_mark_shape()` / ggforce mark-label machinery**,
   not by us. Feeding our placement in is a real integration point (Stage 1 spike pins it),
   not the one-liner the design doc implies.

## Key decisions (confirmed 2026-07-08)
- [x] mascarade becomes a compiled package (Rcpp + BH). **OK.**
- [x] Add `polylabelr` (poles) + `Rcpp` to Imports, `BH`,`Rcpp` to LinkingTo. **OK.**
- [x] Branch `feat/label-placement`; keep prototype until Stage 6. **OK.**
- **`geom_mark_shape()` is a local copy of ggforce** — its labeling section can be replaced
  directly. This de-risks Stage 1(b): no upstream interception needed, just edit the copy.

---

## Stage 0 — Groundwork & regression baseline
- Create branch `feat/label-placement`.
- Capture **golden `layoutScore`** for the reference datasets (COPD orig+r1/r2/r3 + gallery)
  from the current prototype pipeline → `tests/testthat/fixtures/golden_scores.rds`.
  This is the anchor every later stage must not regress against.
- **Commit** ("plan + golden baselines"). **Reread plan + design doc.**

## Stage 1 — De-risk spikes (no package changes yet)
- (a) **BH R-tree spike**: minimal Rcpp+BH program building a Boost.Geometry rtree over a few
  polygons and answering a box-intersect query. If it doesn't compile/run here → fall back to
  the integral-image grid (documented rollback) and note it.
- (b) **Label hook**: trace exactly where `geom_mark_shape()`/`ggforce_mark_label` set label
  positions, and identify where our `(cx,cy)` + leader must be injected.
- (c) **Poles**: confirm `polylabelr::poi` on the mask polygons reproduces the prototype poles.
- Write findings into this file.
- **Commit** ("spike notes"). **Reread plan + design doc.**

## Stage 2 — C++ kernels → `src/`
- DESCRIPTION: add `Rcpp` (Imports), `Rcpp`,`BH` (LinkingTo); `src/Makevars` if needed.
- Port the six kernels, reorganised + renamed:
  `src/geometry.{h,cpp}` (BH rtree box-fit + `bg::` predicates),
  `src/candidates.cpp` (radialCandidates), `src/assign.cpp` (hungarian),
  `src/refine.cpp` (oneMoveSweep, twoMoveBnB = length B&B + **lexicographic pick**),
  `src/polish.cpp` (forcePolish), `src/foreign.cpp` (foreignLength).
- Drop dead kernels (mincon/uncross/greedy_opt/…/core_bt).
- **Verify**: `devtools::load_all()` compiles; kernels reproduce prototype outputs on a fixture.
- **Commit** ("C++ kernels + BH box-fit"). **Reread plan + design doc.**

## Stage 3 — R driver `R/place.R`
- `placeLabels(geom, xlim, ylim, aspect, gp, MU=55, iters=120)` + stage wrappers
  (`genCandidates, seedLayout, refineOneMove, refineTwoMove, polishLayout`), the `prep` glue
  (candidate pool + per-label boundary slot, `rows`/index maps, `GEO` columns), `labelBoxes`
  (real text metrics via systemfonts), and the geom builder (poles + rtree).
- **Verify**: `placeLabels` on reference datasets → `layoutScore` matches golden (0/0/0, no worse).
- **Commit** ("R placement driver"). **Reread plan + design doc.**

## Stage 4 — `fancyMask` draw-stage integration
- `generateMask()` stays a plain data.table; add poles as a companion (attribute or
  second return) so `fancyMask` can build `geom` once.
- `fancyMask`: build `geom` (poles + rtree) once; place labels via a `makeContent()` dynamic
  grob calling `placeLabels`, wired into the label hook found in Stage 1; real metrics from the
  live device.
- **Verify**: a real ggplot renders with placed leaders; re-render at a different size stays
  conflict-free (re-place on view change).
- **Commit** ("fancyMask draw-stage placement"). **Reread plan + design doc.**

## Stage 5 — Score helper + tests
- `tests/testthat/helper-score.R` (`layoutScore`, `scoreBetter`).
- Tests: golden regression (`scoreBetter(new, golden) || equal`), re-place-on-view-change
  conflict-freeness, feasibility fallback (boundary slot / overflow flag present).
- **Verify**: `devtools::test()` green.
- **Commit** ("score helper + tests"). **Reread plan + design doc.**

## Stage 6 — Cleanup, docs, check
- Remove prototype (`label_iters/`, `tmp/`, old reports) per design doc "drop" list.
- `devtools::document()`; NEWS; finalise DESCRIPTION.
- `devtools::check()` clean (or documented residual notes).
- **Commit** ("cleanup + docs + check"). **Reread plan + design doc.**

## Deferred (explicitly out of scope for this pass)
- Overflow / "not enough space" behaviour beyond the boundary-slot fallback + `overflow` flag.
- Placement memoisation and resize warm-start (v2).
- Auto-expand / fit policy in `fancyMask`.
