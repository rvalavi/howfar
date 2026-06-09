# howfar

[![R-CMD-check](https://github.com/rvalavi/howfar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rvalavi/howfar/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

Fast Euclidean distance rasters from spatial features (points, lines, polygons), with a parallel Rust backend.

Given a `SpatVector` and a `SpatRaster` template, `distance_to()` returns a raster where each cell holds the distance to the nearest feature. The same function works for all three geometry types. Distances are returned in **metres by default** (for any CRS); pass `unit = "mapunit"` for raw CRS-unit distances.

```r
library(howfar)
library(terra)

template <- rast(extent = ext(0, 1000, 0, 1000), resolution = 10)
line <- vect("LINESTRING(0 500, 1000 500)")

d <- distance_to(line, template)
plot(d)
```

## Why

`terra::distance()` is the standard tool for vector-to-raster distance, but it can be slow for large rasters. `howfar` uses a different approach: rasterize the features once, then run a parallel exact linear-time Euclidean Distance Transform (Felzenszwalb & Huttenlocher 2012) in Rust with [rayon](https://crates.io/crates/rayon). For large rasters this is generally much faster.

## Important: accuracy caveat — use cautiously

**`howfar` is faster but less accurate than `terra::distance()`.**

Because features are discretized onto the raster grid before computing distance, the result is accurate only to *roughly one cell*. Cells inside or on a feature get distance `0`; cells beyond get the Euclidean distance to the nearest source **cell**, not to the exact line/polygon edge in CRS coordinates.

By contrast, `terra::distance()` computes the true geometric distance to the vector feature, accurate at sub-cell precision.

**Rule of thumb**:

- Use `howfar` when (a) your cell size is small relative to the precision you need, or (b) speed matters more than sub-cell accuracy.
- Use `terra::distance()` when you need exact geometric distance — for example, when comparing distances on the order of a fraction of a cell.

If you're unsure whether the accuracy is good enough for your use case, compare against `terra::distance()` on a small subset of your data first.

## Distance units

`distance_to()` returns distances in **metres by default** (`unit = "meter"`), for any coordinate reference system (CRS):

- **Projected CRS** — the planar Euclidean distance between cell centers, converted from the CRS linear unit to metres (e.g. feet → metres). Exact.
- **Geographic CRS** (longitude/latitude) — the **Haversine** great-circle distance in metres.

Pass `unit = "mapunit"` to instead get the raw distance transform in the template's own CRS units (cell distance × resolution) — faster, and a true metric distance only when the CRS is projected with square cells.

In both modes the nearest source cell is always found by the Euclidean distance transform on the raster grid (in cell/index space, "as is"); for `unit = "meter"` the distance to that nearest cell is then re-measured with the appropriate metric.

> [!NOTE]
> **Lon/lat distances are approximate.** For geographic CRSs the result is an approximation, for two reasons: (1) the nearest source is chosen in grid space (degrees), which can differ from the geodesically nearest source over very large extents or near the poles; and (2) Haversine assumes a spherical Earth. Over regional extents at moderate latitudes this is effectively exact; over continental/global extents or high latitudes it can drift from the true geodesic distance — reproject to a suitable projected CRS first for best accuracy.

## Installation

Requires a Rust toolchain — install via [rustup](https://rustup.rs/):

```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then in R:

```r
# install.packages("devtools")
devtools::install_github("rvalavi/howfar")
```

For local development:

```r
devtools::install()    # builds the Rust crate and installs the R package
devtools::test()       # runs the testthat suite
```

## Usage

```r
library(howfar)
library(terra)

template <- rast(extent = ext(-5.5, 5.5, -5.5, 5.5), resolution = 1)

# Points
pt  <- vect("POINT(0 0)")
d_pt <- distance_to(pt, template)

# Lines
line  <- vect("LINESTRING(0 -5, 0 5)")
d_ln  <- distance_to(line, template)

# Polygons (cells inside the polygon = 0)
poly  <- vect("POLYGON((-1 -1, 1 -1, 1 1, -1 1, -1 -1))")
d_pol <- distance_to(poly, template)

# Distance to polygon boundary only (positive inside and outside)
d_bound <- distance_to(as.lines(poly), template)

# Combine multiple geometry types
d_any <- min(d_pt, d_ln, d_pol)

# Cap parallelism (default uses all logical cores)
d <- distance_to(line, template, n_cores = 4)
```

## API

```r
distance_to(x, template, touches = TRUE, n_cores = NULL, unit = c("meter", "mapunit"))
```

| argument   | type                     | meaning                                                                 |
|------------|--------------------------|-------------------------------------------------------------------------|
| `x`        | `SpatVector` / `sf` / `sfc` | points, lines, *or* polygons (one type per call)                       |
| `template` | `SpatRaster`             | defines output extent, resolution, and CRS                              |
| `touches`  | logical                  | passed to `terra::rasterize()`; `TRUE` marks every cell touched         |
| `n_cores`  | integer or `NULL`        | number of threads; `NULL` = all logical cores; `1` = serial             |
| `unit`     | `"meter"` / `"mapunit"`  | output units; `"meter"` (default) = metres, `"mapunit"` = CRS units     |

Returns a single-layer `SpatRaster` with distances in metres (`unit = "meter"`) or CRS units (`unit = "mapunit"`). See [Distance units](#distance-units).

## Geometry-specific behavior

- **Points** — distance to the nearest cell containing a point.
- **Lines** — distance to the nearest cell crossed by a line.
- **Polygons** — cells inside the polygon get `0`; cells outside get the distance to the nearest cell that overlaps the polygon.
- **Polygon boundary only** (positive both inside and outside) — pass `terra::as.lines(poly)` instead of `poly`.

For mixed geometry types, call `distance_to()` separately on each `SpatVector` and combine with `terra::min()` — a `SpatVector` can hold only one geometry type at a time.

## How it works

1. Rasterize the input features onto `template` with `terra::rasterize(field = 1, background = 0)` — produces a binary mask.
2. Pass the mask to Rust as a row-major numeric vector.
3. Run an exact 2D Euclidean Distance Transform using Felzenszwalb & Huttenlocher's (2012) linear-time lower-envelope-of-parabolas algorithm:
   - 1D EDT along each row, rows processed in parallel via rayon.
   - Transpose.
   - 1D EDT along each row of the transposed matrix (i.e., each column of the original).
   - Transpose back.
4. Take the square root to get distances in cell units, then convert to the requested `unit`: for `"mapunit"`, multiply by cell size (CRS units); for `"meter"`, also track each cell's nearest source cell (the *feature transform*) and measure the metric distance between cell centers — planar Euclidean (scaled to metres) for projected CRSs, or Haversine for lon/lat.

The 1D EDT is `O(n)` per row/column and the work is embarrassingly parallel across rows and across columns.

## Limitations

- `unit = "mapunit"` assumes square cells (warns if not); `unit = "meter"` handles non-square cells correctly.
- Lon/lat distances (`unit = "meter"`) are approximate — see [Distance units](#distance-units).
- One geometry type per call.
- In-memory only — the raster must fit in RAM.
- Accuracy is bounded by the cell size (see caveat above).

## License

MIT.

## References

Felzenszwalb, P. F. & Huttenlocher, D. P. (2012). *Distance Transforms of Sampled Functions*. Theory of Computing, 8(1), 415–428.
