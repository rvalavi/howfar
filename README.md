# howfar

Fast Euclidean distance rasters from spatial features (points, lines, polygons), with a parallel Rust backend.

Given a `SpatVector` and a `SpatRaster` template, `distance_to()` returns a raster where each cell holds the Euclidean distance to the nearest feature. The same function works for all three geometry types.

``` r
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

-   Use `howfar` when (a) your cell size is small relative to the precision you need, or (b) speed matters more than sub-cell accuracy.
-   Use `terra::distance()` when you need exact geometric distance — for example, when comparing distances on the order of a fraction of a cell.

If you're unsure whether the accuracy is good enough for your use case, compare against `terra::distance()` on a small subset of your data first.

## Installation

Requires a Rust toolchain — install via [rustup](https://rustup.rs/):

``` sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then in R:

``` r
# install.packages("devtools")
devtools::install_github("rvalavi/howfar")
```

For local development:

``` r
devtools::install()    # builds the Rust crate and installs the R package
devtools::test()       # runs the testthat suite
```

## Usage

``` r
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

## Geometry-specific behavior

-   **Points** — distance to the nearest cell containing a point.
-   **Lines** — distance to the nearest cell crossed by a line.
-   **Polygons** — cells inside the polygon get `0`; cells outside get the distance to the nearest cell that overlaps the polygon.
-   **Polygon boundary only** (positive both inside and outside) — pass `terra::as.lines(poly)` instead of `poly`.

For mixed geometry types, call `distance_to()` separately on each `SpatVector` and combine with `terra::min()` — a `SpatVector` can hold only one geometry type at a time.

## How it works

1.  Rasterize the input features onto `template` with `terra::rasterize(field = 1, background = 0)` — produces a binary mask.
2.  Pass the mask to Rust as a row-major numeric vector.
3.  Run an exact 2D Euclidean Distance Transform using Felzenszwalb & Huttenlocher's (2012) linear-time lower-envelope-of-parabolas algorithm:
    -   1D EDT along each row, rows processed in parallel via rayon.
    -   Transpose.
    -   1D EDT along each row of the transposed matrix (i.e., each column of the original).
    -   Transpose back.
4.  Take the square root and multiply by cell size to convert from squared cell units to CRS-unit distances.

The 1D EDT is `O(n)` per row/column and the work is embarrassingly parallel across rows and across columns.

## Limitations

-   Assumes square cells (warns if not).
-   One geometry type per call.
-   In-memory only — the raster must fit in RAM.
-   Accuracy is bounded by the cell size (see caveat above).

## License

MIT.

## References

Felzenszwalb, P. F. & Huttenlocher, D. P. (2012). *Distance Transforms of Sampled Functions*. Theory of Computing, 8(1), 415–428.
