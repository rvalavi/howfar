#' Euclidean distance raster from spatial features
#'
#' Returns a raster where each cell value is the shortest Euclidean distance
#' (in CRS units) to the nearest feature in `x`. Points, lines, and polygons
#' are all supported. Features are rasterized onto `template`, then a parallel
#' exact Euclidean distance transform is computed in Rust.
#'
#' @param x A `SpatVector` of points, lines, *or* polygons (one geometry type
#'   per call — a `SpatVector` cannot hold multiple types). `sf`/`sfc` objects
#'   are also accepted and converted with [terra::vect()].
#' @param template A `SpatRaster` defining the extent, resolution, and CRS of
#'   the output raster. Values in `template` are ignored; only its geometry
#'   is used.
#' @param touches Logical. Passed to [terra::rasterize()]: if `TRUE` (default),
#'   every cell touched by a feature is marked as a source. If `FALSE`, only
#'   cells whose center falls on a feature are marked. For polygons, `FALSE`
#'   marks cells whose center lies inside the polygon. Has no effect for
#'   points.
#' @param n_cores Integer or `NULL`. Number of threads to use for the parallel
#'   distance transform. `NULL` (default) uses all logical cores available
#'   (rayon's default). A positive integer caps parallelism at that value;
#'   `1` runs serially.
#'
#' @return A single-layer `SpatRaster` with the distance to the nearest feature
#'   in each cell.
#'
#' @details
#' Distance assumes square cells. If the template has non-square cells, a
#' warning is issued and the x-resolution is used as the scale factor.
#'
#' **Geometry-specific semantics:**
#' * **Points**: distance is to the nearest cell containing a point.
#' * **Lines**: distance is to the nearest cell crossed by a line (with
#'   `touches = TRUE`).
#' * **Polygons**: cells inside the polygon get distance `0`; cells outside
#'   get the distance to the nearest cell that overlaps the polygon. For
#'   distance to the polygon *boundary* only (positive both inside and
#'   outside), convert polygons to lines first with [terra::as.lines()] and
#'   pass that.
#'
#' To compute distance to a mix of geometry types, call `distance_to()`
#' separately for each `SpatVector` and combine with `pmin()` or
#' [terra::min()].
#'
#' If no feature intersects any cell, every cell in the result will be `Inf`.
#'
#' @examples
#' \dontrun{
#' library(terra)
#' template <- rast(extent = ext(-5.5, 5.5, -5.5, 5.5), resolution = 1)
#'
#' # Distance to a line
#' line <- vect("LINESTRING(0 -5, 0 5)")
#' d_line <- distance_to(line, template)
#'
#' # Distance to a point
#' point <- vect("POINT(0 0)")
#' d_point <- distance_to(point, template)
#'
#' # Distance to a polygon (cells inside the polygon = 0)
#' poly <- vect("POLYGON((-1 -1, 1 -1, 1 1, -1 1, -1 -1))")
#' d_poly <- distance_to(poly, template)
#'
#' # Distance to polygon boundary only (positive both inside and outside)
#' d_bound <- distance_to(as.lines(poly), template)
#'
#' # Distance to combined point + line features
#' d_combined <- min(d_line, d_point)
#'
#' # Limit to 4 threads
#' d <- distance_to(line, template, n_cores = 4)
#' }
#'
#' @export
distance_to <- function(x, template, touches = TRUE, n_cores = NULL) {
    if (!inherits(template, "SpatRaster")) {
        stop("`template` must be a SpatRaster.")
    }

    if (inherits(x, "sf") || inherits(x, "sfc")) {
        x <- terra::vect(x)
    }
    if (!inherits(x, "SpatVector")) {
        stop("`x` must be a SpatVector, sf, or sfc object.")
    }

    if (is.null(n_cores)) {
        n_cores <- 0L  # 0 signals "use rayon default" on the Rust side
    } else {
        n_cores <- as.integer(n_cores)
        if (length(n_cores) != 1L || is.na(n_cores) || n_cores < 0L) {
            stop("`n_cores` must be NULL or a non-negative integer.")
        }
    }

    mask_rast <- terra::rasterize(
        x, template,
        field = 1, background = 0, touches = touches
    )

    mask_vec <- as.numeric(terra::values(mask_rast))
    mask_vec[is.na(mask_vec)] <- 0

    nrows <- as.integer(terra::nrow(mask_rast))
    ncols <- as.integer(terra::ncol(mask_rast))

    edt_cells <- rust_edt(mask_vec, nrows, ncols, n_cores)

    cell_x <- terra::xres(template)
    cell_y <- terra::yres(template)
    if (abs(cell_x - cell_y) > 1e-10 * max(cell_x, cell_y)) {
        warning(
            "Non-square cells (xres=", cell_x, ", yres=", cell_y,
            "); distance assumes square cells (using xres)."
        )
    }

    result <- terra::rast(template)
    terra::values(result) <- edt_cells * cell_x
    names(result) <- "distance"
    result
}
