#' Distance raster from spatial features
#'
#' Returns a raster where each cell value is the shortest distance to the
#' nearest feature in `x`. Points, lines, and polygons are all supported.
#' Features are rasterized onto `template`, then a parallel exact Euclidean
#' distance transform is computed in Rust to locate, for every cell, its
#' nearest source cell. The distance to that cell is then measured with the
#' chosen `method`.
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
#' @param method Character, one of `"auto"` (default), `"planar"`, or
#'   `"haversine"`. How the distance between a cell and its nearest source cell
#'   is measured. `"auto"` uses `"haversine"` for geographic (longitude/
#'   latitude) CRSs and `"planar"` otherwise. See Details.
#' @param unit Character, one of `"m"` (metres, default) or `"km"`. Output unit
#'   for *metric* results. Ignored when the result is in map units — i.e. a
#'   `"planar"` distance on a geographic or CRS-less template (see Details).
#'
#' @return A single-layer `SpatRaster` of distances to the nearest feature, in
#'   metres/kilometres for metric results, or in the template's map units for a
#'   planar distance on a geographic/CRS-less template.
#'
#' @details
#' **Methods.**
#' * **`"planar"`** — Euclidean distance between cell centers in the template's
#'   own coordinates. For a projected CRS this is a true metric distance,
#'   returned in metres (converted from the CRS linear unit via
#'   [terra::linearUnits()], so feet etc. become metres); `unit` selects m/km.
#'   For a geographic CRS the coordinates are degrees, so the result is in
#'   **map units (degrees)** and `unit` is ignored. Non-square cells are handled
#'   correctly (the actual center-to-center distance is measured).
#' * **`"haversine"`** — great-circle distance in metres on a sphere of radius
#'   6,371,000 m. The nearest source is found by a true spherical
#'   nearest-neighbour search (correct at any latitude and across the
#'   antimeridian), and the Haversine distance to it is returned. Requires a
#'   geographic CRS; on a projected CRS it errors (reproject to lon/lat first).
#'   `unit` selects m/km.
#' * **`"auto"`** (default) — `"haversine"` for geographic CRSs, `"planar"`
#'   otherwise. Distances therefore default to metres for both projected and
#'   geographic inputs, and to map units only when the CRS is missing.
#'
#' **Nearest source.** For `"planar"` the nearest source cell is found by the
#' exact Euclidean distance transform on the raster grid. For `"haversine"` it
#' is found by a great-circle nearest-neighbour search, so the reported distance
#' is to the genuinely closest source on the sphere (not merely the grid-nearest
#' one).
#'
#' **Accuracy (`"haversine"`).** Because the nearest source is found correctly
#' by great-circle distance, the only approximations are (1) features are
#' rasterized to cell centers (accurate to ~1 cell, as everywhere in `howfar`),
#' and (2) Haversine assumes a spherical Earth (vs. an ellipsoidal geodesic,
#' typically <0.5%). For sub-cell or ellipsoidal precision, use
#' [terra::distance()].
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
#'
#' # Projected template (metres) -> planar distance in metres
#' template <- rast(extent = ext(0, 1000, 0, 1000), resolution = 10,
#'                  crs = "EPSG:3857")
#' line <- vect("LINESTRING(500 0, 500 1000)", crs = "EPSG:3857")
#' d <- distance_to(line, template)              # metres (auto -> planar)
#' d_km <- distance_to(line, template, unit = "km")
#'
#' # Geographic template -> great-circle distance in metres
#' tmpl_ll <- rast(extent = ext(-2.5, 2.5, -2.5, 2.5), resolution = 0.1,
#'                 crs = "EPSG:4326")
#' pt <- vect("POINT(0 0)", crs = "EPSG:4326")
#' d_ll <- distance_to(pt, tmpl_ll)              # metres (auto -> haversine)
#'
#' # Force planar on a lon/lat template -> distance in degrees (map units)
#' d_deg <- distance_to(pt, tmpl_ll, method = "planar")
#'
#' # Limit to 4 threads
#' d <- distance_to(line, template, n_cores = 4)
#' }
#'
#' @export
distance_to <- function(x, template, touches = TRUE, n_cores = NULL,
                        method = c("auto", "planar", "haversine"),
                        unit = c("m", "km")) {
    unit_supplied <- !missing(unit)
    method <- match.arg(method)
    unit <- match.arg(unit)

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

    # is.lonlat() warns on an unknown CRS; we handle NA explicitly below.
    lonlat <- suppressWarnings(terra::is.lonlat(template))

    if (method == "auto") {
        method <- if (isTRUE(lonlat)) "haversine" else "planar"
    }

    xmin <- terra::xmin(template)
    ymax <- terra::ymax(template)
    xres <- terra::xres(template)
    yres <- terra::yres(template)

    if (method == "haversine") {
        if (is.na(lonlat)) {
            warning(
                "`template` has no CRS; assuming lon/lat coordinates for the ",
                "Haversine distance."
            )
        } else if (!isTRUE(lonlat)) {
            stop(
                "method = \"haversine\" requires a geographic (lon/lat) CRS. ",
                "Reproject `template` to lon/lat, or use method = \"planar\"."
            )
        }
        dist <- rust_geo_distance(
            mask_vec, nrows, ncols, n_cores,
            xmin, ymax, xres, yres
        )
        if (unit == "km") dist <- dist / 1000
    } else {
        # planar: Euclidean distance between cell centers, in CRS coordinates.
        metric <- FALSE
        m_per_unit <- 1
        if (isTRUE(lonlat)) {
            # lon/lat coordinates are degrees -> map units; unit not applicable.
            if (unit_supplied && unit == "km") {
                message(
                    "`unit` is ignored for a planar distance on a lon/lat CRS; ",
                    "the result is in map units (degrees)."
                )
            }
        } else if (is.na(lonlat)) {
            # Unknown CRS -> planar distance is in map units.
        } else {
            lu <- terra::linearUnits(template)
            if (length(lu) == 1L && !is.na(lu) && lu > 0) {
                m_per_unit <- lu  # convert CRS linear unit to metres
                metric <- TRUE
            } else {
                warning(
                    "Could not determine the CRS linear unit; planar distances ",
                    "are in map units."
                )
            }
        }
        dist <- rust_edt_meters(
            mask_vec, nrows, ncols, n_cores,
            xmin, ymax, xres, yres, FALSE, m_per_unit
        )
        if (metric && unit == "km") dist <- dist / 1000
    }

    result <- terra::rast(template)
    terra::values(result) <- dist
    names(result) <- "distance"
    result
}
