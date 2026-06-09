test_that("rust_edt: zero at source, correct distances elsewhere", {
    # 5x5 grid, single source at row 3 col 3 (1-indexed).
    # Row-major flat index for (r, c) 1-indexed: (r-1) * 5 + c
    mask <- rep(0, 25)
    mask[(3 - 1) * 5 + 3] <- 1  # index 13

    d <- rust_edt(mask, 5L, 5L, 0L)

    # Source cell distance = 0
    expect_equal(d[13], 0)

    # Cell (3, 1): 2 cells left of source
    expect_equal(d[(3 - 1) * 5 + 1], 2)

    # Cell (1, 1): top-left corner, distance sqrt(8)
    expect_equal(d[1], sqrt(8), tolerance = 1e-10)

    # Cell (5, 5): bottom-right corner, distance sqrt(8)
    expect_equal(d[25], sqrt(8), tolerance = 1e-10)
})

test_that("rust_edt: vertical line gives column-distance raster", {
    nrows <- 7L
    ncols <- 11L
    mask <- rep(0, nrows * ncols)
    # Mark column 6 (1-indexed) of every row
    for (r in seq_len(nrows)) {
        mask[(r - 1) * ncols + 6] <- 1
    }

    d <- rust_edt(mask, nrows, ncols, 0L)
    m <- matrix(d, nrow = nrows, ncol = ncols, byrow = TRUE)

    # Every row should have the same column-distance pattern: |c - 6|
    expected_row <- abs(seq_len(ncols) - 6)
    for (r in seq_len(nrows)) {
        expect_equal(as.numeric(m[r, ]), as.numeric(expected_row),
                     tolerance = 1e-10)
    }
})

test_that("rust_edt: n_cores=1 gives same answer as default (parallel)", {
    mask <- rep(0, 25)
    mask[13] <- 1
    expect_equal(rust_edt(mask, 5L, 5L, 1L),
                 rust_edt(mask, 5L, 5L, 0L))
})

test_that("distance_to: vertical line on a small raster", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-5.5, 5.5, -5.5, 5.5),
        resolution = 1,
        crs = ""
    )
    line <- terra::vect("LINESTRING(0 -5, 0 5)", crs = "")

    d <- distance_to(line, template)

    expect_s4_class(d, "SpatRaster")
    expect_equal(terra::nrow(d), 11)
    expect_equal(terra::ncol(d), 11)

    m <- terra::as.matrix(d, wide = TRUE)
    middle_row <- as.numeric(m[6, ])
    expect_equal(min(middle_row), 0)
    expect_equal(max(middle_row), 5, tolerance = 1)
})

test_that("distance_to: single point", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-5.5, 5.5, -5.5, 5.5),
        resolution = 1,
        crs = ""
    )
    pt <- terra::vect("POINT(0 0)", crs = "")

    d <- distance_to(pt, template)

    expect_s4_class(d, "SpatRaster")
    m <- terra::as.matrix(d, wide = TRUE)

    # Cell containing (0, 0) is row 6, col 6 -> distance 0
    expect_equal(m[6, 6], 0)
    # Top-left corner (-5, 5) -> distance sqrt(50) ~ 7.07
    expect_equal(m[1, 1], sqrt(50), tolerance = 1e-6)
    # Bottom-right corner (5, -5) -> distance sqrt(50)
    expect_equal(m[11, 11], sqrt(50), tolerance = 1e-6)
})

test_that("distance_to: polygon (inside = 0, outside = positive)", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-5.5, 5.5, -5.5, 5.5),
        resolution = 1,
        crs = ""
    )
    # Square polygon from (-1, -1) to (1, 1)
    poly <- terra::vect("POLYGON((-1 -1, 1 -1, 1 1, -1 1, -1 -1))",
                        crs = "")

    d <- distance_to(poly, template)

    expect_s4_class(d, "SpatRaster")
    m <- terra::as.matrix(d, wide = TRUE)

    # Cell at polygon center (0, 0) -> inside, distance 0
    expect_equal(m[6, 6], 0)
    # Cell at polygon corner (-1, 1), row 5 col 5 -> on boundary, distance 0
    expect_equal(m[5, 5], 0)
    # Far corner (-5, 5) -> nearest marked cell is (-1, 1), distance sqrt(32)
    expect_equal(m[1, 1], sqrt(32), tolerance = 1e-6)
})

test_that("distance_to: polygon boundary via terra::as.lines()", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-5.5, 5.5, -5.5, 5.5),
        resolution = 1,
        crs = ""
    )
    poly <- terra::vect("POLYGON((-1 -1, 1 -1, 1 1, -1 1, -1 -1))",
                        crs = "")

    # Distance to polygon boundary (interior cells get positive distance)
    d_boundary <- distance_to(terra::as.lines(poly), template)

    m <- terra::as.matrix(d_boundary, wide = TRUE)
    # Cell at polygon center (0, 0) -> NOT on boundary -> distance > 0
    expect_gt(m[6, 6], 0)
    # Cell on polygon corner -> on boundary, distance 0
    expect_equal(m[5, 5], 0)
})

test_that("distance_to: n_cores argument is accepted", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-5.5, 5.5, -5.5, 5.5),
        resolution = 1,
        crs = ""
    )
    line <- terra::vect("LINESTRING(0 -5, 0 5)", crs = "")

    d1 <- distance_to(line, template, n_cores = 1)
    d2 <- distance_to(line, template, n_cores = NULL)
    expect_equal(terra::values(d1), terra::values(d2))

    expect_error(distance_to(line, template, n_cores = -1),
                 "non-negative")
})

test_that("distance_to: planar on a projected CRS returns metres (and km)", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(0, 1000, 0, 1000),
        resolution = 10,
        crs = "EPSG:3857"  # Web Mercator: projected, linear unit = metre
    )
    line <- terra::vect("LINESTRING(505 0, 505 1000)", crs = "EPSG:3857")

    d_m  <- distance_to(line, template, method = "planar", unit = "m")
    d_km <- distance_to(line, template, method = "planar", unit = "km")

    vm <- terra::values(d_m)[, 1]
    vk <- terra::values(d_km)[, 1]

    expect_equal(min(vm), 0)                        # source cells
    expect_true(all(is.finite(vm)))
    expect_equal(vk, vm / 1000, tolerance = 1e-12)  # km = m / 1000
})

test_that("distance_to: auto picks haversine for lon/lat (metres)", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-2.5, 2.5, -2.5, 2.5),
        resolution = 1,
        crs = "EPSG:4326"
    )
    pt <- terra::vect("POINT(0 0)", crs = "EPSG:4326")

    d  <- distance_to(pt, template)                       # auto -> haversine
    dh <- distance_to(pt, template, method = "haversine")
    expect_equal(terra::values(d), terra::values(dh))

    m <- terra::as.matrix(d, wide = TRUE)
    one_degree <- 6371000 * pi / 180  # ~111195 m
    expect_equal(m[3, 3], 0)                              # the point's own cell
    expect_equal(m[3, 4], one_degree, tolerance = 1e-4)  # 1 deg lon at equator
    expect_equal(m[2, 3], one_degree, tolerance = 1e-4)  # 1 deg lat

    d_km <- distance_to(pt, template, unit = "km")
    expect_equal(terra::values(d_km), terra::values(d) / 1000, tolerance = 1e-9)
})

test_that("distance_to: planar on lon/lat gives map units (degrees)", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-2.5, 2.5, -2.5, 2.5),
        resolution = 1,
        crs = "EPSG:4326"
    )
    pt <- terra::vect("POINT(0 0)", crs = "EPSG:4326")

    d <- distance_to(pt, template, method = "planar")
    m <- terra::as.matrix(d, wide = TRUE)

    # Distances are in degrees: 1 cell = 1 degree here.
    expect_equal(m[3, 3], 0)
    expect_equal(m[3, 4], 1, tolerance = 1e-9)        # 1 degree east
    expect_equal(m[1, 1], sqrt(8), tolerance = 1e-9)  # corner: (2, 2) cells away
})

test_that("distance_to: haversine on a projected CRS errors", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(0, 1000, 0, 1000),
        resolution = 10,
        crs = "EPSG:3857"
    )
    line <- terra::vect("LINESTRING(505 0, 505 1000)", crs = "EPSG:3857")

    expect_error(distance_to(line, template, method = "haversine"),
                 "geographic")
})

test_that("distance_to: haversine without a CRS warns (assumes lon/lat)", {
    skip_if_not_installed("terra")

    template <- terra::rast(
        extent = terra::ext(-2.5, 2.5, -2.5, 2.5),
        resolution = 1,
        crs = ""
    )
    pt <- terra::vect("POINT(0 0)", crs = "")

    expect_warning(distance_to(pt, template, method = "haversine"), "lon/lat")
})

test_that("distance_to: haversine finds the great-circle-nearest source", {
    skip_if_not_installed("terra")

    # cols centered at lon 0..10, rows centered at lat 81..75.
    template <- terra::rast(
        extent = terra::ext(-0.5, 10.5, 74.5, 81.5),
        resolution = 1,
        crs = "EPSG:4326"
    )
    # Two point sources. For the cell at (lon 0, lat 80):
    #   A = (5, 80): 5 deg lon  -> ~96.5 km on the ground (the TRUE nearest)
    #   B = (0, 78): 2 deg lat  -> ~222 km (nearer only in raw degrees)
    src <- rbind(
        terra::vect("POINT(5 80)", crs = "EPSG:4326"),
        terra::vect("POINT(0 78)", crs = "EPSG:4326")
    )
    d <- distance_to(src, template, method = "haversine")
    m <- terra::as.matrix(d, wide = TRUE)

    hav <- function(lon1, lat1, lon2, lat2) {
        R <- 6371000; d2r <- pi / 180
        dlat <- (lat2 - lat1) * d2r; dlon <- (lon2 - lon1) * d2r
        a <- sin(dlat / 2)^2 + cos(lat1 * d2r) * cos(lat2 * d2r) * sin(dlon / 2)^2
        2 * R * asin(pmin(1, sqrt(a)))
    }

    val <- m[2, 1]  # cell at (lon 0, lat 80)
    expect_equal(val, hav(0, 80, 5, 80), tolerance = 1e-6)  # distance to A
    expect_lt(val, hav(0, 80, 0, 78))                       # NOT the degree-nearest B
})
