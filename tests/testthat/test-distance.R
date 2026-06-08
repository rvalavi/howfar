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
