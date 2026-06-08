//! Parallel exact Euclidean distance transform for binary rasters.
//!
//! Implements Felzenszwalb & Huttenlocher (2012) linear-time EDT
//! ("Distance Transforms of Sampled Functions"). The algorithm is applied
//! once along rows, then once along columns, with rows/columns processed
//! in parallel via rayon. Result: O(n_cells) sequential work, embarrassingly
//! parallel across rows/columns.

use extendr_api::prelude::*;
use rayon::prelude::*;

/// 1D exact squared EDT of a sampled function.
///
/// Given an input array `f` of length `n`, computes
/// `d[q] = min over p in 0..n of ((q - p)^2 + f[p])` for each `q`.
///
/// Values of `f` greater than or equal to `inf` are treated as "no source"
/// (infinite height parabola — never wins). This avoids NaN from infinity
/// arithmetic in the lower-envelope construction.
fn dt_1d(f: &[f64], d: &mut [f64], inf: f64) {
    let n = f.len();
    debug_assert_eq!(d.len(), n);
    if n == 0 {
        return;
    }

    // Find the first non-background index.
    let mut first = 0usize;
    while first < n && f[first] >= inf {
        first += 1;
    }

    if first == n {
        // No source in this slice — all-infinity output.
        for x in d.iter_mut() {
            *x = inf;
        }
        return;
    }

    // Lower envelope: parabolas indexed by their apex location.
    //   v[k] is the apex of the k-th parabola in the envelope.
    //   z[k] is the x-coordinate where parabolas k-1 and k cross.
    //   z[0] = -inf, z[k_max+1] = +inf bracket the envelope.
    let mut v: Vec<usize> = vec![0; n];
    let mut z: Vec<f64> = vec![0.0; n + 1];

    v[0] = first;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    let mut k: usize = 0;

    for q in (first + 1)..n {
        let fq = f[q];
        if fq >= inf {
            // Infinite-height parabola never enters the lower envelope.
            continue;
        }
        let qf = q as f64;
        let q_sq = qf * qf;

        // Pop parabolas from the back of the envelope while the new one
        // dominates them. Loop terminates: at k=0, z[0]=-inf so s>z[0] for
        // any finite s (and we've filtered the inf case above).
        loop {
            let vk = v[k];
            let fvk = f[vk];
            let vkf = vk as f64;
            let s = ((fq + q_sq) - (fvk + vkf * vkf)) / (2.0 * (qf - vkf));

            if s <= z[k] {
                k -= 1;
            } else {
                k += 1;
                v[k] = q;
                z[k] = s;
                z[k + 1] = f64::INFINITY;
                break;
            }
        }
    }

    // Backward pass: for each output index q, walk z[] to find which
    // parabola wins at q, then evaluate it.
    let mut kr: usize = 0;
    for q in 0..n {
        let qf = q as f64;
        while z[kr + 1] < qf {
            kr += 1;
        }
        let vk = v[kr];
        let dq = qf - vk as f64;
        d[q] = dq * dq + f[vk];
    }
}

/// Transpose a row-major matrix of shape (in_rows × in_cols) into a
/// row-major matrix of shape (in_cols × in_rows). Parallel by output row.
fn transpose(input: &[f64], in_rows: usize, in_cols: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; input.len()];
    out.par_chunks_mut(in_rows)
        .enumerate()
        .for_each(|(c, dst)| {
            // dst is row c of the transposed matrix.
            // dst[r] = input[r * in_cols + c]
            for (r, slot) in dst.iter_mut().enumerate() {
                *slot = input[r * in_cols + c];
            }
        });
    out
}

/// Compute the squared exact Euclidean distance transform of a binary mask.
///
/// `mask` is row-major (cell-order from `terra::values()`) of size
/// `nrows × ncols`. A cell is a source if `mask[i] > 0.5`. The result is
/// the squared Euclidean distance to the nearest source, in *cell units*,
/// in the same row-major layout.
fn edt_squared(mask: &[f64], nrows: usize, ncols: usize) -> Vec<f64> {
    assert_eq!(mask.len(), nrows * ncols, "mask size mismatch");

    // Pick INF large enough that any "infinity" sentinel beats every real
    // squared distance, but small enough that arithmetic in the intersection
    // formula stays in a precise f64 range.
    //
    // Max squared 2D distance is nrows^2 + ncols^2; we use 4x that as a
    // safe sentinel.
    let max_dim = nrows.max(ncols) as f64;
    let inf = 4.0 * max_dim * max_dim + 4.0;

    // Initialize: sources -> 0, background -> inf.
    let init: Vec<f64> = mask
        .par_iter()
        .map(|&m| if m > 0.5 { 0.0 } else { inf })
        .collect();

    // Phase 1: 1D EDT along rows (rows are contiguous in row-major layout).
    let mut row_pass = vec![0.0f64; init.len()];
    row_pass
        .par_chunks_mut(ncols)
        .zip(init.par_chunks(ncols))
        .for_each(|(out, inp)| dt_1d(inp, out, inf));

    // Phase 2: 1D EDT along columns. Transpose first so columns become
    // contiguous rows, run DT, transpose back.
    let transposed = transpose(&row_pass, nrows, ncols);
    let mut col_pass = vec![0.0f64; transposed.len()];
    col_pass
        .par_chunks_mut(nrows)
        .zip(transposed.par_chunks(nrows))
        .for_each(|(out, inp)| dt_1d(inp, out, inf));

    transpose(&col_pass, ncols, nrows)
}

/// Exact Euclidean distance transform of a binary mask.
///
/// Returns the Euclidean distance (not squared) from each cell to the
/// nearest source cell, in cell units. R-callable entry point.
///
/// @param mask Numeric vector, row-major, length nrows*ncols. Source cells
///   have value > 0.5; background cells have value <= 0.5 (or NA, but the
///   R wrapper should replace NAs first).
/// @param nrows Integer. Number of raster rows.
/// @param ncols Integer. Number of raster columns.
/// @param n_cores Integer. Number of threads to use. 0 means "use rayon's
///   default" (all logical cores). Positive values cap parallelism at that
///   thread count; 1 runs serially. Negative values are treated as 0.
/// @return Numeric vector of distances, row-major, same length as `mask`.
#[extendr]
fn rust_edt(mask: Vec<f64>, nrows: i32, ncols: i32, n_cores: i32) -> Vec<f64> {
    let nr = nrows as usize;
    let nc = ncols as usize;
    let n_threads = n_cores.max(0) as usize;

    let compute = || {
        let sq = edt_squared(&mask, nr, nc);
        sq.into_par_iter().map(f64::sqrt).collect::<Vec<f64>>()
    };

    if n_threads == 0 {
        // Use rayon's global thread pool (default: all logical cores).
        compute()
    } else {
        // Build a local thread pool capped at n_threads and run the EDT
        // inside it. install() routes all par_iter calls in this closure
        // through the local pool.
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .expect("failed to build rayon thread pool");
        pool.install(compute)
    }
}

// extendr macro: generates R_init_howfar_extendr and registers
// `rust_edt` as a callable native routine.
extendr_module! {
    mod howfar;
    fn rust_edt;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn dt_1d_single_source() {
        let inf = 100.0;
        let f = vec![inf, 0.0, inf, inf, inf];
        let mut d = vec![0.0; 5];
        dt_1d(&f, &mut d, inf);
        // squared distances from index 1: [1, 0, 1, 4, 9]
        assert_eq!(d, vec![1.0, 0.0, 1.0, 4.0, 9.0]);
    }

    #[test]
    fn dt_1d_two_sources() {
        let inf = 100.0;
        let f = vec![0.0, inf, inf, 0.0, inf];
        let mut d = vec![0.0; 5];
        dt_1d(&f, &mut d, inf);
        // squared distances: [0, 1, 1, 0, 1]
        assert_eq!(d, vec![0.0, 1.0, 1.0, 0.0, 1.0]);
    }

    #[test]
    fn dt_1d_all_background() {
        let inf = 100.0;
        let f = vec![inf, inf, inf];
        let mut d = vec![0.0; 3];
        dt_1d(&f, &mut d, inf);
        assert_eq!(d, vec![inf, inf, inf]);
    }

    #[test]
    fn edt_2d_single_source() {
        // 5x5 grid with single source at row 2, col 2 (0-indexed center).
        let mut mask = vec![0.0; 25];
        mask[2 * 5 + 2] = 1.0;
        let sq = edt_squared(&mask, 5, 5);

        // Expected squared distances: (r-2)^2 + (c-2)^2 for each (r, c).
        for r in 0..5 {
            for c in 0..5 {
                let expected = ((r as f64 - 2.0).powi(2)) + ((c as f64 - 2.0).powi(2));
                let got = sq[r * 5 + c];
                assert!(
                    approx_eq(got, expected, 1e-9),
                    "at (r={}, c={}): expected {}, got {}", r, c, expected, got
                );
            }
        }
    }

    #[test]
    fn edt_2d_vertical_line() {
        // 7x11 grid with a vertical line along column 5.
        let nrows = 7;
        let ncols = 11;
        let mut mask = vec![0.0; nrows * ncols];
        for r in 0..nrows {
            mask[r * ncols + 5] = 1.0;
        }
        let sq = edt_squared(&mask, nrows, ncols);

        // Squared distance to vertical line at c=5 is (c-5)^2 (independent of row).
        for r in 0..nrows {
            for c in 0..ncols {
                let expected = (c as f64 - 5.0).powi(2);
                let got = sq[r * ncols + c];
                assert!(
                    approx_eq(got, expected, 1e-9),
                    "at (r={}, c={}): expected {}, got {}", r, c, expected, got
                );
            }
        }
    }
}
