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

/// Like [`dt_1d`], but also records, for each output index `q`, the apex
/// (source location) of the parabola that wins at `q` in `arg[q]`. When the
/// input slice has no source (all values `>= inf`) the distances are set to
/// `inf` and `arg` is left untouched — callers treat an `inf` distance as
/// "no source" and ignore the corresponding `arg`.
fn dt_1d_arg(f: &[f64], d: &mut [f64], arg: &mut [usize], inf: f64) {
    let n = f.len();
    debug_assert_eq!(d.len(), n);
    debug_assert_eq!(arg.len(), n);
    if n == 0 {
        return;
    }

    let mut first = 0usize;
    while first < n && f[first] >= inf {
        first += 1;
    }

    if first == n {
        for x in d.iter_mut() {
            *x = inf;
        }
        return;
    }

    let mut v: Vec<usize> = vec![0; n];
    let mut z: Vec<f64> = vec![0.0; n + 1];

    v[0] = first;
    z[0] = f64::NEG_INFINITY;
    z[1] = f64::INFINITY;
    let mut k: usize = 0;

    for q in (first + 1)..n {
        let fq = f[q];
        if fq >= inf {
            continue;
        }
        let qf = q as f64;
        let q_sq = qf * qf;

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

    let mut kr: usize = 0;
    for q in 0..n {
        let qf = q as f64;
        while z[kr + 1] < qf {
            kr += 1;
        }
        let vk = v[kr];
        let dq = qf - vk as f64;
        d[q] = dq * dq + f[vk];
        arg[q] = vk;
    }
}

/// Transpose a row-major matrix of shape (in_rows × in_cols) into a
/// row-major matrix of shape (in_cols × in_rows) of `T`. Parallel by output row.
fn transpose<T: Copy + Default + Send + Sync>(input: &[T], in_rows: usize, in_cols: usize) -> Vec<T> {
    let mut out = vec![T::default(); input.len()];
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

/// Like [`edt_squared`], but also returns the *feature transform*: for each
/// cell, the row-major linear index of the nearest source cell (the same
/// planar-Euclidean nearest used for the distance). Cells with no reachable
/// source (only when the mask has zero sources) get index `usize::MAX` and
/// squared distance `inf`.
///
/// Returns `(squared_distance, nearest_source_index)`, both row-major.
fn edt_squared_ft(mask: &[f64], nrows: usize, ncols: usize) -> (Vec<f64>, Vec<usize>) {
    assert_eq!(mask.len(), nrows * ncols, "mask size mismatch");

    let max_dim = nrows.max(ncols) as f64;
    let inf = 4.0 * max_dim * max_dim + 4.0;

    let init: Vec<f64> = mask
        .par_iter()
        .map(|&m| if m > 0.5 { 0.0 } else { inf })
        .collect();

    // Phase 1 (rows): row_d2[r,c] = squared distance to the nearest source in
    // row r; row_arg[r,c] = the column of that source.
    let mut row_d2 = vec![0.0f64; init.len()];
    let mut row_arg = vec![0usize; init.len()];
    row_d2
        .par_chunks_mut(ncols)
        .zip(row_arg.par_chunks_mut(ncols))
        .zip(init.par_chunks(ncols))
        .for_each(|((dout, aout), inp)| dt_1d_arg(inp, dout, aout, inf));

    // Phase 2 (columns): transpose so columns are contiguous, run the 1D DT on
    // the row-pass distances. col_arg_t holds the source *row* r' that wins.
    let row_d2_t = transpose(&row_d2, nrows, ncols);
    let mut col_d2_t = vec![0.0f64; row_d2_t.len()];
    let mut col_arg_t = vec![0usize; row_d2_t.len()];
    col_d2_t
        .par_chunks_mut(nrows)
        .zip(col_arg_t.par_chunks_mut(nrows))
        .zip(row_d2_t.par_chunks(nrows))
        .for_each(|((dout, aout), inp)| dt_1d_arg(inp, dout, aout, inf));

    let d2 = transpose(&col_d2_t, ncols, nrows);
    let col_arg = transpose(&col_arg_t, ncols, nrows); // col_arg[r,c] = source row r'

    // Combine: the nearest source is in row r' = col_arg[r,c], and within that
    // row its column is row_arg[r', c]. (`d2 >= inf` only when the whole mask
    // has no source.)
    let feat: Vec<usize> = (0..nrows * ncols)
        .into_par_iter()
        .map(|i| {
            if d2[i] >= inf {
                usize::MAX
            } else {
                let c = i % ncols;
                let rp = col_arg[i];
                let cp = row_arg[rp * ncols + c];
                rp * ncols + cp
            }
        })
        .collect();

    (d2, feat)
}

/// Euclidean distance between two points in a projected coordinate system.
/// `hypot` avoids intermediate overflow/underflow while staying accurate.
fn euclidean(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
    let dx = x2 - x1;
    let dy = y2 - y1;
    dx.hypot(dy)
}

/// Haversine (great-circle) distance in metres between two longitude/latitude
/// points (in degrees) on a sphere of mean Earth radius.
fn haversine(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let r = 6_371_000.0; // mean Earth radius in metres
    let (lat1, lon1) = (lat1.to_radians(), lon1.to_radians());
    let (lat2, lon2) = (lat2.to_radians(), lon2.to_radians());

    let dlat = lat2 - lat1;
    let dlon = lon2 - lon1;

    let a = (dlat / 2.0).sin().powi(2)
        + lat1.cos() * lat2.cos() * (dlon / 2.0).sin().powi(2);
    let c = 2.0 * a.sqrt().atan2((1.0 - a).sqrt());

    r * c
}

/// Convert a longitude/latitude (degrees) to a unit vector on the sphere.
/// Nearest neighbour by 3D chord distance is exactly nearest neighbour by
/// great-circle distance (chord length is monotonic in the central angle), so
/// a Euclidean nearest-neighbour search on these vectors gives the truly
/// closest source by Haversine distance — at any latitude, across the
/// antimeridian, and over the poles.
fn geo_to_unit(lon_deg: f64, lat_deg: f64) -> [f64; 3] {
    let lat = lat_deg.to_radians();
    let lon = lon_deg.to_radians();
    let cos_lat = lat.cos();
    [cos_lat * lon.cos(), cos_lat * lon.sin(), lat.sin()]
}

/// Squared Euclidean distance between two 3D points.
fn sq_dist3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

const KD_NONE: u32 = u32::MAX;

/// A node of [`KdTree3`].
#[derive(Clone, Copy)]
struct KdNode {
    point: [f64; 3],
    payload: u32,
    left: u32,
    right: u32,
    axis: u8,
}

/// Minimal balanced 3D k-d tree for exact nearest-neighbour queries.
///
/// Built once; `nearest` is read-only (and `Sync`), so queries run in
/// parallel. Used to find, for each cell, the source cell closest on the unit
/// sphere — i.e. the source nearest in great-circle distance.
struct KdTree3 {
    nodes: Vec<KdNode>,
    root: u32,
}

impl KdTree3 {
    fn build(mut items: Vec<([f64; 3], u32)>) -> Option<KdTree3> {
        if items.is_empty() {
            return None;
        }
        let mut nodes = Vec::with_capacity(items.len());
        let root = kd_build(&mut items[..], 0, &mut nodes);
        Some(KdTree3 { nodes, root })
    }

    /// Returns `(payload, squared_distance)` of the nearest stored point.
    fn nearest(&self, q: &[f64; 3]) -> (u32, f64) {
        let mut best = (KD_NONE, f64::INFINITY);
        self.search(self.root, q, &mut best);
        best
    }

    fn search(&self, node: u32, q: &[f64; 3], best: &mut (u32, f64)) {
        if node == KD_NONE {
            return;
        }
        let nd = &self.nodes[node as usize];
        let d2 = sq_dist3(&nd.point, q);
        if d2 < best.1 {
            *best = (nd.payload, d2);
        }
        let axis = nd.axis as usize;
        let diff = q[axis] - nd.point[axis];
        let (near, far) = if diff < 0.0 {
            (nd.left, nd.right)
        } else {
            (nd.right, nd.left)
        };
        self.search(near, q, best);
        // Only descend the far branch if it could hold a closer point.
        if diff * diff < best.1 {
            self.search(far, q, best);
        }
    }
}

/// Recursively build a balanced subtree from `items` by median-splitting on a
/// cycling axis; appends nodes to `nodes` and returns this subtree's root index.
fn kd_build(items: &mut [([f64; 3], u32)], depth: usize, nodes: &mut Vec<KdNode>) -> u32 {
    let axis = depth % 3;
    let mid = items.len() / 2;
    items.select_nth_unstable_by(mid, |a, b| {
        a.0[axis]
            .partial_cmp(&b.0[axis])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let node_idx = nodes.len() as u32;
    nodes.push(KdNode {
        point: items[mid].0,
        payload: items[mid].1,
        left: KD_NONE,
        right: KD_NONE,
        axis: axis as u8,
    });

    let (left, right_incl) = items.split_at_mut(mid);
    let right = &mut right_incl[1..]; // exclude the median itself

    if !left.is_empty() {
        let lc = kd_build(left, depth + 1, nodes);
        nodes[node_idx as usize].left = lc;
    }
    if !right.is_empty() {
        let rc = kd_build(right, depth + 1, nodes);
        nodes[node_idx as usize].right = rc;
    }
    node_idx
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

/// Exact Euclidean distance transform returning distance in **metres**.
///
/// Finds, for each cell, the planar-Euclidean nearest source cell (via the
/// feature transform), then measures the metric distance between cell centres:
/// `haversine` when `lonlat` is true, otherwise `euclidean` scaled by
/// `m_per_unit` (metres per CRS linear unit). Cell centres are derived from the
/// grid geometry: cell `(r, c)` has centre `(xmin + (c + 0.5) * xres,
/// ymax - (r + 0.5) * yres)`. Cells with no source return infinity.
///
/// @param mask Numeric vector, row-major, length nrows*ncols (>0.5 = source).
/// @param nrows Integer. Number of raster rows.
/// @param ncols Integer. Number of raster columns.
/// @param n_cores Integer. Threads to use; 0 = rayon default, 1 = serial.
/// @param xmin Double. Left edge (x) of the raster extent.
/// @param ymax Double. Top edge (y) of the raster extent.
/// @param xres Double. Cell size in x.
/// @param yres Double. Cell size in y.
/// @param lonlat Logical. If true, use Haversine (lon/lat in degrees).
/// @param m_per_unit Double. Metres per CRS linear unit (projected only).
/// @return Numeric vector of distances in metres, row-major.
#[extendr]
fn rust_edt_meters(
    mask: Vec<f64>,
    nrows: i32,
    ncols: i32,
    n_cores: i32,
    xmin: f64,
    ymax: f64,
    xres: f64,
    yres: f64,
    lonlat: bool,
    m_per_unit: f64,
) -> Vec<f64> {
    let nr = nrows as usize;
    let nc = ncols as usize;
    let n_threads = n_cores.max(0) as usize;

    let compute = || {
        let (_d2, feat) = edt_squared_ft(&mask, nr, nc);

        // Cell-centre coordinates from a row-major linear index.
        let cell_centre = |idx: usize| -> (f64, f64) {
            let r = idx / nc;
            let c = idx % nc;
            (
                xmin + (c as f64 + 0.5) * xres,
                ymax - (r as f64 + 0.5) * yres,
            )
        };

        (0..nr * nc)
            .into_par_iter()
            .map(|i| {
                let src = feat[i];
                if src == usize::MAX {
                    f64::INFINITY // no source anywhere in the mask
                } else if src == i {
                    0.0 // this cell is itself a source
                } else {
                    let (xq, yq) = cell_centre(i);
                    let (xs, ys) = cell_centre(src);
                    if lonlat {
                        haversine(xq, yq, xs, ys)
                    } else {
                        euclidean(xq, yq, xs, ys) * m_per_unit
                    }
                }
            })
            .collect::<Vec<f64>>()
    };

    if n_threads == 0 {
        compute()
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .expect("failed to build rayon thread pool");
        pool.install(compute)
    }
}

/// Exact great-circle (Haversine) distance to the nearest source, in metres.
///
/// Unlike the grid distance transform, this finds the source cell that is
/// genuinely nearest in great-circle distance — via a 3D nearest-neighbour
/// query over the source cells mapped onto the unit sphere — and returns the
/// Haversine distance to it. Correct at any latitude and across the
/// antimeridian; the only remaining approximation is the rasterisation of
/// features to cell centres.
///
/// Cell centres are lon/lat (degrees) from the grid geometry, using the same
/// indexing convention as `rust_edt_meters`. Cells with no source return
/// infinity.
///
/// @param mask Numeric vector, row-major, length nrows*ncols (>0.5 = source).
/// @param nrows Integer. Number of raster rows.
/// @param ncols Integer. Number of raster columns.
/// @param n_cores Integer. Threads to use; 0 = rayon default, 1 = serial.
/// @param xmin Double. Left edge (longitude) of the extent.
/// @param ymax Double. Top edge (latitude) of the extent.
/// @param xres Double. Cell size in longitude (degrees).
/// @param yres Double. Cell size in latitude (degrees).
/// @return Numeric vector of great-circle distances in metres, row-major.
#[extendr]
fn rust_geo_distance(
    mask: Vec<f64>,
    nrows: i32,
    ncols: i32,
    n_cores: i32,
    xmin: f64,
    ymax: f64,
    xres: f64,
    yres: f64,
) -> Vec<f64> {
    let nr = nrows as usize;
    let nc = ncols as usize;
    let n_threads = n_cores.max(0) as usize;

    let compute = || {
        let cell_lonlat = |idx: usize| -> (f64, f64) {
            let r = idx / nc;
            let c = idx % nc;
            (
                xmin + (c as f64 + 0.5) * xres,
                ymax - (r as f64 + 0.5) * yres,
            )
        };

        // Build a 3D nearest-neighbour index over the source cells.
        let mut items: Vec<([f64; 3], u32)> = Vec::new();
        for i in 0..nr * nc {
            if mask[i] > 0.5 {
                let (lon, lat) = cell_lonlat(i);
                items.push((geo_to_unit(lon, lat), i as u32));
            }
        }
        let tree = KdTree3::build(items);

        (0..nr * nc)
            .into_par_iter()
            .map(|i| {
                if mask[i] > 0.5 {
                    return 0.0; // this cell is a source
                }
                match &tree {
                    None => f64::INFINITY, // no source anywhere
                    Some(t) => {
                        let (lon, lat) = cell_lonlat(i);
                        let (src, _d2) = t.nearest(&geo_to_unit(lon, lat));
                        let (slon, slat) = cell_lonlat(src as usize);
                        haversine(lon, lat, slon, slat)
                    }
                }
            })
            .collect::<Vec<f64>>()
    };

    if n_threads == 0 {
        compute()
    } else {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .expect("failed to build rayon thread pool");
        pool.install(compute)
    }
}

// extendr macro: generates R_init_howfar_extendr and registers the native
// routines callable from R.
extendr_module! {
    mod howfar;
    fn rust_edt;
    fn rust_edt_meters;
    fn rust_geo_distance;
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

    #[test]
    fn edt_ft_single_source() {
        // 5x5 grid, single source at (2, 2): every cell's nearest is that cell.
        let mut mask = vec![0.0; 25];
        mask[2 * 5 + 2] = 1.0;
        let (d2, feat) = edt_squared_ft(&mask, 5, 5);

        for r in 0..5 {
            for c in 0..5 {
                let i = r * 5 + c;
                let expected = (r as f64 - 2.0).powi(2) + (c as f64 - 2.0).powi(2);
                assert!(approx_eq(d2[i], expected, 1e-9));
                assert_eq!(feat[i], 2 * 5 + 2, "nearest source is the only source");
            }
        }
    }

    #[test]
    fn edt_ft_vertical_line() {
        // Vertical line at column 5: nearest source for (r, c) is (r, 5),
        // with no horizontal ambiguity.
        let nrows = 7;
        let ncols = 11;
        let mut mask = vec![0.0; nrows * ncols];
        for r in 0..nrows {
            mask[r * ncols + 5] = 1.0;
        }
        let (d2, feat) = edt_squared_ft(&mask, nrows, ncols);

        for r in 0..nrows {
            for c in 0..ncols {
                let i = r * ncols + c;
                assert!(approx_eq(d2[i], (c as f64 - 5.0).powi(2), 1e-9));
                assert_eq!(feat[i], r * ncols + 5, "nearest is (r, 5)");
            }
        }
    }

    #[test]
    fn edt_ft_no_source() {
        let mask = vec![0.0; 9];
        let (_d2, feat) = edt_squared_ft(&mask, 3, 3);
        assert!(feat.iter().all(|&f| f == usize::MAX));
    }

    #[test]
    fn euclidean_3_4_5() {
        assert!(approx_eq(euclidean(0.0, 0.0, 3.0, 4.0), 5.0, 1e-12));
        assert!(approx_eq(euclidean(1.0, 1.0, 1.0, 1.0), 0.0, 1e-12));
    }

    #[test]
    fn haversine_one_degree_latitude() {
        // One degree of latitude on the mean-radius sphere.
        let expected = 6_371_000.0 * std::f64::consts::PI / 180.0;
        assert!(approx_eq(haversine(0.0, 0.0, 0.0, 1.0), expected, 1e-6));
        assert!(approx_eq(haversine(10.0, 20.0, 10.0, 20.0), 0.0, 1e-9));
    }

    #[test]
    fn kdtree_matches_brute_force() {
        // xorshift64* PRNG for deterministic pseudo-random points in [0, 1).
        let mut s: u64 = 0x9E3779B97F4A7C15;
        let mut rng = || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let pts: Vec<([f64; 3], u32)> =
            (0..300).map(|i| ([rng(), rng(), rng()], i as u32)).collect();
        let tree = KdTree3::build(pts.clone()).unwrap();

        for _ in 0..100 {
            let q = [rng(), rng(), rng()];
            let (_payload, d2) = tree.nearest(&q);
            let brute = pts
                .iter()
                .map(|(p, _)| sq_dist3(p, &q))
                .fold(f64::INFINITY, f64::min);
            assert!(approx_eq(d2, brute, 1e-12), "kd {} vs brute {}", d2, brute);
        }
    }

    #[test]
    fn geo_nearest_respects_great_circle() {
        // At latitude 80, degrees of longitude are compressed. For a cell at
        // (lon=0, lat=80):
        //   A = (5, 80): 5 deg of longitude -> ~96.5 km on the ground
        //   B = (0, 78): 2 deg of latitude  -> ~222 km
        // In raw lon/lat degrees B looks nearer (2 < 5), but on the sphere A
        // wins — which a great-circle nearest-neighbour search must reflect.
        let cell = geo_to_unit(0.0, 80.0);
        let tree = KdTree3::build(vec![
            (geo_to_unit(5.0, 80.0), 0), // A
            (geo_to_unit(0.0, 78.0), 1), // B
        ])
        .unwrap();
        let (winner, _) = tree.nearest(&cell);
        assert_eq!(winner, 0, "great-circle nearest should be A, not B");

        let d_a = haversine(0.0, 80.0, 5.0, 80.0);
        let d_b = haversine(0.0, 80.0, 0.0, 78.0);
        assert!(d_a < d_b);
        assert!((d_a - 96_515.0).abs() < 100.0, "haversine to A ~96.5km, got {}", d_a);
    }
}
