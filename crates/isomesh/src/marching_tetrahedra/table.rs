//! The six-tetrahedron decomposition of a cube, and the 16-case table for one
//! tetrahedron — both **derived at compile time**.
//!
//! # Why this is generated rather than transcribed
//!
//! Same reason as [`crate::marching_cubes::table`], and here the case is
//! stronger: the original method is Doi & Koide, *An efficient method of
//! triangulating equi-valued surfaces by using tetrahedral cells*, IEICE Trans.
//! E74-D(1), 214–224 (1991) — which has **no DOI, is not in the local corpus and
//! could not be obtained**. There is nothing to transcribe from even if
//! transcription were wanted, and inventing a table is what rule 5 forbids.
//!
//! So both structures are constructed, and then checked. Neither is a lookup.
//!
//! # The decomposition
//!
//! **Kuhn's**, also called Freudenthal's: the six tetrahedra are the six
//! *monotone paths* from corner `0` at `(0,0,0)` to corner `7` at `(1,1,1)`
//! along cube edges, one per ordering of the three axes. Every tetrahedron
//! therefore contains the main diagonal `0–7`, and the six tile the cube exactly.
//!
//! ```text
//! axes (x,y,z) -> 0,1,3,7      axes (y,x,z) -> 0,2,3,7      axes (z,x,y) -> 0,4,5,7
//! axes (x,z,y) -> 0,1,5,7      axes (y,z,x) -> 0,2,6,7      axes (z,y,x) -> 0,4,6,7
//! ```
//!
//! # Why *this* decomposition, and not the five-tetrahedron one
//!
//! Five tetrahedra also tile a cube and are cheaper, but the five-tet tiling is
//! **chiral**: neighbouring cells must alternate its handedness like a
//! checkerboard or their shared face is split by two different diagonals and the
//! surfaces do not meet. That is a crack, and it is the same failure `✗11` is
//! about.
//!
//! Kuhn's needs no alternation, and the reason is worth stating because it is
//! the property that makes this safe rather than a convention that happens to
//! work. Each cube face carries exactly two of the tetrahedra's triangles, split
//! by one diagonal:
//!
//! ```text
//! −x face: 0–6     +x face: 1–7
//! −y face: 0–5     +y face: 2–7
//! −z face: 0–3     +z face: 4–7
//! ```
//!
//! Take two cells adjacent along `x`. The first splits its `+x` face on `1–7`,
//! which in world terms is the segment from `origin + (1,0,0)` to
//! `origin + (1,1,1)`. The second splits its `−x` face on `0–6`, which is
//! `origin' + (0,0,0)` to `origin' + (0,1,1)` — and `origin' = origin + (1,0,0)`,
//! so those are **the same two points**. The same holds on `y` and `z`. Every
//! cell in the lattice uses the same corner numbering, so every shared face is
//! split the same way by both of its cells, with nothing to coordinate.
//! `every_shared_face_is_split_the_same_way_by_both_cells` checks all three axes.
//!
//! # The case table
//!
//! A tetrahedron has four corners, so sixteen sign configurations, and **no
//! ambiguity**: the linear interpolant on a tetrahedron is determined by its
//! four values, and there is no face saddle and no body saddle to decide. That
//! is the whole reason this algorithm is the topological reference the others
//! are compared against.
//!
//! One inside corner cuts the three edges meeting it and gives a triangle; two
//! inside corners cut four edges and give a quad, fanned into two triangles;
//! three inside is one-inside complemented. Winding is not guessed either — it
//! is computed from the corner coordinates at compile time, in integers, and
//! flipped where the cross product points the wrong way.

use crate::cube::corner_offset;

/// A cube splits into six tetrahedra.
pub const TET_COUNT: usize = 6;

/// A tetrahedron has six edges.
pub const TET_EDGE_COUNT: usize = 6;

/// Most triangles one tetrahedron can contribute: the two-inside quad.
pub const MAX_TET_TRIANGLES: usize = 2;

/// The six tetrahedra, as cube-corner indices.
///
/// Built by walking the axes in each of the six orders and accumulating bits, so
/// the corners of every tetrahedron are ordered by inclusion: `TETS[t][0]` is
/// always corner `0` and `TETS[t][3]` always corner `7`. That ordering is relied
/// on by [`tet_edge_corners`] — every edge runs from fewer bits to more, so its
/// two corner offsets differ by a `0/1` step on each axis and never by `-1`.
pub static TETS: [[u8; 4]; TET_COUNT] = build_tets();

const fn build_tets() -> [[u8; 4]; TET_COUNT] {
    // The six orderings of three axes, written out because `const fn` has no
    // iterator to permute with. Checked against a generated permutation in
    // `the_six_tetrahedra_are_the_six_axis_orderings`.
    let orders = [
        [0u8, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut out = [[0u8; 4]; TET_COUNT];
    let mut t = 0usize;
    while t < TET_COUNT {
        let mut corner = 0u8;
        out[t][0] = 0;
        let mut step = 0usize;
        while step < 3 {
            corner |= 1 << orders[t][step];
            out[t][step + 1] = corner;
            step += 1;
        }
        t += 1;
    }
    out
}

/// The two corners of each edge of a tetrahedron, as indices into `TETS[t]`.
///
/// Lower index first, which given `TETS`' inclusion ordering also means fewer
/// bits first.
pub const TET_EDGES: [[u8; 2]; TET_EDGE_COUNT] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];

/// The cube corners an edge of tetrahedron `t` joins.
#[must_use]
pub const fn tet_edge_corners(t: usize, edge: usize) -> [u8; 2] {
    let [a, b] = TET_EDGES[edge];
    [TETS[t][a as usize], TETS[t][b as usize]]
}

/// The triangulation of one tetrahedron for one sign configuration.
#[derive(Clone, Copy, Debug)]
pub struct TetCase {
    /// Number of triangles.
    pub count: u8,
    /// Each triangle as three tetrahedron-edge indices.
    pub triangles: [[u8; 3]; MAX_TET_TRIANGLES],
}

/// All sixteen configurations for each of the six tetrahedra.
///
/// Indexed `[tetrahedron][case]`, where the case is a four-bit mask over
/// `TETS[t]` — bit `i` set means `TETS[t][i]` is inside. It depends on the
/// tetrahedron as well as the case because the **winding** does: two tetrahedra
/// with the same sign pattern sit at different handedness inside the cube.
pub static TET_CASES: [[TetCase; 16]; TET_COUNT] = build_tet_cases();

const fn build_tet_cases() -> [[TetCase; 16]; TET_COUNT] {
    let empty = TetCase {
        count: 0,
        triangles: [[0u8; 3]; MAX_TET_TRIANGLES],
    };
    let mut out = [[empty; 16]; TET_COUNT];
    let mut t = 0usize;
    while t < TET_COUNT {
        let mut case = 0usize;
        while case < 16 {
            out[t][case] = build_tet_case(t, case as u8);
            case += 1;
        }
        t += 1;
    }
    out
}

/// Doubled coordinates of the midpoint of a tetrahedron edge.
///
/// Doubled so the midpoint of two `0/1` corners is an integer and the whole
/// orientation test can run in `i32` at compile time. A float there would be
/// exact for these values anyway, but integers make that obvious rather than
/// something to check.
const fn edge_midpoint_x2(t: usize, edge: usize) -> [i32; 3] {
    let [a, b] = tet_edge_corners(t, edge);
    let pa = corner_offset(a);
    let pb = corner_offset(b);
    [
        (pa[0] + pb[0]) as i32,
        (pa[1] + pb[1]) as i32,
        (pa[2] + pb[2]) as i32,
    ]
}

const fn build_tet_case(t: usize, case: u8) -> TetCase {
    let mut triangles = [[0u8; 3]; MAX_TET_TRIANGLES];
    let mut count = 0usize;

    // Which edges are cut: one end inside, the other outside.
    let mut cut = [false; TET_EDGE_COUNT];
    let mut cut_count = 0usize;
    let mut e = 0usize;
    while e < TET_EDGE_COUNT {
        let [a, b] = TET_EDGES[e];
        let a_in = case & (1 << a) != 0;
        let b_in = case & (1 << b) != 0;
        cut[e] = a_in != b_in;
        if cut[e] {
            cut_count += 1;
        }
        e += 1;
    }

    // Three cut edges is one triangle; four is a quad. Nothing else can occur on
    // a tetrahedron, which is the absence of ambiguity made concrete.
    if cut_count == 3 {
        let mut tri = [0u8; 3];
        let mut n = 0usize;
        let mut e = 0usize;
        while e < TET_EDGE_COUNT {
            if cut[e] {
                tri[n] = e as u8;
                n += 1;
            }
            e += 1;
        }
        triangles[0] = orient(t, tri, case);
        count = 1;
    } else if cut_count == 4 {
        let mut quad = [0u8; 4];
        let mut n = 0usize;
        let mut e = 0usize;
        while e < TET_EDGE_COUNT {
            if cut[e] {
                quad[n] = e as u8;
                n += 1;
            }
            e += 1;
        }
        // The four cut edges in index order are not a cycle round the quad: two
        // of them are diagonally opposite. Reorder so consecutive entries share
        // a tetrahedron corner, which is what makes the fan a quad and not a
        // bowtie.
        let quad = order_quad(quad);
        triangles[0] = orient(t, [quad[0], quad[1], quad[2]], case);
        triangles[1] = orient(t, [quad[0], quad[2], quad[3]], case);
        count = 2;
    }

    TetCase {
        count: count as u8,
        triangles,
    }
}

/// Do these two tetrahedron edges meet at a corner?
///
/// A free function rather than a closure because `const fn` cannot call one.
const fn edges_share_a_corner(x: u8, y: u8) -> bool {
    let [xa, xb] = TET_EDGES[x as usize];
    let [ya, yb] = TET_EDGES[y as usize];
    xa == ya || xa == yb || xb == ya || xb == yb
}

/// Put four cut edges into an order that walks the quad's boundary.
///
/// Two cut edges are adjacent on the quad when they share an endpoint *corner*
/// of the tetrahedron; the two that share neither are the diagonal. Starting
/// anywhere and repeatedly taking an unused edge that shares a corner with the
/// last one traverses the ring.
const fn order_quad(quad: [u8; 4]) -> [u8; 4] {
    let mut out = [quad[0], 0, 0, 0];
    let mut used = [true, false, false, false];
    let mut n = 1usize;
    while n < 4 {
        let mut i = 1usize;
        while i < 4 {
            if !used[i] && edges_share_a_corner(out[n - 1], quad[i]) {
                out[n] = quad[i];
                used[i] = true;
                break;
            }
            i += 1;
        }
        n += 1;
    }
    out
}

/// Wind a triangle so its normal points **away from the solid**.
///
/// Computed, not tabulated: take the cross product of the triangle as built,
/// dot it with a vector pointing from an inside corner toward the triangle, and
/// swap two vertices if the sign is wrong. All in doubled integer coordinates,
/// so it is exact and runs at compile time.
const fn orient(t: usize, tri: [u8; 3], case: u8) -> [u8; 3] {
    let p0 = edge_midpoint_x2(t, tri[0] as usize);
    let p1 = edge_midpoint_x2(t, tri[1] as usize);
    let p2 = edge_midpoint_x2(t, tri[2] as usize);

    let u = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let v = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let normal = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];

    // Any inside corner will do: the triangle separates the inside corners from
    // the outside ones, so "away from an inside corner" is unambiguous.
    let mut inside = 0usize;
    while inside < 4 {
        if case & (1 << inside) != 0 {
            break;
        }
        inside += 1;
    }
    let c = corner_offset(TETS[t][inside]);
    let from_inside = [
        p0[0] - 2 * c[0] as i32,
        p0[1] - 2 * c[1] as i32,
        p0[2] - 2 * c[2] as i32,
    ];
    let facing =
        normal[0] * from_inside[0] + normal[1] * from_inside[1] + normal[2] * from_inside[2];

    if facing >= 0 {
        tri
    } else {
        [tri[0], tri[2], tri[1]]
    }
}
