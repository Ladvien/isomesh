//! The exact convex wedge, shared by the experiments that need a crease.
//!
//! Ticket: R-008, extracted from R-006's `experiment_p13` so there is one
//! wedge rather than two. A fixture copied between benches is a fixture that
//! can drift, and R-008's whole question is about geometry R-006 measured — the
//! two have to be looking at the same object.

use isomesh::Sdf;

/// The scalar the wedge experiments run in.
///
/// `f64`, because the quantity being measured is an angle at a discontinuity
/// and the question is whether it moves in the third decimal.
pub(crate) type Scalar = f64;

/// An exact convex wedge: the set of points within `±θ/2` of the `+x` axis in
/// the `xy` plane, extruded along `z`.
///
/// # Why this is exact and `max(d₁, d₂)` is not
///
/// Inside, the distance to a convex region is the distance to the nearest
/// bounding plane, which is `max(d₁, d₂)` — correct. Outside, a point beyond
/// both planes is nearest the **edge**, and `max` returns the distance to the
/// further plane instead, which is smaller. Phase 11 names that object: a
/// Pseudo-SDF, eikonal almost everywhere and wrong at the seam. Here the seam is
/// the whole experiment, so the exterior is computed as the distance to the two
/// bounding rays, clamped at the apex.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Wedge {
    /// Half the dihedral, in radians.
    pub(crate) half: Scalar,
    /// World position of the apex; the crease is the line `x = apex.x`,
    /// `y = apex.y` running along `z`.
    pub(crate) apex: [Scalar; 2],
    /// Rotation of the whole wedge about the crease, in radians.
    pub(crate) rotation: Scalar,
}

impl Wedge {
    pub(crate) fn new(dihedral_deg: f64, apex: [Scalar; 2], rotation_deg: f64) -> Self {
        Self {
            half: dihedral_deg.to_radians() / 2.0,
            apex,
            rotation: rotation_deg.to_radians(),
        }
    }

    /// Outward normals of the two bounding half-planes, in the `xy` plane.
    pub(crate) fn plane_normals(&self) -> [[Scalar; 2]; 2] {
        let (s, c) = (self.half.sin(), self.half.cos());
        [[-s, c], [-s, -c]]
    }

    /// Directions of the two bounding rays, from the apex.
    fn ray_directions(&self) -> [[Scalar; 2]; 2] {
        let (s, c) = (self.half.sin(), self.half.cos());
        [[c, s], [c, -s]]
    }

    /// `p` relative to the apex and rotated into the wedge's own frame.
    pub(crate) fn local(&self, p: [Scalar; 3]) -> [Scalar; 2] {
        let (x, y) = (p[0] - self.apex[0], p[1] - self.apex[1]);
        let (s, c) = (self.rotation.sin(), self.rotation.cos());
        [x * c + y * s, -x * s + y * c]
    }

    /// A direction in the wedge's frame, back in world coordinates.
    pub(crate) fn unrotate(&self, v: [Scalar; 2]) -> [Scalar; 3] {
        let (s, c) = (self.rotation.sin(), self.rotation.cos());
        [v[0] * c - v[1] * s, v[0] * s + v[1] * c, 0.0]
    }

    /// The nearest point on ray `dir` to `q`, and the vector from it to `q`.
    fn to_ray(q: [Scalar; 2], dir: [Scalar; 2]) -> ([Scalar; 2], Scalar) {
        let t = (q[0] * dir[0] + q[1] * dir[1]).max(0.0);
        let away = [q[0] - dir[0] * t, q[1] - dir[1] * t];
        (away, (away[0] * away[0] + away[1] * away[1]).sqrt())
    }
}

impl Sdf for Wedge {
    type Scalar = Scalar;

    fn sample(&self, p: [Scalar; 3]) -> Scalar {
        let q = self.local(p);
        let [n0, n1] = self.plane_normals();
        let d0 = q[0] * n0[0] + q[1] * n0[1];
        let d1 = q[0] * n1[0] + q[1] * n1[1];
        if d0 <= 0.0 && d1 <= 0.0 {
            // Inside a convex region: the distance to the nearest boundary.
            d0.max(d1)
        } else {
            let [r0, r1] = self.ray_directions();
            let (_, e0) = Self::to_ray(q, r0);
            let (_, e1) = Self::to_ray(q, r1);
            e0.min(e1)
        }
    }

    fn gradient(&self, p: [Scalar; 3]) -> [Scalar; 3] {
        let q = self.local(p);
        let [n0, n1] = self.plane_normals();
        let d0 = q[0] * n0[0] + q[1] * n0[1];
        let d1 = q[0] * n1[0] + q[1] * n1[1];
        if d0 <= 0.0 && d1 <= 0.0 {
            // The active plane's outward normal. A tie is the crease seen from
            // inside; `>=` picks one deterministically, which is the honest
            // answer where the field has no single one.
            let n = if d0 >= d1 { n0 } else { n1 };
            self.unrotate(n)
        } else {
            let [r0, r1] = self.ray_directions();
            let (a0, e0) = Self::to_ray(q, r0);
            let (a1, e1) = Self::to_ray(q, r1);
            let (away, e) = if e0 <= e1 { (a0, e0) } else { (a1, e1) };
            // **The threshold is relative, and a bare `e > 0.0` was a real bug
            // (M-289).** `away = q − dir·t` subtracts two vectors of magnitude
            // `|q|`, so on a point that is on a ray it is not zero but a
            // cancellation residue of order `ε·|q|` — and normalising that
            // returns a **random unit vector**. Every Marching Cubes vertex is
            // on the surface to within an ulp and roughly half of them land
            // epsilon-*outside*, which is how R-006 came to report thousands of
            // "normals pointing into the solid": they were compared against
            // noise. The limit of the exterior gradient approaching a ray along
            // the surface is the plane normal, so the fallback is the right
            // answer and not a fudge.
            let qlen = (q[0] * q[0] + q[1] * q[1]).sqrt();
            if e > 1e-10 * (1.0 + qlen) {
                self.unrotate([away[0] / e, away[1] / e])
            } else {
                let n = if e0 <= e1 { n0 } else { n1 };
                self.unrotate(n)
            }
        }
    }
}
