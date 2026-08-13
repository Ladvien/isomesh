//! Writing an extracted surface into a Bevy [`Mesh`](bevy_mesh::Mesh).

use bevy_asset::RenderAssetUsages;
use bevy_mesh::{Indices, Mesh, PrimitiveTopology};
use isomesh::{MeshBuffer, MeshSink};

/// A [`MeshSink`] whose buffers are exactly the arrays a Bevy [`Mesh`] wants.
///
/// An extractor writes into this directly, so the vertex data is never copied on
/// its way into the asset — [`into_mesh`](Self::into_mesh) hands the `Vec`s over
/// by move. That is the whole point of the type; if you already have an
/// [`isomesh::MeshBuffer`] you are reusing across chunks, use [`to_bevy_mesh`]
/// instead and accept the copy.
///
/// # Example
///
/// ```
/// use bevy_isomesh::MeshBuilder;
/// use isomesh::fields::Sphere;
/// use isomesh::marching_cubes::MarchingCubes;
/// use isomesh::RuntimeShape3;
///
/// let mut builder = MeshBuilder::new();
/// let mut mc = MarchingCubes::<f32>::new();
/// let shape = RuntimeShape3::new([33; 3]).expect("valid shape");
/// mc.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut builder).expect("extraction");
///
/// let mesh = builder.into_mesh();
/// assert!(mesh.count_vertices() > 0);
/// ```
#[derive(Clone, Debug)]
pub struct MeshBuilder {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    uvs: Vec<[f32; 2]>,
    colors: Vec<[f32; 4]>,
    indices: Vec<u32>,
    uv_scale: f32,
}

impl MeshBuilder {
    /// An empty builder.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            positions: Vec::new(),
            normals: Vec::new(),
            uvs: Vec::new(),
            colors: Vec::new(),
            indices: Vec::new(),
            uv_scale: 1.0,
        }
    }

    /// How many world units one texture repeat spans. Default `1.0`.
    #[must_use]
    pub const fn with_uv_scale(mut self, scale: f32) -> Self {
        self.uv_scale = scale;
        self
    }

    /// Truncate without releasing capacity, so one builder can serve many chunks.
    pub fn reset(&mut self) {
        self.positions.clear();
        self.normals.clear();
        self.uvs.clear();
        self.colors.clear();
        self.indices.clear();
    }

    /// Vertices written so far.
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.positions.len()
    }

    /// Triangles written so far.
    #[must_use]
    pub fn triangle_count(&self) -> usize {
        self.indices.len() / 3
    }

    /// The positions written so far.
    ///
    /// Present so a caller can run the core crate's validity harness on exactly
    /// the data it is about to hand to the renderer, rather than on a
    /// reconstruction of it.
    #[must_use]
    pub fn positions(&self) -> &[[f32; 3]] {
        &self.positions
    }

    /// The normals written so far.
    #[must_use]
    pub fn normals(&self) -> &[[f32; 3]] {
        &self.normals
    }

    /// The flat index triples written so far.
    #[must_use]
    pub fn indices(&self) -> &[u32] {
        &self.indices
    }

    /// The per-vertex colour buffer, to be filled or left alone.
    ///
    /// Handed out mutably so
    /// [`isomesh::paint::shade`](isomesh::paint::shade) can write straight into
    /// the array the [`Mesh`] will own, on the same no-copy reasoning as the
    /// rest of this type:
    ///
    /// ```ignore
    /// isomesh::paint::shade(builder.positions(), &world, builder.colors_mut());
    /// ```
    ///
    /// (that borrow needs splitting in real code — read the positions into the
    /// call, or shade into a scratch buffer and swap.)
    ///
    /// A mesh either carries colours or does not, exactly as it either carries
    /// UVs or does not: leave this empty and
    /// [`into_mesh`](Self::into_mesh) omits the attribute entirely. That is one
    /// path writing one attribute, not a fallback — there is no second way to
    /// get a colour onto a vertex here.
    pub fn colors_mut(&mut self) -> &mut Vec<[f32; 4]> {
        &mut self.colors
    }

    /// The per-vertex colours written so far.
    #[must_use]
    pub fn colors(&self) -> &[[f32; 4]] {
        &self.colors
    }

    /// Hand the arrays to a Bevy [`Mesh`], by move.
    ///
    /// Consuming rather than borrowing is deliberate: a `Mesh` owns its vertex
    /// data, so the choice is between transferring ownership and copying, and
    /// this is the transfer. Use [`to_bevy_mesh`] when you would rather keep the
    /// buffer and pay for a copy.
    ///
    /// # Colours are present or absent
    ///
    /// [`Mesh::ATTRIBUTE_COLOR`] is inserted when
    /// [`colors_mut`](Self::colors_mut) was filled and omitted when it was not.
    /// `StandardMaterial` multiplies its `base_color` by the vertex colour, so
    /// padding every mesh with opaque white would be a no-op that still costs
    /// 16 bytes a vertex on the twenty-odd examples that do not paint anything.
    ///
    /// A mismatched length is a caller error and is **not** silently repaired —
    /// Bevy rejects an attribute whose length disagrees with the positions, and
    /// that is the right place for it to fail.
    #[must_use]
    pub fn into_mesh(self) -> Mesh {
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            // The default: the asset stays in main memory as well as being
            // uploaded, which is what a collider baker or a validity check
            // needs. Drop MAIN_WORLD once nothing reads it back.
            RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, self.positions);
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, self.normals);
        mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, self.uvs);
        if !self.colors.is_empty() {
            mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, self.colors);
        }
        mesh.insert_indices(Indices::U32(self.indices));
        mesh
    }
}

impl Default for MeshBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl MeshSink for MeshBuilder {
    type Scalar = f32;

    #[inline]
    fn vertex(&mut self, position: [f32; 3], normal: [f32; 3]) -> u32 {
        let index = self.positions.len();
        assert!(
            index < u32::MAX as usize,
            "MeshBuilder exceeded the u32 index space at {index} vertices"
        );
        self.uvs.push(triplanar_uv(position, normal, self.uv_scale));
        self.positions.push(position);
        self.normals.push(normal);
        index as u32
    }

    #[inline]
    fn triangle(&mut self, a: u32, b: u32, c: u32) {
        self.indices.push(a);
        self.indices.push(b);
        self.indices.push(c);
    }

    #[inline]
    fn reserve(&mut self, vertices: usize, triangles: usize) {
        self.positions.reserve(vertices);
        self.normals.reserve(vertices);
        self.uvs.reserve(vertices);
        self.indices.reserve(triangles * 3);
    }
}

/// Dominant-axis planar projection.
///
/// An isosurface has no natural parameterisation, so there is no correct UV to
/// emit — the good answer is triplanar blending in the shader, which needs no UV
/// attribute at all. This is the cheap stand-in that lets a plain
/// `StandardMaterial` show something sensible: project along whichever axis the
/// normal points most strongly along, and take the other two coordinates.
///
/// It seams visibly where the dominant axis changes. That is inherent to
/// per-vertex planar projection and is the reason shader-side triplanar exists.
#[inline]
fn triplanar_uv(position: [f32; 3], normal: [f32; 3], scale: f32) -> [f32; 2] {
    let n = [normal[0].abs(), normal[1].abs(), normal[2].abs()];
    let (u, v) = if n[0] >= n[1] && n[0] >= n[2] {
        (position[1], position[2])
    } else if n[1] >= n[2] {
        (position[2], position[0])
    } else {
        (position[0], position[1])
    };
    [u / scale, v / scale]
}

/// Copy an [`isomesh::MeshBuffer`] into a Bevy [`Mesh`].
///
/// This **copies**, because the buffer is the thing you are reusing across
/// chunks and a `Mesh` needs to own its data. When the copy matters, extract
/// into a [`MeshBuilder`] instead and let the extractor write straight into the
/// arrays the asset will take.
///
/// UVs are not emitted here — see this module's `triplanar_uv` note for why
/// there is no correct
/// answer, and note that a buffer carries no information this function could
/// use to invent one.
#[must_use]
pub fn to_bevy_mesh(buffer: &MeshBuffer<f32>) -> Mesh {
    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, buffer.positions.clone());
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, buffer.normals.clone());
    mesh.insert_indices(Indices::U32(buffer.indices.clone()));
    mesh
}

#[cfg(test)]
mod tests {
    use super::*;
    use isomesh::RuntimeShape3;
    use isomesh::fields::Sphere;
    use isomesh::marching_cubes::MarchingCubes;

    fn sphere_builder() -> MeshBuilder {
        let mut builder = MeshBuilder::new();
        let mut mc = MarchingCubes::<f32>::new();
        let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
        mc.extract(
            &Sphere::<f32>::canonical(),
            &shape,
            [-2.0; 3],
            4.0 / 16.0,
            &mut builder,
        )
        .expect("extraction");
        builder
    }

    #[test]
    fn extraction_writes_straight_into_the_mesh_arrays() {
        let builder = sphere_builder();
        let vertices = builder.vertex_count();
        let triangles = builder.triangle_count();
        assert!(vertices > 0 && triangles > 0);

        let mesh = builder.into_mesh();
        assert_eq!(mesh.count_vertices(), vertices);
        assert_eq!(
            mesh.indices().map(bevy_mesh::Indices::len),
            Some(triangles * 3)
        );
        assert_eq!(mesh.primitive_topology(), PrimitiveTopology::TriangleList);
        assert!(mesh.attribute(Mesh::ATTRIBUTE_POSITION).is_some());
        assert!(mesh.attribute(Mesh::ATTRIBUTE_NORMAL).is_some());
        assert!(mesh.attribute(Mesh::ATTRIBUTE_UV_0).is_some());
    }

    #[test]
    fn indices_are_u32() {
        let mesh = sphere_builder().into_mesh();
        assert!(matches!(mesh.indices(), Some(Indices::U32(_))));
    }

    /// A mesh nobody painted carries no colour attribute — the twenty-odd
    /// examples that predate E-208 do not pay 16 bytes a vertex for white.
    #[test]
    fn an_unpainted_mesh_has_no_color_attribute() {
        let mesh = sphere_builder().into_mesh();
        assert!(mesh.attribute(Mesh::ATTRIBUTE_COLOR).is_none());
    }

    /// E-208's bridge: `paint::shade` writes into the array the `Mesh` ends up
    /// owning, and the attribute appears because it was filled.
    #[test]
    fn shading_the_builder_puts_colors_on_the_mesh() {
        use isomesh::paint::{Edit, PaintStack, Splat};

        let mut builder = sphere_builder();
        let vertices = builder.vertex_count();

        let log: [Edit<Sphere<f32>, Sphere<f32>, f32>; 1] = [Edit::Spray(Splat {
            shape: Sphere {
                center: [0.0, 1.0, 0.0],
                radius: 0.6,
            },
            color: [1.0, 0.0, 0.0, 1.0],
            softness: 0.1,
            depth: 0.1,
        })];
        let world = PaintStack {
            base: Sphere::<f32>::canonical(),
            edits: &log,
            background: [0.5, 0.5, 0.5, 1.0],
        };

        // The borrow has to be split: shade into a scratch buffer, then swap it
        // into the builder so the Mesh still gets the array by move.
        let mut colors = Vec::new();
        isomesh::paint::shade(builder.positions(), &world, &mut colors);
        core::mem::swap(builder.colors_mut(), &mut colors);
        assert_eq!(builder.colors().len(), vertices);

        let mesh = builder.into_mesh();
        let attribute = mesh
            .attribute(Mesh::ATTRIBUTE_COLOR)
            .expect("colour attribute");
        assert_eq!(attribute.len(), vertices);
    }

    /// The core crate's own validity harness, run on what Bevy is about to
    /// render. If the bridge reordered or dropped anything this would fail.
    #[test]
    fn the_bridged_mesh_is_still_a_closed_surface() {
        let builder = sphere_builder();
        let positions: Vec<[f32; 3]> = builder.positions.clone();
        let indices = builder.indices.clone();
        let report = isomesh::validate::validate_indexed(
            &positions,
            &indices,
            &isomesh::validate::ValidateConfig::from_cell_size(f64::from(4.0f32 / 16.0))
                .expect("valid cell size"),
        );
        assert!(report.is_closed(), "{report}");
        assert_eq!(report.euler_characteristic, 2, "{report}");
    }

    #[test]
    fn reset_keeps_capacity() {
        let mut builder = sphere_builder();
        let capacity = builder.positions.capacity();
        assert!(capacity > 0);
        builder.reset();
        assert_eq!(builder.vertex_count(), 0);
        assert_eq!(builder.positions.capacity(), capacity);
    }

    #[test]
    fn buffer_conversion_matches_the_builder() {
        let mut buffer = MeshBuffer::<f32>::new();
        let mut mc = MarchingCubes::<f32>::new();
        let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
        mc.extract(
            &Sphere::<f32>::canonical(),
            &shape,
            [-2.0; 3],
            4.0 / 16.0,
            &mut buffer,
        )
        .expect("extraction");
        let mesh = to_bevy_mesh(&buffer);
        assert_eq!(mesh.count_vertices(), buffer.vertex_count());
        assert_eq!(
            mesh.indices().map(bevy_mesh::Indices::len),
            Some(buffer.indices.len())
        );
    }

    #[test]
    fn uvs_follow_the_dominant_axis() {
        // A +x-facing normal projects the y and z coordinates.
        assert_eq!(
            triplanar_uv([1.0, 2.0, 3.0], [1.0, 0.0, 0.0], 1.0),
            [2.0, 3.0]
        );
        // A +y-facing normal projects z and x.
        assert_eq!(
            triplanar_uv([1.0, 2.0, 3.0], [0.0, 1.0, 0.0], 1.0),
            [3.0, 1.0]
        );
        // A +z-facing normal projects x and y.
        assert_eq!(
            triplanar_uv([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 1.0),
            [1.0, 2.0]
        );
        // Scale divides.
        assert_eq!(
            triplanar_uv([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 2.0),
            [0.5, 1.0]
        );
    }
}
