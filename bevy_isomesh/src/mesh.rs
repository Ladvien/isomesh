//! Writing an extracted surface into a Bevy [`Mesh`].

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
/// use isomesh::mc::MarchingCubes;
/// use isomesh::RuntimeShape3;
///
/// let mut builder = MeshBuilder::new();
/// let mut mc = MarchingCubes::<f32>::new();
/// let shape = RuntimeShape3::new([33; 3]);
/// mc.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut builder);
///
/// let mesh = builder.into_mesh();
/// assert!(mesh.count_vertices() > 0);
/// ```
#[derive(Clone, Debug)]
pub struct MeshBuilder {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    uvs: Vec<[f32; 2]>,
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

    /// Hand the arrays to a Bevy [`Mesh`], by move.
    ///
    /// Consuming rather than borrowing is deliberate: a `Mesh` owns its vertex
    /// data, so the choice is between transferring ownership and copying, and
    /// this is the transfer. Use [`to_bevy_mesh`] when you would rather keep the
    /// buffer and pay for a copy.
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
/// UVs are not emitted here — see [`triplanar_uv`] for why there is no correct
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
    use isomesh::mc::MarchingCubes;

    fn sphere_builder() -> MeshBuilder {
        let mut builder = MeshBuilder::new();
        let mut mc = MarchingCubes::<f32>::new();
        let shape = RuntimeShape3::new([17; 3]);
        mc.extract(
            &Sphere::<f32>::canonical(),
            &shape,
            [-2.0; 3],
            4.0 / 16.0,
            &mut builder,
        );
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
            &isomesh::validate::ValidateConfig::from_cell_size(f64::from(4.0f32 / 16.0)),
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
        let shape = RuntimeShape3::new([17; 3]);
        mc.extract(
            &Sphere::<f32>::canonical(),
            &shape,
            [-2.0; 3],
            4.0 / 16.0,
            &mut buffer,
        );
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
