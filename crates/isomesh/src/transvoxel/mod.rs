//! Transvoxel — stitching two blocks meshed at different resolutions.
//!
//! Field-derived LOD (G-004) meshes a distant block at half resolution, and
//! where a half-resolution block meets a full-resolution one the two surfaces do
//! not meet: the coarse side places one vertex on an edge where the fine side
//! places two, and the gap between them is a crack you can see the sky through.
//! Lengyel 2010 §4.3 is blunt about the alternatives — *"existing stitching
//! methods can be carried over to voxel-based terrain only under very limited
//! conditions, and none of them are capable of handling all of the cases that can
//! arise in an unrestricted voxel map."*
//!
//! Transvoxel closes it with a **transition cell** on the boundary, and the
//! dissertation's key move is to split that cell in two so the interesting half
//! depends on nine sample values rather than thirteen, twenty, or a million. See
//! [`table`] for the geometry, the derivation, and the number that checks it.
//!
//! # Status
//!
//! A-011a — the case table — is here, and so is [`cell`], which places a
//! transition cell's crossings in world space and asserts the identity the seam
//! rests on: a half-resolution crossing lands **bit-identically** on a vertex the
//! coarse neighbour's Marching Cubes pass produced.
//!
//! A-011b's triangulation is here too — [`cell::TransitionCell::emit`] winds
//! the cycles [`table::transition_links`] gives into a centroid fan, and
//! `transition_cells_close_the_gap_between_two_resolutions` asserts zero
//! seam-plane boundary edges across a real full-resolution/half-resolution
//! chunk pair. A-011c's inset taper lives in [`inset`], with its own
//! real-width acceptance test.

pub mod cell;
pub mod inset;
pub mod table;
