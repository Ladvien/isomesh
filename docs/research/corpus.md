# Corpus: the home-still reading list, by domain

Created by loop iteration 1 (2026-10-30). Each stem is a home-still `doc_id`, found with `distill_search` on that date; open one with `markdown_read` or `doc_link`. A domain's list is what the corpus returned for it, filtered by hand to the documents actually on topic, so it is a starting point rather than a census. Each pass that researches adds what it found under the domain it belongs to; a document the corpus lacks goes under the research note's `## Pending index` instead.

The longer, older reading notes are the dated files beside this one (`2026-08-10-meshing-algorithm-catalog-v2.md`, `2026-08-15-sdf-corpus-buildout.md`, `2026-08-18-corpus-audit-and-procurement.md`); this file is the per-domain index the loop keeps current.

## Marching Cubes and its ambiguities

- `10.1145_37401.37422`: Lorensen and Cline 1987, *Marching Cubes*.
- `10.1109_visual.1991.175782`: Nielson and Hamann 1991, the asymptotic decider.
- `10.1109_tvcg.2003.1207437`: *On Marching Cubes* (Nielson 2003).
- `10.1016_j.cag.2013.04.004`: *Practical considerations on Marching Cubes 33 topological correctness* (2013).
- `10.1016_j.cag.2006.07.021`: Newman and Yi 2006, a survey of the Marching Cubes family.
- `10.1007_s00371-021-02139-w`: *A parallel dual marching cubes approach to quad only surface reconstruction* (2021).

## Dual Contouring and the QEF

- `SchaeferWarren2`: Schaefer and Warren, *Dual Contouring: "The Secret Sauce"*.
- `10.48550_arXiv.2604.00157`: *Dual Contouring of Signed Distance Data* (2026).
- `10.1111_j.1467-8659.2005.00879.x`: *Cubical Marching Squares* (2005).
- `10.48550_arxiv.2308.05371`: *FlexiCubes* (2023).
- `10.48550_arxiv.2202.01999`: *Neural Dual Contouring* (2022); added from `catalog_recent`, iteration 2.

## Surface Nets

- `10.1007_bfb0056277`: Gibson 1998, *Constrained elastic surface nets*.

## Transvoxel, LOD and seams

- `transvoxel_dissertation_lengyel2010`: Lengyel 2010, *Voxel-Based Terrain for Real-Time Virtual Simulations*; deployed in the C4 Engine (p. 99).
- `transvoxel_paper_lengyel2010`: *The Transvoxel Algorithm* (case tables).
- `10.1145_1015706.1015799_hhoppe`: Losasso and Hoppe 2004, *Geometry Clipmaps*.
- `10.1109_visual.1998.745282`: Hoppe 1998, *Smooth view-dependent level-of-detail control and its application to terrain rendering*.

## Signed distance fields and CSG

- `10.14733_cadaps.2023.1154-1174`: *Operations on Signed Distance Function Estimates* (2023).
- `10.1145_3730845`: *Sphere Carving: Bounding Volumes for Signed Distance Fields* (2025).
- `10.1016_j.cag.2010.07.003`: *Fast reliable interrogation of procedurally defined implicit surfaces using extended revised affine arithmetic* (2010).
- `10.1016_j.gmod.2021.101098`: *Hybrid Function Representation for Heterogeneous Objects* (2021).
- `10.48550_arxiv.2205.14151`: *Interactive and Robust Mesh Booleans* (2022); mesh-side CSG with exact predicates, the alternative to field-side CSG. Added iteration 2.

## Voxel, cave and destructible terrain in games

- `transvoxel_dissertation_lengyel2010`: C4 Engine voxel terrain; p. 28 says a block is a renderable node only where the surface crosses it, Appendix A has the editing tools.
- `10.1145_3102163.3102216`: *Procedural feature generation for volumetric terrains* (2017).
- `10.48550_arXiv.2505.02017`: *Aokana: A GPU-Driven Voxel Rendering Framework for Open World Games* (2025); describes Minecraft's chunking.
- `s2025-advances-materials`: SIGGRAPH 2025 *Advances in Real-Time Rendering in Games* course.
- `10.48550_arxiv.2306.09310`: *Infinite Photorealistic Worlds using Procedural Generation* (Infinigen, 2023); SDF terrain with caves, meshed per view. Added iteration 2.
- `10.48550_arxiv.2210.11463`: *Breaking Bad: A Dataset for Geometric Fracture and Reassembly* (2022); fracture of meshes, the destructible-object side. Added iteration 2.
- `10.48550_arxiv.2205.02961`: *Approximate Convex Decomposition for 3D Meshes with Collision-Aware Concavity and Tree Search* (2022); collider meshes from a surface. Added iteration 2.

## Mesh validity and repair

- `scanpaper_final`: Ju 2004, *Robust Repair of Polygonal Models*: a polygon soup in, an intermediate volume grid, the output surface contoured from it.
- `10.1145_2431211.2431214`: *Polygon Mesh Repairing: An Application Perspective* (2013).
- `10.1007_s00371-010-0416-3`: *A lightweight approach to repairing digitized polygon meshes* (2010).
- `10.1109_tvcg.2003.1196006`: *Simplification and repair of polygonal models using volumetric techniques* (2003).
- `sig2026_Subgrid_Marching_Tetrahedra`: *Subgrid Marching Tetrahedra* (2026).

## GPU meshing and the render path

- `10.1007_s13173-012-0097-z`: *Marching cubes technique for volumetric visualization accelerated with graphics processing units* (2012).
- `book_real_time_rendering_4th_akenine`: *Real-Time Rendering*, 4th ed.; pp. 816-818 on draw-call and state-change overhead.
- `WebGPU_Sourcebook`: *The WebGPU Sourcebook* (2024).
