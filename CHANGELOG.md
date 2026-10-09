# Changelog

## 0.3.0

Blends between curved surfaces, exact B-spline surfaces and spindle
tori. Version 0.2.0 was tagged in this repository but never published
to crates.io; its changes (blend tools for straight and circular edges,
`Options::should_stop`) are first published here.

- `Surface::Torus` accepts spindle and horn tori (major radius above
  zero, no longer above the minor radius), written to STEP as
  `DEGENERATE_TOROIDAL_SURFACE(..., .T.)` when the major radius is the
  smaller. Blend tools build convex rim fillets up to the rim's radius.
  A caller that relied on `reconstruct` refusing such a torus as
  malformed sees it accepted now.
- `blend::tools`: the edges of a sphere corner take the least margin of
  every tangent chain one of them belongs to (before, a corner edge
  chained to an arc could make an open tool).
- B-spline surfaces: `Surface::BSpline` and `BSplineSurface` (clamped,
  optionally rational), reconstructed, validated, measured and written
  to STEP, and the `spline` module (evaluation and projection, curve and
  surface fitting, canal and ruled surfaces). Breaking for a caller
  that matches `Surface` exhaustively.
- `Contact::Boundary`: a B-spline patch touching another surface along
  a side of its domain, in `Report::tangencies` and `find_tangencies`
  (one entry per side). Breaking for a caller that matches `Contact`
  exhaustively.
- `Tolerances::surface_fit` (default 1e-7): how far a B-spline patch's
  side may be from the surface it touches. Breaking for a caller that
  builds `Tolerances` with a struct literal (use `..Default::default()`).
- Blends between curved surfaces with no common axis: `blend::Path::Curve`
  (an edge's points, and per face the triangles near it that the tool is
  conformed to). The ball's centre is marched along the intersection of
  the two faces' offsets, and the blend is a canal (fillet) or ruled
  (chamfer) B-spline surface fitted to it, two patches for a closed curve;
  its ends are `End::Plane` or `End::Open`. Breaking for a caller that
  matches `Path` exhaustively.
- `blend::Tool::fit`: how far a tool's B-spline blends may be from the
  faces they meet (0 when every blend is a quadric or torus); a caller
  reconstructing the result sets `Tolerances::surface_fit` to at least
  that. Breaking for a caller that builds `Tool` with a struct literal.
- `blend::curve_sections` (a curve's blend across it at fractions of its
  length) and `blend::check` (what `tools` would refuse, without fitting
  or meshing the curves' blends), for size checks run many times.
- Reconstruction: a parameter-space curve along a tangent contact is
  refined until its image stays on the edge's own curve too, not only
  on the other face (across a tangent contact the distance to that face
  changes only to second order, so a coarse curve passed 1e-4 off it).
  Files with tangent blends grow (a plate's straight and rim fillets:
  137 KB to 325 KB). A mesh chain along a B-spline patch's contact side
  is accepted within a fifth of the patch's width (before: 1e-4 of the
  chain's length), measured after moving it onto the other face; and a
  corner where contact sides of different patches meet is taken where
  they meet.
- `spline::Evaluator::project_extended`: projection that may run a
  little past the patch's sides onto the surface its end spans continue
  (what reconstruction's implicit form uses), for measuring how far
  points beside a side stand off the patch.
- Reconstruction of B-spline faces is several times faster: the corner
  and edge solver projects once per step for both the value and the
  gradient, stops (on a patch only) once the residual stops halving and
  keeps the best point, and an edge fit on a patch that does not
  converge is given up after three doublings without progress. Points
  that solves left wandering at the projection's rounding are now the
  best ones found, so fitted edges change in their last digits.
- A sphere face none of whose bounding circles can be parallels (a box
  corner's patch) is framed with its poles square to the face and as
  many bounding circles as possible meridians, whose parameter-space
  curves are segments: a rounded box's STEP file shrinks (a plate's
  straight and rim fillets: 325 KB to 198 KB). A periodic face's seam
  avoids passing within 1e-6 of the model's size from a vertex of
  another loop, which split off a degenerate edge.

## 0.2.0 (not published to crates.io)

- **Blend tools** (`blend`): solids that fillet or chamfer a B-rep's
  edges when a mesh kernel subtracts them (convex edges) or adds them
  (concave ones), with every triangle tagged with its exact surface, so
  the boolean's result reconstructs with true blend faces. Straight edges
  between planes at any angle and parallel plane–cylinder and
  cylinder–cylinder pairs give cylinder fillets and plane chamfers; ends
  are cut, mitred, run on, or meet a sphere patch where three fillets
  meet. Circular edges between surfaces of revolution about one axis give
  torus fillets and cone chamfers, sectioned at the caller's polygon so
  they meet a faceted face without creases. Edges that run on into each
  other come out as one solid. `blend::section` gives an edge's
  cross-section for checks before anything is built.
- **Cancellation**: `Options::should_stop` is polled between the stages
  of reconstruction and once per face; a stop returns `Error::Stopped`.

**Breaking:** `Options` is no longer `Copy` (it is `Clone`) and its
`PartialEq` compares the stop signal by identity; `Error` has a new
variant, `Stopped`.


## 0.1.0

First release: exact B-rep reconstruction from a triangle mesh tagged
with planes, cylinders, cones, spheres and tori, with exact line, circle
and ellipse edges and B-spline edges elsewhere; tangency, touching-point
and sliver handling; a faceted fallback for untagged triangles; a
structural and geometric validator; exact volume and area; and
deterministic STEP AP214 output, the same bytes on every platform,
including wasm32.
