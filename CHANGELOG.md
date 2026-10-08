# Changelog

## 0.2.0

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
