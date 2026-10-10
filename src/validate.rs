//! A structural and geometric validator for [`Brep`], so that tests (and
//! callers) can tell a valid solid from a broken one without a CAD kernel.
//!
//! Validity alone is not enough: testing found a file that read back as
//! a valid solid with the wrong volume. Callers should also compare
//! [`crate::measure`] with an independent volume (the mesh's).

use crate::bspline;
use crate::curve::{self, CurveEval};
use crate::math::*;
use crate::measure::{face_param, shell_volumes};
use crate::model::Brep;

/// What [`validate`] found.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Validation {
    /// Every problem found, in words. Empty for a valid B-rep.
    pub errors: Vec<String>,
    /// The genus from the Euler–Poincaré formula
    /// `V − E + 2F − L = 2 (S − G)`, if that gives a whole number.
    pub genus: Option<i64>,
    /// The signed volume of each shell (negative for voids), when it could
    /// be integrated.
    pub shell_volumes: Vec<f64>,
    /// Findings that are not errors, in words: the B-rep has fewer
    /// handles than the input mesh.
    pub notes: Vec<String>,
    /// The faces the errors are on (an edge's error is on the faces that
    /// use it), sorted, each once. With `Report::face_triangles` they say
    /// which input triangles to build some other way. Errors of the whole
    /// B-rep (its genus) name no face.
    pub error_faces: Vec<u32>,
}

impl Validation {
    /// Whether no problem was found.
    pub fn is_valid(&self) -> bool {
        self.errors.is_empty()
    }
}

/// Checks `brep`:
///
/// - every index is in range;
/// - every loop is closed (each coedge ends where the next starts);
/// - every edge is used by exactly two coedges, in opposite directions
///   (orientation consistency), by two faces, or twice by one face if it
///   is a seam;
/// - every face is in exactly one shell, and shells are edge-connected;
/// - the Euler–Poincaré formula gives a whole, non-negative genus, no more
///   than the input mesh's when the report has it (fewer is a note);
/// - no face's boundary crosses itself, or touches itself with a corner
///   inside one of its edges (within `tolerance`);
/// - geometry, within `tolerance` (model units): curve ends at their
///   vertices, edges on both their faces' surfaces, parameter-space curves
///   present on curved faces with ends at the edge's vertices;
/// - every shell encloses a positive volume, or a negative one if it is a
///   void.
pub fn validate(brep: &Brep, tolerance: f64) -> Validation {
    let mut out = validate_unsorted(brep, tolerance);
    out.error_faces.sort_unstable();
    out.error_faces.dedup();
    out
}

fn validate_unsorted(brep: &Brep, tolerance: f64) -> Validation {
    let mut out = Validation::default();
    let (nv, ne, nf) = (brep.vertices.len(), brep.edges.len(), brep.faces.len());
    // The faces using each edge, for saying where an edge's error is.
    let mut edge_faces: Vec<Vec<u32>> = vec![Vec::new(); ne];
    for (fi, f) in brep.faces.iter().enumerate() {
        for c in f.loops.iter().flat_map(|l| &l.coedges) {
            if let Some(v) = edge_faces.get_mut(c.edge as usize)
                && !v.contains(&(fi as u32))
            {
                v.push(fi as u32);
            }
        }
    }
    let at_face = |f: usize| vec![f as u32];
    let at_edge = |e: usize| edge_faces.get(e).cloned().unwrap_or_default();
    let err = |out: &mut Validation, s: String, faces: Vec<u32>| {
        if out.errors.len() < 100 {
            out.errors.push(s);
        }
        out.error_faces.extend(faces);
    };
    for (i, e) in brep.edges.iter().enumerate() {
        if e.start as usize >= nv || e.end as usize >= nv {
            err(
                &mut out,
                format!("edge {i}: vertex out of range"),
                at_edge(i),
            );
        }
        if e.range[1] <= e.range[0] {
            err(
                &mut out,
                format!("edge {i}: empty parameter range"),
                at_edge(i),
            );
        }
    }
    let mut uses: Vec<Vec<(usize, bool)>> = vec![Vec::new(); ne];
    let mut nloops = 0i64;
    for (fi, f) in brep.faces.iter().enumerate() {
        if f.loops.is_empty() {
            err(&mut out, format!("face {fi}: no loops"), at_face(fi));
        }
        let outer = f.loops.iter().filter(|l| l.outer).count();
        if outer > 1 {
            err(
                &mut out,
                format!("face {fi}: more than one outer loop"),
                at_face(fi),
            );
        } else if outer == 0
            && !f.loops.is_empty()
            && matches!(
                f.surface,
                crate::model::Surface::Plane { .. } | crate::model::Surface::Faceted
            )
        {
            // A bounded planar region has an outer loop. Without one the
            // face's loops all wind the wrong way for its normal: a face
            // folded through a tunnel of no thickness in the mesh (BOSL2
            // `attachments__084`, which OCCT read back as unorientable).
            err(
                &mut out,
                format!("face {fi}: no outer loop (its loops wind against its normal)"),
                at_face(fi),
            );
        }
        for (li, lp) in f.loops.iter().enumerate() {
            nloops += 1;
            if lp.coedges.is_empty() {
                err(&mut out, format!("face {fi} loop {li}: empty"), at_face(fi));
                continue;
            }
            for c in &lp.coedges {
                if c.edge as usize >= ne {
                    err(
                        &mut out,
                        format!("face {fi} loop {li}: edge out of range"),
                        at_face(fi),
                    );
                }
            }
            if out.errors.iter().any(|e| e.contains("out of range")) {
                continue;
            }
            let ends = |c: &crate::model::Coedge| {
                let e = &brep.edges[c.edge as usize];
                if c.forward {
                    (e.start, e.end)
                } else {
                    (e.end, e.start)
                }
            };
            // A closed edge (a whole circle or ellipse) is a loop by
            // itself. Inside a longer loop it is a hole that touches the
            // boundary at a corner: OCCT splits it off as a wire of its
            // own and rejects the face for wires that nest wrongly (BOSL2
            // `threading__048` written partly as facets, where a facet's
            // plane cut a cone's tip off in an ellipse through one of the
            // facet's corners). On a periodic face a closed edge shares a
            // loop with the seam (a cylinder's rims), which is right.
            let planar = matches!(
                f.surface,
                crate::model::Surface::Plane { .. } | crate::model::Surface::Faceted
            );
            if planar && lp.coedges.len() > 1 {
                for c in &lp.coedges {
                    let e = &brep.edges[c.edge as usize];
                    if e.start == e.end && !e.seam {
                        err(
                            &mut out,
                            format!(
                                "face {fi} loop {li}: closed edge {} inside a longer loop (a hole touching the boundary)",
                                c.edge
                            ),
                            at_edge(c.edge as usize),
                        );
                    }
                }
            }
            // Two straight edges make a loop with no area: they overlap
            // rather than cross, so the crossing check below cannot see it.
            if lp.coedges.len() == 2
                && lp.coedges.iter().all(|c| {
                    matches!(
                        brep.edges[c.edge as usize].curve,
                        crate::model::Curve::Line { .. }
                    )
                })
            {
                err(
                    &mut out,
                    format!("face {fi} loop {li}: two straight edges along one line (no area)"),
                    at_face(fi),
                );
            } else if lp.coedges.len() == 2 {
                // Two curves that run along each other within the
                // tolerance make a sliver of no area too, on any surface
                // (a cone between two facets whose planes cut it almost
                // alike: BOSL2 `threading__048` written partly as facets,
                // which OCCT read back as an unorientable face).
                let (a, b) = (
                    &brep.edges[lp.coedges[0].edge as usize],
                    &brep.edges[lp.coedges[1].edge as usize],
                );
                let same_ends =
                    (a.start, a.end) == (b.start, b.end) || (a.start, a.end) == (b.end, b.start);
                // Not a seam, which a lone sphere's or torus's face runs
                // along both ways.
                let distinct = lp.coedges[0].edge != lp.coedges[1].edge && !a.seam && !b.seam;
                if distinct && same_ends && a.start != a.end && run_along(a, b, tolerance) {
                    err(
                        &mut out,
                        format!("face {fi} loop {li}: two edges along one curve (no area)"),
                        at_face(fi),
                    );
                }
            }
            for k in 0..lp.coedges.len() {
                let a = ends(&lp.coedges[k]).1;
                let b = ends(&lp.coedges[(k + 1) % lp.coedges.len()]).0;
                if a != b {
                    err(
                        &mut out,
                        format!("face {fi} loop {li}: open after coedge {k}"),
                        at_face(fi),
                    );
                }
                let c = &lp.coedges[k];
                uses[c.edge as usize].push((fi, c.forward));
            }
        }
    }
    if out.errors.iter().any(|e| e.contains("out of range")) {
        return out;
    }
    for (i, u) in uses.iter().enumerate() {
        let seam = brep.edges[i].seam;
        if u.len() != 2 {
            err(
                &mut out,
                format!("edge {i}: used {} times, not twice", u.len()),
                at_edge(i),
            );
        } else if u[0].1 == u[1].1 {
            err(
                &mut out,
                format!("edge {i}: both uses run the same way"),
                at_edge(i),
            );
        } else if (u[0].0 == u[1].0) != seam {
            err(
                &mut out,
                format!(
                    "edge {i}: {}",
                    if seam {
                        "a seam used by two faces"
                    } else {
                        "used twice by one face but not a seam"
                    }
                ),
                at_edge(i),
            );
        }
    }
    // Shells.
    let mut shell_of = vec![usize::MAX; nf];
    for (si, sh) in brep.shells.iter().enumerate() {
        for &f in &sh.faces {
            if f as usize >= nf {
                err(
                    &mut out,
                    format!("shell {si}: face out of range"),
                    Vec::new(),
                );
            } else if shell_of[f as usize] != usize::MAX {
                err(&mut out, format!("face {f}: in two shells"), vec![f]);
            } else {
                shell_of[f as usize] = si;
            }
        }
    }
    if let Some(f) = shell_of.iter().position(|&s| s == usize::MAX) {
        err(&mut out, format!("face {f}: in no shell"), at_face(f));
    }
    for (i, u) in uses.iter().enumerate() {
        if u.len() == 2 && shell_of[u[0].0] != shell_of[u[1].0] {
            err(&mut out, format!("edge {i}: joins two shells"), at_edge(i));
        }
    }
    let mut uf = UnionFind::new(nf);
    for u in &uses {
        if u.len() == 2 {
            uf.join(u[0].0, u[1].0);
        }
    }
    for (si, sh) in brep.shells.iter().enumerate() {
        if let Some(&f0) = sh.faces.first() {
            let r = uf.find(f0 as usize);
            if sh.faces.iter().any(|&f| uf.find(f as usize) != r) {
                err(
                    &mut out,
                    format!("shell {si}: not connected"),
                    sh.faces.clone(),
                );
            }
        }
    }
    // Euler–Poincaré.
    let mut used_v = vec![false; nv];
    for e in &brep.edges {
        used_v[e.start as usize] = true;
        used_v[e.end as usize] = true;
    }
    let v_count = used_v.iter().filter(|&&u| u).count() as i64;
    let chi = v_count - ne as i64 + 2 * nf as i64 - nloops;
    if chi % 2 != 0 {
        err(
            &mut out,
            format!("Euler–Poincaré: V − E + 2F − L = {chi} is odd"),
            Vec::new(),
        );
    } else {
        let g = brep.shells.len() as i64 - chi / 2;
        out.genus = Some(g);
        if g < 0 {
            err(
                &mut out,
                format!("Euler–Poincaré: genus {g} is negative"),
                Vec::new(),
            );
        }
        // More handles than the mesh means the B-rep joined what the
        // mesh keeps apart (corners merged that should not be): an error.
        // Fewer can be right. Rounding can leave a mesh with a tunnel of
        // no thickness, where a plane meets one that is equal within
        // the merge tolerance but whose triangles straddle it (a cube
        // standing on a rotated prism, BOSL2 `attachments__079`): the
        // two sides of the tunnel are one exact face, so the B-rep has
        // no handle there and is the exact model's. A tunnel with any
        // volume would change the volume, which callers must compare
        // with the mesh's anyway (see the module documentation).
        let r = &brep.report;
        if r.mesh_components == brep.shells.len() && r.mesh_components > 0 {
            if g > r.mesh_genus {
                err(
                    &mut out,
                    format!("genus {g} is more than the input mesh's {}", r.mesh_genus),
                    Vec::new(),
                );
            } else if g < r.mesh_genus {
                out.notes.push(format!(
                    "genus {g} is less than the input mesh's {} (a tunnel of no thickness in the mesh)",
                    r.mesh_genus
                ));
            }
        }
    }
    // Geometry.
    for (i, e) in brep.edges.iter().enumerate() {
        if e.start != e.end {
            let mut prev: Option<V> = None;
            let len: f64 = curve::samples(&e.curve, e.range, 8)
                .filter_map(|q| {
                    let step = prev.map(|p| (q - p).len());
                    prev = Some(q);
                    step
                })
                .sum();
            if len < tolerance {
                err(
                    &mut out,
                    format!("edge {i}: degenerate (length {len:.2e})"),
                    at_edge(i),
                );
            }
        }
        let (a, b) = (
            curve::eval(&e.curve, e.range[0]),
            curve::eval(&e.curve, e.range[1]),
        );
        let (pa, pb) = (
            V::from(brep.vertices[e.start as usize]),
            V::from(brep.vertices[e.end as usize]),
        );
        let d = (a - pa).len().max((b - pb).len());
        if d > tolerance {
            err(
                &mut out,
                format!("edge {i}: curve ends {d:.2e} from its vertices"),
                at_edge(i),
            );
        }
    }
    for (fi, f) in brep.faces.iter().enumerate() {
        let Some(p) = face_param(f) else {
            err(
                &mut out,
                format!("face {fi}: unsupported surface {}", f.surface.kind()),
                at_face(fi),
            );
            continue;
        };
        for lp in &f.loops {
            for c in &lp.coedges {
                let e = &brep.edges[c.edge as usize];
                let n = curve::sample_count(&e.curve, e.range).min(256);
                let worst = curve::samples(&e.curve, e.range, n)
                    .map(|q| p.s.f(q).abs())
                    .fold(0.0, f64::max);
                if worst > tolerance {
                    err(
                        &mut out,
                        format!("face {fi}: edge {} is {worst:.2e} off the surface", c.edge),
                        at_edge(c.edge as usize),
                    );
                }
                if !p.has_pcurves() {
                    continue;
                }
                let Some(pc) = &c.pcurve else {
                    err(
                        &mut out,
                        format!("face {fi}: edge {} has no parameter-space curve", c.edge),
                        at_face(fi),
                    );
                    continue;
                };
                let [k0, k1] = [pc.knots[0], *pc.knots.last().unwrap_or(&0.0)];
                let (q0, q1) = (bspline::eval(pc, k0), bspline::eval(pc, k1));
                let (s0, s1) = (p.eval(q0[0], q0[1]), p.eval(q1[0], q1[1]));
                let ev = CurveEval::new(&e.curve);
                let (c0, c1) = (ev.at(e.range[0]), ev.at(e.range[1]));
                let d = (s0 - c0).len().max((s1 - c1).len());
                if d > tolerance {
                    err(
                        &mut out,
                        format!(
                            "face {fi}: edge {}'s parameter-space curve ends {d:.2e} from the edge",
                            c.edge
                        ),
                        at_face(fi),
                    );
                }
            }
        }
    }
    // A loop of a planar face narrower than the tolerance (a needle
    // triangle of the mesh written as a facet, or a slit whose two sides
    // are not both lines) has no trustworthy area or orientation: its
    // sides overlap within the tolerance. OCCT reads such a wire as
    // crossing itself, or a hole as badly oriented (BOSL2 `distributors`
    // and `example017.scad`, written partly as facets, read back invalid
    // with them). Its width is taken as 4 · area / perimeter, a
    // triangle's height on its long side, a slit's half-width.
    for (fi, f) in brep.faces.iter().enumerate() {
        let planar = matches!(
            f.surface,
            crate::model::Surface::Plane { .. } | crate::model::Surface::Faceted
        );
        if !planar {
            continue;
        }
        let Some(lines) = face_polylines(brep, fi, 1) else {
            continue;
        };
        for (li, pts) in lines.iter().enumerate() {
            if pts.len() < 2 {
                continue;
            }
            let (mut area, mut perimeter) = (0.0, 0.0);
            for k in 0..pts.len() {
                let (p, q) = (pts[k], pts[(k + 1) % pts.len()]);
                area += p[0] * q[1] - p[1] * q[0];
                perimeter += ((q[0] - p[0]).powi(2) + (q[1] - p[1]).powi(2)).sqrt();
            }
            let width = 2.0 * area.abs() / perimeter.max(1e-300);
            if width < tolerance {
                err(
                    &mut out,
                    format!(
                        "face {fi} loop {li}: {width:.2e} wide, narrower than the tolerance (no area)"
                    ),
                    at_face(fi),
                );
            }
        }
        // Each hole inside the outer loop and outside the other holes.
        // OCCT rejects a face whose wires nest otherwise ("invalid
        // imbrication of wires"), and no other check here sees it when no
        // two loops cross. None of the corpora's files has had one, but a
        // face rebuilt partly from facets could. A hole may touch the outer
        // loop or another hole, so a corner within the tolerance of the
        // other loop is not counted. Holes are compared only where their
        // boxes meet: a plate with 400 holes would otherwise compare every
        // pair point by point.
        let Some(outer) = f.loops.iter().position(|l| l.outer) else {
            continue;
        };
        let boxes: Vec<[f64; 4]> = lines
            .iter()
            .map(|q| {
                q.iter().fold(
                    [
                        f64::INFINITY,
                        f64::INFINITY,
                        f64::NEG_INFINITY,
                        f64::NEG_INFINITY,
                    ],
                    |b, p| {
                        [
                            b[0].min(p[0]),
                            b[1].min(p[1]),
                            b[2].max(p[0]),
                            b[3].max(p[1]),
                        ]
                    },
                )
            })
            .collect();
        let meet = |i: usize, j: usize| {
            let (a, b) = (boxes[i], boxes[j]);
            a[0] <= b[2] + tolerance
                && b[0] <= a[2] + tolerance
                && a[1] <= b[3] + tolerance
                && b[1] <= a[3] + tolerance
        };
        for (li, pts) in lines.iter().enumerate() {
            if li == outer {
                continue;
            }
            let strictly = |q: &[[f64; 2]], inside: bool| {
                pts.iter()
                    .any(|&p| polyline_distance(q, p) > tolerance && winds(q, p) == inside)
            };
            let other_holes = lines
                .iter()
                .enumerate()
                .any(|(lj, q)| lj != li && lj != outer && meet(li, lj) && strictly(q, true));
            if strictly(&lines[outer], false) || other_holes {
                err(
                    &mut out,
                    format!(
                        "face {fi} loop {li}: a hole not inside its face's outer loop, or inside another hole"
                    ),
                    at_face(fi),
                );
            }
        }
    }
    if out.errors.is_empty() {
        // Reconstruction made the crossing check already; the touch check
        // is made again, since `tolerance` is the caller's and may be
        // larger than reconstruction's.
        let crossed = !brep.report.faces_uncrossed;
        for fi in 0..nf {
            let crossing = if crossed {
                face_crossing(brep, fi)
            } else {
                None
            };
            if let Some(e) = crossing.or_else(|| face_touch(brep, fi, tolerance)) {
                err(&mut out, e, at_face(fi));
            }
        }
    }
    if out.errors.is_empty() {
        match shell_volumes(brep) {
            Ok(vols) => {
                for (si, (&v, sh)) in vols.iter().zip(&brep.shells).enumerate() {
                    if (v > 0.0) == sh.void || v == 0.0 {
                        err(
                            &mut out,
                            format!(
                                "shell {si}: encloses volume {v:.6}, {}",
                                if sh.void {
                                    "but is a void"
                                } else {
                                    "inside out"
                                }
                            ),
                            sh.faces.clone(),
                        );
                    }
                }
                out.shell_volumes = vols;
            }
            Err(e) => err(&mut out, format!("volume: {e}"), Vec::new()),
        }
    }
    out
}

/// Whether edge `b` stays within `tol` of edge `a` everywhere (sampled
/// densely, so the chords' sagitta is far below `tol` on the short
/// curves where this matters).
fn run_along(a: &crate::model::Edge, b: &crate::model::Edge, tol: f64) -> bool {
    let pa = curve::sample(&a.curve, a.range, 512);
    let pb = curve::sample(&b.curve, b.range, 64);
    pb.iter().all(|&p| {
        pa.windows(2)
            .map(|w| {
                let d = w[1] - w[0];
                let l2 = d.dot(d);
                let t = if l2 > 0.0 {
                    ((p - w[0]).dot(d) / l2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (p - (w[0] + d * t)).len()
            })
            .fold(f64::INFINITY, f64::min)
            < tol
    })
}

/// Whether the closed polyline `q` winds around `p` (crossing parity).
fn winds(q: &[[f64; 2]], p: [f64; 2]) -> bool {
    let mut inside = false;
    for k in 0..q.len() {
        let (a, b) = (q[k], q[(k + 1) % q.len()]);
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
            if x > p[0] {
                inside = !inside;
            }
        }
    }
    inside
}

/// The distance from `p` to the closed polyline `q`.
fn polyline_distance(q: &[[f64; 2]], p: [f64; 2]) -> f64 {
    let mut best = f64::INFINITY;
    for k in 0..q.len() {
        let (a, b) = (q[k], q[(k + 1) % q.len()]);
        let d = [b[0] - a[0], b[1] - a[1]];
        let l2 = d[0] * d[0] + d[1] * d[1];
        let t = if l2 > 0.0 {
            (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / l2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let (x, y) = (a[0] + d[0] * t - p[0], a[1] + d[1] * t - p[1]);
        best = best.min((x * x + y * y).sqrt());
    }
    best
}

/// A face's loops as polylines in its parameter plane (plane coordinates
/// for planes, the parameter-space curves otherwise), `k` times the usual
/// sampling. `None` if a curved face lacks parameter-space curves.
fn face_polylines(b: &Brep, fi: usize, k: usize) -> Option<Vec<Vec<[f64; 2]>>> {
    let f = &b.faces[fi];
    let p = face_param(f)?;
    let mut out = Vec::with_capacity(f.loops.len());
    for lp in &f.loops {
        let mut pts: Vec<[f64; 2]> = Vec::new();
        for c in &lp.coedges {
            let e = &b.edges[c.edge as usize];
            let mut s: Vec<[f64; 2]> = if p.has_pcurves() {
                let pc = c.pcurve.as_ref()?;
                let (t0, t1) = (pc.knots[0], *pc.knots.last()?);
                let n = k * (4 * bspline::spans(pc).len()).clamp(2, 512);
                (0..=n)
                    .map(|i| bspline::eval(pc, t0 + (t1 - t0) * i as f64 / n as f64))
                    .collect()
            } else {
                let n = match e.curve {
                    crate::model::Curve::Line { .. } => 1,
                    _ => k * curve::sample_count(&e.curve, e.range),
                };
                curve::sample(&e.curve, e.range, n)
                    .into_iter()
                    .map(|q| {
                        let d = q - p.o;
                        [d.dot(p.x), d.dot(p.y)]
                    })
                    .collect()
            };
            if !c.forward {
                s.reverse();
            }
            if !pts.is_empty() {
                s.remove(0);
            }
            pts.extend(s);
        }
        out.push(pts);
    }
    Some(out)
}

/// Where a face's boundary crosses itself, in its parameter plane: a face
/// folded over, which is how a mesh whose topology differs from the exact
/// model's shows up once its corners sit at their exact positions (a sliver
/// the exact geometry would have removed). Loops may touch (a pinch at a
/// shared vertex), but not cross.
pub(crate) fn face_crossing(b: &Brep, fi: usize) -> Option<String> {
    // Chords of nearly tangent arcs can cross where the arcs do not, so a
    // crossing found at the usual sampling is confirmed at a finer one.
    let found = |k: usize| -> Option<String> {
        let lines = face_polylines(b, fi, k)?;
        crossing_in(&lines).map(|(a, bb)| {
            format!(
                "face {fi}: boundary crosses itself near ({:.4}, {:.4}) and ({:.4}, {:.4})",
                a[0], a[1], bb[0], bb[1]
            )
        })
    };
    if !small_face_may_cross(b, fi) {
        return None;
    }
    found(1)?;
    found(8)
}

/// Whether [`face_crossing`] can find a crossing on face `fi`: `false`
/// only for a planar face of a few straight edges whose boundary does not
/// cross itself by `crossing_in`'s own test (the same points, tolerance
/// and pairs of segments), made with no allocation. On a faceted model
/// of hundreds of thousands of such faces, building the general test's
/// polylines and sorting their segments was a fifth of reconstruction.
/// A crossing found, or any other face, is `true`: the general test then
/// decides, and says where.
fn small_face_may_cross(b: &Brep, fi: usize) -> bool {
    const SMALL: usize = 16;
    let f = &b.faces[fi];
    let Some(p) = face_param(f) else {
        return true;
    };
    let coedges = || f.loops.iter().flat_map(|l| &l.coedges);
    if p.has_pcurves()
        || f.loops.len() > SMALL
        || coedges().count() > SMALL
        || coedges().any(|c| {
            !matches!(
                b.edges[c.edge as usize].curve,
                crate::model::Curve::Line { .. }
            )
        })
    {
        return true;
    }
    // `face_polylines(.., 1)`: per coedge on a line, its two ends in plane
    // coordinates in loop order, the first dropped after the loop's
    // first coedge. Segments are (start, end, loop, index in loop).
    let mut segs = [([0.0f64; 2], [0.0f64; 2], 0usize, 0usize); 2 * SMALL];
    let mut ns = 0;
    let mut lens = [0usize; SMALL];
    let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
    for (li, lp) in f.loops.iter().enumerate() {
        let mut pts = [[0.0f64; 2]; SMALL + 1];
        let mut n = 0;
        for c in &lp.coedges {
            let e = &b.edges[c.edge as usize];
            // `curve::sample(.., 1)`'s parameters, to the bit.
            let t = |i: usize| e.range[0] + (e.range[1] - e.range[0]) * i as f64 / 1.0;
            let mut s = [t(0), t(1)].map(|t| {
                let d = curve::eval(&e.curve, t) - p.o;
                [d.dot(p.x), d.dot(p.y)]
            });
            if !c.forward {
                s.reverse();
            }
            let from = usize::from(n > 0);
            for q in &s[from..] {
                pts[n] = *q;
                n += 1;
            }
        }
        lens[li] = n;
        for i in 0..n {
            let (a, c) = (pts[i], pts[(i + 1) % n]);
            for d in 0..2 {
                lo[d] = lo[d].min(a[d]);
                hi[d] = hi[d].max(a[d]);
            }
            segs[ns] = (a, c, li, i);
            ns += 1;
        }
    }
    let segs = &segs[..ns];
    let ext = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    let eps = 1e-9 * ext;
    let orient = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };
    let minx = |s: &([f64; 2], [f64; 2], usize, usize)| s.0[0].min(s.1[0]);
    for (i, si) in segs.iter().enumerate() {
        for sj in &segs[i + 1..] {
            // `crossing_in` sorts by the smaller x and stops at the first
            // segment starting past the earlier one's end: the pairs it
            // compares, in either order.
            let (first, second) = if minx(sj) < minx(si) {
                (sj, si)
            } else {
                (si, sj)
            };
            if minx(second) > first.0[0].max(first.1[0]) + eps {
                continue;
            }
            let ((p1, p2, li, ii), (q1, q2, lj, jj)) = (*si, *sj);
            if li == lj {
                let n = lens[li];
                if (ii + 1) % n == jj || (jj + 1) % n == ii || ii == jj {
                    continue;
                }
            }
            let lp = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2)).sqrt();
            let lq = ((q2[0] - q1[0]).powi(2) + (q2[1] - q1[1]).powi(2)).sqrt();
            let (o1, o2) = (orient(p1, p2, q1), orient(p1, p2, q2));
            let (o3, o4) = (orient(q1, q2, p1), orient(q1, q2, p2));
            let (ep, eq) = (eps * lp, eps * lq);
            if ((o1 > ep && o2 < -ep) || (o1 < -ep && o2 > ep))
                && ((o3 > eq && o4 < -eq) || (o3 < -eq && o4 > eq))
            {
                return true;
            }
        }
    }
    false
}

/// Where a corner of a face lies inside another edge of the same face,
/// within `tol` in space: a boundary touching itself in the middle of an
/// edge. Corners that coincide (a pinch) are allowed.
///
/// Manifold keeps two bodies that touch along a line apart with duplicated
/// vertices, and after a rotation its rounding can leave the duplicates on
/// different sides of each other, so that the mesh joins faces the exact
/// model only pinches. Once the vertices sit at their exact positions, a
/// corner lies on an edge of its own face. Counted by its topology the
/// B-rep is still a closed 2-manifold, but a reader that joins by position
/// sees a self-intersecting wire: the rotated Menger sponge
/// (`example024.scad`) read back from OCCT as an open shell with 60 free
/// edges because of this. The test is in space, not in the parameter
/// plane, whose `u` (an angle) and `v` (a length) are not comparable.
pub(crate) fn face_touch(b: &Brep, fi: usize, tol: f64) -> Option<String> {
    if !small_face_touches(b, &b.faces[fi], tol) {
        return None;
    }
    face_touch_search(b, fi, tol)
}

/// [`face_touch`]'s search on any face.
fn face_touch_search(b: &Brep, fi: usize, tol: f64) -> Option<String> {
    let f = &b.faces[fi];
    // (vertex, loop) and (edge, loop) pairs.
    let mut verts: Vec<(u32, usize)> = Vec::new();
    let mut edges: Vec<(u32, usize)> = Vec::new();
    for (li, lp) in f.loops.iter().enumerate() {
        for c in &lp.coedges {
            let e = &b.edges[c.edge as usize];
            verts.push((e.start, li));
            verts.push((e.end, li));
            edges.push((c.edge, li));
        }
    }
    verts.sort_unstable();
    verts.dedup();
    edges.sort_unstable();
    edges.dedup();
    if verts.len() < 3 {
        return None;
    }
    let at = |v: u32| V::from(b.vertices[v as usize]);
    // The corners sorted by x, so that each edge looks only at those in
    // its box: faces of thousands of corners are common (text, `$fn`).
    let mut by_x: Vec<(f64, u32, usize)> = verts.iter().map(|&(v, l)| (at(v).x, v, l)).collect();
    by_x.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
    // One buffer of samples for all the face's edges (faceted faces are
    // a few lines each, and there are hundreds of thousands of them).
    let mut pts: Vec<V> = Vec::new();
    for &(ei, el) in &edges {
        let e = &b.edges[ei as usize];
        let n = curve::sample_count(&e.curve, e.range);
        pts.clear();
        pts.extend((0..=n).map(|i| {
            curve::eval(
                &e.curve,
                e.range[0] + (e.range[1] - e.range[0]) * i as f64 / n as f64,
            )
        }));
        let pts = &pts[..];
        let mut lo = pts[0];
        let mut hi = pts[0];
        for p in pts {
            lo = v(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = v(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        // The samples' box misses the curve between them by at most the
        // sagitta, which a twentieth of the spacing bounds at this
        // sampling; the pad covers that and the tolerance.
        let pad = tol + 0.05 * (hi - lo).len() / n as f64;
        let (s, t) = (at(e.start), at(e.end));
        let first = by_x.partition_point(|&(x, _, _)| x < lo.x - pad);
        for &(x, vi, vl) in &by_x[first..] {
            if x > hi.x + pad {
                break;
            }
            if vi == e.start || vi == e.end {
                continue;
            }
            let p = at(vi);
            if p.y < lo.y - pad
                || p.z < lo.z - pad
                || p.y > hi.y + pad
                || p.z > hi.z + pad
                || (p - s).len() <= tol
                || (p - t).len() <= tol
            {
                continue;
            }
            if within(&e.curve, e.range, pts, p, tol) {
                return Some(format!(
                    "face {fi}: boundary touches itself, a corner at ({:.4}, {:.4}, {:.4}) lies on edge {ei} ({})",
                    p.x,
                    p.y,
                    p.z,
                    if vl == el {
                        "same loop"
                    } else {
                        "another loop"
                    }
                ));
            }
        }
    }
    None
}

/// Whether [`face_touch`] can find anything on `f`: `false` only for a
/// face of a few straight edges where no corner lies on another edge, by
/// [`face_touch`]'s own test (the same samples, boxes and distances),
/// made with no sorting and no allocation. Faceted and polygonal models
/// are hundreds of thousands of such faces, and setting up the general
/// search (corner and edge lists, sorted, deduplicated) took most of the
/// check's time on them. Any other face, or a touch found, is `true`:
/// the general search then decides, and names what it finds as before.
fn small_face_touches(b: &Brep, f: &crate::model::Face, tol: f64) -> bool {
    const SMALL: usize = 16;
    let n_coedges: usize = f.loops.iter().map(|l| l.coedges.len()).sum();
    if n_coedges > SMALL
        || f.loops.iter().flat_map(|l| &l.coedges).any(|c| {
            !matches!(
                b.edges[c.edge as usize].curve,
                crate::model::Curve::Line { .. }
            )
        })
    {
        return true;
    }
    let at = |v: u32| V::from(b.vertices[v as usize]);
    // Every corner, with repeats (a repeat finds what the first finds).
    let corners = || {
        f.loops.iter().flat_map(|l| &l.coedges).flat_map(|c| {
            let e = &b.edges[c.edge as usize];
            [e.start, e.end]
        })
    };
    let mut distinct = [u32::MAX; 2 * SMALL];
    let mut nd = 0;
    for v in corners() {
        if !distinct[..nd].contains(&v) {
            distinct[nd] = v;
            nd += 1;
        }
    }
    if nd < 3 {
        return false;
    }
    for c in f.loops.iter().flat_map(|l| &l.coedges) {
        let e = &b.edges[c.edge as usize];
        // `sample_count` of a line.
        const N: usize = 4;
        let pts: [V; N + 1] = std::array::from_fn(|i| {
            curve::eval(
                &e.curve,
                e.range[0] + (e.range[1] - e.range[0]) * i as f64 / N as f64,
            )
        });
        let mut lo = pts[0];
        let mut hi = pts[0];
        for p in &pts {
            lo = v(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = v(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let pad = tol + 0.05 * (hi - lo).len() / N as f64;
        let (s, t) = (at(e.start), at(e.end));
        for &vi in &distinct[..nd] {
            if vi == e.start || vi == e.end {
                continue;
            }
            let p = at(vi);
            // `face_touch`'s window in x, then its box and its ends.
            if p.x < lo.x - pad
                || p.x > hi.x + pad
                || p.y < lo.y - pad
                || p.z < lo.z - pad
                || p.y > hi.y + pad
                || p.z > hi.z + pad
                || (p - s).len() <= tol
                || (p - t).len() <= tol
            {
                continue;
            }
            if within(&e.curve, e.range, &pts, p, tol) {
                return true;
            }
        }
    }
    false
}

/// Whether `p` is within `tol` of a curve sampled as `pts` over `range`,
/// by [`distance_to_curve`].
///
/// For a line, a corner farther than `tol` from the whole infinite line
/// (plus a margin far above the rounding of either distance) is not
/// within `tol` of any point the search evaluates, so the answer is
/// known without it: the search costs 160 evaluations, and the touch
/// checks ask it of most corners of a faceted model's thin faces (their
/// boxes overlap their edges'). Any other case is the search's answer.
fn within(c: &crate::model::Curve, range: [f64; 2], pts: &[V], p: V, tol: f64) -> bool {
    if let crate::model::Curve::Line { origin, direction } = c {
        let (o, d) = (V::from(*origin), V::from(*direction));
        let dl = d.len();
        if dl > 0.0 {
            let off = (p - o).cross(d).len() / dl;
            let reach = o.len() + p.len() + dl * range[0].abs().max(range[1].abs());
            if off > tol + 1e-12 * reach {
                return false;
            }
        }
    }
    distance_to_curve(c, range, pts, p) <= tol
}

/// The distance from `p` to a curve sampled as `pts` over `range`: from
/// the nearest sample's neighbourhood, refined by golden-section search.
fn distance_to_curve(c: &crate::model::Curve, range: [f64; 2], pts: &[V], p: V) -> f64 {
    let n = pts.len() - 1;
    let k = (0..=n)
        .min_by(|&i, &j| (pts[i] - p).len().total_cmp(&(pts[j] - p).len()))
        .expect("samples");
    let step = (range[1] - range[0]) / n as f64;
    let (mut a, mut z) = (
        range[0] + step * k.saturating_sub(1) as f64,
        range[0] + step * (k + 1).min(n) as f64,
    );
    let d = |t: f64| (curve::eval(c, t) - p).len();
    let g = 0.5 * (5f64.sqrt() - 1.0);
    for _ in 0..80 {
        let (m1, m2) = (z - g * (z - a), a + g * (z - a));
        if d(m1) < d(m2) {
            z = m2;
        } else {
            a = m1;
        }
    }
    d(0.5 * (a + z)).min((pts[k] - p).len())
}

fn crossing_in(lines: &[Vec<[f64; 2]>]) -> Option<([f64; 2], [f64; 2])> {
    // Segments (start, end, loop, index), closed per loop.
    let mut segs: Vec<([f64; 2], [f64; 2], usize, usize)> = Vec::new();
    let mut lo = [f64::INFINITY; 2];
    let mut hi = [f64::NEG_INFINITY; 2];
    for (li, l) in lines.iter().enumerate() {
        let n = l.len();
        for i in 0..n {
            let (a, c) = (l[i], l[(i + 1) % n]);
            for d in 0..2 {
                lo[d] = lo[d].min(a[d]);
                hi[d] = hi[d].max(a[d]);
            }
            segs.push((a, c, li, i));
        }
    }
    let ext = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(1e-300);
    let eps = 1e-9 * ext;
    let len_of = |li: usize| lines[li].len();
    // A stable sort by the smaller x, made as an unstable one by that x
    // and then the position (the same order), with the key computed once:
    // curved faces have thousands of segments, and sorting them was a
    // third of this check.
    let mut order: Vec<(f64, u32)> = segs
        .iter()
        .enumerate()
        .map(|(i, s)| (s.0[0].min(s.1[0]), i as u32))
        .collect();
    order.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let segs: Vec<_> = order.iter().map(|&(_, i)| segs[i as usize]).collect();
    let orient = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };
    for i in 0..segs.len() {
        let (p1, p2, li, ii) = segs[i];
        let maxx = p1[0].max(p2[0]);
        for &(q1, q2, lj, jj) in &segs[i + 1..] {
            if q1[0].min(q2[0]) > maxx + eps {
                break;
            }
            if li == lj {
                let n = len_of(li);
                if (ii + 1) % n == jj || (jj + 1) % n == ii || ii == jj {
                    continue;
                }
            }
            let lp = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2)).sqrt();
            let lq = ((q2[0] - q1[0]).powi(2) + (q2[1] - q1[1]).powi(2)).sqrt();
            let (o1, o2) = (orient(p1, p2, q1), orient(p1, p2, q2));
            let (o3, o4) = (orient(q1, q2, p1), orient(q1, q2, p2));
            let (ep, eq) = (eps * lp, eps * lq);
            if ((o1 > ep && o2 < -ep) || (o1 < -ep && o2 > ep))
                && ((o3 > eq && o4 < -eq) || (o3 < -eq && o4 > eq))
            {
                return Some((p1, q1));
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Coedge, Curve, Edge, Face, Frame, Loop, Surface};

    /// One planar face on z = 0 bounded by `loops` (each a polygon of
    /// points, a vertex of its own per point, so a repeated point is two
    /// vertices at one position, as Manifold keeps touching bodies apart),
    /// with straight edges.
    fn polygon_face(loops: &[Vec<[f64; 2]>]) -> Brep {
        let mut b = Brep::default();
        let mut face_loops = Vec::new();
        for (li, pts) in loops.iter().enumerate() {
            let first = b.vertices.len() as u32;
            b.vertices.extend(pts.iter().map(|p| [p[0], p[1], 0.0]));
            let mut coedges = Vec::new();
            for k in 0..pts.len() {
                let (s, e) = (first + k as u32, first + ((k + 1) % pts.len()) as u32);
                let (p, q) = (b.vertices[s as usize], b.vertices[e as usize]);
                let d = [q[0] - p[0], q[1] - p[1], 0.0];
                let l = (d[0] * d[0] + d[1] * d[1]).sqrt();
                coedges.push(Coedge {
                    edge: b.edges.len() as u32,
                    forward: true,
                    pcurve: None,
                });
                b.edges.push(Edge {
                    start: s,
                    end: e,
                    curve: Curve::Line {
                        origin: p,
                        direction: [d[0] / l, d[1] / l, 0.0],
                    },
                    range: [0.0, l],
                    seam: false,
                });
            }
            face_loops.push(Loop {
                coedges,
                outer: li == 0,
            });
        }
        b.faces.push(Face {
            surface: Surface::Plane {
                origin: [0.0; 3],
                normal: [0.0, 0.0, 1.0],
            },
            frame: Frame {
                origin: [0.0; 3],
                z: [0.0, 0.0, 1.0],
                x: [1.0, 0.0, 0.0],
            },
            ref_radius: 0.0,
            same_sense: true,
            loops: face_loops,
            faceted: true,
        });
        b
    }

    /// Polygons of 4 and 5 corners on a 3 × 3 grid, sheared a little so
    /// that not every edge is axis-aligned: bow ties, corners on other
    /// edges, repeated points, plain quads.
    fn polygons() -> Vec<Vec<[f64; 2]>> {
        let grid: Vec<[f64; 2]> = (0..9)
            .map(|i| {
                let (x, y) = ((i % 3) as f64, (i / 3) as f64);
                [x + 0.125 * y, y]
            })
            .collect();
        let mut out = Vec::new();
        for a in 0..9 {
            for b in 0..9 {
                for c in 0..9 {
                    for d in 0..9 {
                        if a == b || b == c || c == d || d == a {
                            continue;
                        }
                        out.push(vec![grid[a], grid[b], grid[c], grid[d]]);
                        let e = (a + b + c + d) % 9;
                        if e != d && e != a {
                            out.push(vec![grid[a], grid[b], grid[c], grid[d], grid[e]]);
                        }
                    }
                }
            }
        }
        out
    }

    /// The quick tests in front of the crossing and touch searches answer
    /// as the searches do on the small straight-edged faces they take: a
    /// quick "no" is never a missed crossing or touch, and a "yes" is one
    /// the search finds (and names).
    #[test]
    fn the_quick_face_checks_agree_with_the_searches() {
        let hole = vec![[0.75, 0.75], [1.25, 0.75], [1.0, 1.25]];
        let (mut crossings, mut touches, mut faces) = (0, 0, 0);
        for (i, outer) in polygons().into_iter().enumerate() {
            let mut loops = vec![outer];
            if i % 7 == 0 {
                loops.push(hole.clone());
            }
            let b = polygon_face(&loops);
            faces += 1;
            let search = face_polylines(&b, 0, 1).and_then(|l| crossing_in(&l));
            assert_eq!(
                small_face_may_cross(&b, 0),
                search.is_some(),
                "crossing on {loops:?}"
            );
            crossings += usize::from(search.is_some());
            for tol in [1e-6, 0.05] {
                let search = face_touch_search(&b, 0, tol);
                assert_eq!(
                    small_face_touches(&b, &b.faces[0], tol),
                    search.is_some(),
                    "touch at {tol} on {loops:?}"
                );
                touches += usize::from(search.is_some());
            }
        }
        // Both answers, many times each.
        assert!(
            crossings > faces / 10 && crossings < faces - faces / 10,
            "{crossings} of {faces}"
        );
        assert!(touches > faces / 10, "{touches} of {faces}");
    }

    /// The line shortcut in front of the distance search answers as the
    /// search does: on lines of every direction and length, far from the
    /// origin and near it, with corners on the line, at the tolerance
    /// from it (either side), and far off, at a tolerance as small as the
    /// rounding and at the usual ones.
    #[test]
    fn the_line_shortcut_agrees_with_the_distance_search() {
        let mut seed: u64 = 0x9e37_79b9_7f4a_7c15;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let (mut near, mut far, mut cases) = (0, 0, 0);
        for i in 0..20_000 {
            let scale = [1e-3, 1.0, 1e3, 1e6][i % 4];
            let o = v(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5) * scale;
            let d = v(rnd() - 0.5, rnd() - 0.5, rnd() - 0.5).norm();
            let len = scale * (0.001 + rnd());
            let curve = crate::model::Curve::Line {
                origin: o.arr(),
                direction: d.arr(),
            };
            let range = [0.0, len];
            const N: usize = 4;
            let pts: Vec<V> = (0..=N)
                .map(|k| curve::eval(&curve, len * k as f64 / N as f64))
                .collect();
            let tol = [1e-12 * scale, 1e-9 * scale, 1e-6, 1e-3 * scale][(i / 4) % 4];
            // A point along the line (inside and beyond the range),
            // pushed off it by a multiple of the tolerance.
            let along = o + d * (len * (1.4 * rnd() - 0.2));
            let side = d.perp().norm();
            let off = [0.0, 0.5, 1.0 - 1e-9, 1.0, 1.0 + 1e-9, 2.0, 1e3][(i / 16) % 7];
            let p = along + side * (off * tol);
            let search = distance_to_curve(&curve, range, &pts, p) <= tol;
            assert_eq!(
                within(&curve, range, &pts, p, tol),
                search,
                "{o:?} {d:?} {len} {p:?} {tol}"
            );
            cases += 1;
            near += usize::from(search);
            far += usize::from(!search);
        }
        assert!(near > cases / 10 && far > cases / 10, "{near} {far}");
    }

    /// `validate` makes the crossing check unless the B-rep says that
    /// reconstruction made it: a bow tie is found, and not looked for on a
    /// B-rep whose report vouches for it.
    #[test]
    fn the_crossing_check_is_skipped_only_when_reconstruction_made_it() {
        // A pentagon one of whose edges crosses another, and the same
        // face reversed: every edge used twice, in opposite directions,
        // one shell of genus 0, so the crossing check is reached.
        let mut b = polygon_face(&[vec![
            [0.0, 0.0],
            [4.0, 0.0],
            [4.0, 3.0],
            [1.0, -1.0],
            [0.0, 3.0],
        ]]);
        let mut back = b.faces[0].clone();
        for l in &mut back.loops {
            l.coedges.reverse();
            for c in &mut l.coedges {
                c.forward = !c.forward;
            }
        }
        b.faces.push(back);
        b.shells.push(crate::model::Shell {
            faces: vec![0, 1],
            void: false,
        });
        let errors = |b: &Brep| validate(b, 1e-6).errors;
        let found = errors(&b);
        assert!(
            found.iter().all(|e| e.contains("crosses itself")) && found.len() == 2,
            "{found:?}"
        );
        // Vouched for, it is not looked for (the next check, of the
        // shells' volumes, then finds the faces enclose none).
        b.report.faces_uncrossed = true;
        let found = errors(&b);
        assert!(
            found.iter().all(|e| !e.contains("crosses itself")) && !found.is_empty(),
            "{found:?}"
        );
    }
}
