//! The B-rep while it is being built: mutable, with the surfaces and mesh
//! chains each part came from, which the later passes (seams, parameter
//! curves) need and the published [`crate::Brep`] does not carry.

use crate::curve;
use crate::math::*;
use crate::model::{BSpline, Curve};
use crate::surf::{Param, Surf};

#[derive(Clone, Debug)]
pub(crate) struct TEdge {
    pub v0: usize,
    pub v1: usize,
    pub curve: Curve,
    pub range: [f64; 2],
    /// The faces on either side (the face using it forward first). Equal
    /// for a seam.
    pub faces: [usize; 2],
    /// The mesh points the edge was built from (empty for seams).
    pub chain: Vec<V>,
    pub seam: bool,
    /// Largest distance of samples from either surface.
    pub dev: f64,
    /// How far the mesh chain stands off the curve, across the surfaces
    /// (see `Report::max_chain_deviation`); zero where it does not count.
    pub chain_dev: f64,
}

impl TEdge {
    pub fn closed(&self) -> bool {
        self.v0 == self.v1
    }

    /// Whether the edge is a spike of no width: its two corners are one
    /// point (within `tol`), and its mesh chain runs out and back along
    /// itself, point for point. At a vertex where four blends meet, all
    /// tangent (the end of an L-bracket rounded in two passes, where a
    /// convex arc's torus touches the wall at the corner of its concave
    /// blend), the polygonal tools leave the wall a needle between two
    /// strips of the torus, its two sides one chain of mesh edges out and
    /// back. The wall meets the torus in a figure of eight through that
    /// corner, and the edge fitted between the chain's ends, which solve
    /// to the same point, ran round a whole lobe of it: a face boundary
    /// OCCT rejects (`BadOrientationOfSubshape`, fillet corpus seed 2
    /// models 19 and 1009). A needle bounds nothing; the edge goes as a
    /// short one does.
    pub fn retraced(&self, verts: &[V], tol: f64) -> bool {
        let c = &self.chain;
        if c.len() < 3 || (verts[self.v0] - verts[self.v1]).len() >= tol {
            return false;
        }
        (0..c.len() / 2).all(|i| (c[i] - c[c.len() - 1 - i]).len() < tol)
    }

    /// Whether the edge (of length `l`) is a contact the exact model does
    /// not have: its two corners are one point (within `tol`), and the
    /// curve fitted between them is a wisp of under a thousand times
    /// `tol` while the mesh's chain runs on for twice that or more. Where
    /// two faces touch beside a third at one point (a rotated box's
    /// rounded edge and a cone's rim, fillet corpus seed 2 model 192) the
    /// polygons meet along a chain 0.27 long whose ends solve to the same
    /// corner; the B-spline fitted between them was 3.5e-7 long, and the
    /// validator rightly called it degenerate. A closed curve the model has
    /// is followed by its fit all round.
    pub fn stray(&self, verts: &[V], l: f64, tol: f64) -> bool {
        let c = &self.chain;
        if c.len() < 2 || (verts[self.v0] - verts[self.v1]).len() >= tol || l >= 1e3 * tol {
            return false;
        }
        let along: f64 = c.windows(2).map(|w| (w[1] - w[0]).len()).sum();
        along >= 2.0 * l
    }
}

#[derive(Clone, Debug)]
pub(crate) struct TFace {
    pub surf: Surf,
    pub same_sense: bool,
    pub faceted: bool,
    /// Loops of (edge, forward).
    pub loops: Vec<Vec<(usize, bool)>>,
    /// Set by the parametrisation pass.
    pub param: Option<Param>,
    /// Per loop, per coedge: the parameter-space curve (edge direction).
    pub pcurves: Vec<Vec<Option<BSpline<2>>>>,
    pub outer: Vec<bool>,
    /// A torus face's mesh triangles (empty for other surfaces): which
    /// tube and axis angles the face covers decides where its frame puts
    /// the cuts of both periodic coordinates.
    pub tris: Vec<[V; 3]>,
    /// The input triangles the face was built from, so that a failure can
    /// be reported where it is (`Report::face_triangles`).
    pub source: Vec<u32>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Topo {
    pub verts: Vec<V>,
    /// Vertices a seam ends at: they may no longer move.
    pub pinned: Vec<bool>,
    pub edges: Vec<TEdge>,
    pub faces: Vec<TFace>,
    pub notes: Vec<String>,
}

/// How far `p` is from the closed-form edge `e` (a line, circle or
/// ellipse over its range); `None` for a B-spline.
fn distance_to(e: &TEdge, p: V) -> Option<f64> {
    if matches!(e.curve, Curve::BSpline(_)) {
        return None;
    }
    let mut t = curve::param_of(&e.curve, p);
    if !matches!(e.curve, Curve::Line { .. }) {
        // Angles come in (-π, π]; the range may be anywhere.
        while t < e.range[0] {
            t += TAU;
        }
        while t - TAU >= e.range[0] {
            t -= TAU;
        }
    }
    let ends = [e.range[0], e.range[1]].map(|r| (curve::eval(&e.curve, r) - p).len());
    Some(if t <= e.range[1] {
        (curve::eval(&e.curve, t) - p)
            .len()
            .min(ends[0])
            .min(ends[1])
    } else {
        ends[0].min(ends[1])
    })
}

/// Whether two edges with the same ends run along each other within
/// `tol`: samples of each lie on the other where it has a closed form
/// (two B-splines are not compared, and so never coincide here).
fn coincide(a: &TEdge, b: &TEdge, tol: f64) -> bool {
    let along = |x: &TEdge, y: &TEdge| -> Option<bool> {
        let pts = curve::sample(&y.curve, y.range, 16);
        pts.iter()
            .map(|&p| distance_to(x, p).map(|d| d < tol))
            .try_fold(true, |acc, ok| ok.map(|ok| acc && ok))
    };
    match (along(a, b), along(b, a)) {
        (Some(x), Some(y)) => x && y,
        (Some(x), None) | (None, Some(x)) => x,
        (None, None) => false,
    }
}

impl Topo {
    pub fn add_vertex(&mut self, p: V) -> usize {
        self.verts.push(p);
        self.pinned.push(false);
        self.verts.len() - 1
    }

    /// The vertex a coedge starts at.
    pub fn start(&self, (e, fwd): (usize, bool)) -> usize {
        let e = &self.edges[e];
        if fwd { e.v0 } else { e.v1 }
    }

    /// Splits edge `e` at curve parameter `t` (strictly inside its range)
    /// with a new vertex at `p`, updating every loop that uses it. Returns
    /// the new vertex.
    pub fn split_edge(&mut self, e: usize, t: f64, p: V) -> usize {
        let nv = self.add_vertex(p);
        let mut second = self.edges[e].clone();
        second.v0 = nv;
        second.range[0] = t;
        second.chain = Vec::new();
        self.edges[e].v1 = nv;
        self.edges[e].range[1] = t;
        let ne = self.edges.len();
        let near = self.edges[e].faces;
        self.edges.push(second);
        // The faces the edge names first: it is used once on each side,
        // so finding both uses there is the usual case, and rebuilding
        // every loop of the model for each split made seam placement
        // quadratic (thousands of seams on the BOSL2 fractal tree). Any
        // other count, or names out of date, and every loop is searched.
        let mut found = 0;
        let mut seen = [usize::MAX; 2];
        for (k, f) in near.into_iter().enumerate() {
            if f < self.faces.len() && !seen.contains(&f) {
                seen[k] = f;
                found += Self::split_uses(&mut self.faces[f], e, ne);
            }
        }
        if found != 2 {
            for (fi, f) in self.faces.iter_mut().enumerate() {
                if !seen.contains(&fi) {
                    Self::split_uses(f, e, ne);
                }
            }
        }
        nv
    }

    /// Replaces each use of edge `e` in `f`'s loops by `e` then `ne` (in
    /// the use's direction). Returns how many uses there were.
    fn split_uses(f: &mut TFace, e: usize, ne: usize) -> usize {
        let mut n = 0;
        for lp in &mut f.loops {
            if !lp.iter().any(|&(x, _)| x == e) {
                continue;
            }
            let mut out = Vec::with_capacity(lp.len() + 1);
            for &(x, fwd) in lp.iter() {
                if x == e {
                    n += 1;
                    if fwd {
                        out.push((e, true));
                        out.push((ne, true));
                    } else {
                        out.push((ne, false));
                        out.push((e, false));
                    }
                } else {
                    out.push((x, fwd));
                }
            }
            *lp = out;
        }
        n
    }

    /// Splits line edge `e` at existing vertex `v`, which lies on it
    /// strictly inside its range, updating every loop that uses it. Returns
    /// the new edge (the part from `v` to the old end).
    fn split_line_at(&mut self, e: usize, v: usize) -> usize {
        let t = curve::param_of(&self.edges[e].curve, self.verts[v]);
        let mut second = self.edges[e].clone();
        second.v0 = v;
        second.range[0] = t;
        second.chain = Vec::new();
        self.edges[e].v1 = v;
        self.edges[e].range[1] = t;
        self.edges[e].chain = Vec::new();
        let ne = self.edges.len();
        self.edges.push(second);
        for f in &mut self.faces {
            for lp in &mut f.loops {
                let mut out = Vec::with_capacity(lp.len() + 1);
                for &(x, fwd) in lp.iter() {
                    if x == e {
                        if fwd {
                            out.push((e, true));
                            out.push((ne, true));
                        } else {
                            out.push((ne, false));
                            out.push((e, false));
                        }
                    } else {
                        out.push((x, fwd));
                    }
                }
                *lp = out;
            }
        }
        ne
    }

    /// Removes needles: a face whose loop runs out along a straight edge
    /// and straight back along another, a spike of no width between two
    /// other faces. Where a blend ends square on another blend (a block on
    /// a plate rounded in two passes: the vertical edge's convex blend
    /// ends at the top of the concave blend round the block's foot), three
    /// faces meet along one line: the step the convex blend's end leaves,
    /// the concave blend tangent to the wall, and the wall. The mesh gives
    /// the wall a strip of no width along it, and the wall's boundary
    /// came out touching itself (fillet corpus seed 2 models 1404, 1823,
    /// 1899): a corner of the wall on another of its edges, which the
    /// validator refuses. The two faces either side of the needle meet
    /// along its line instead: the longer edge is split at the shorter's
    /// far corner, and the needle's two edges become one between them.
    /// Returns how many were removed.
    pub fn remove_needles(&mut self, tol: f64) -> usize {
        let mut gone = vec![false; self.edges.len()];
        let mut n = 0;
        let end = |t: &Topo, (e, fwd): (usize, bool)| {
            let e = &t.edges[e];
            if fwd { e.v1 } else { e.v0 }
        };
        let users = |t: &Topo, e: usize| -> Vec<usize> {
            (0..t.faces.len())
                .filter(|&f| t.faces[f].loops.iter().flatten().any(|c| c.0 == e))
                .collect()
        };
        let line =
            |t: &Topo, e: usize| matches!(t.edges[e].curve, Curve::Line { .. }) && !t.edges[e].seam;
        while n < 1000 {
            // The first needle in face, loop and coedge order.
            let mut found = None;
            'search: for f in 0..self.faces.len() {
                for l in 0..self.faces[f].loops.len() {
                    let lp = &self.faces[f].loops[l];
                    let k = lp.len();
                    if k < 3 {
                        continue;
                    }
                    for i in 0..k {
                        let (c1, c2) = (lp[i], lp[(i + 1) % k]);
                        if c1.0 == c2.0 || !line(self, c1.0) || !line(self, c2.0) {
                            continue;
                        }
                        let (a, b, c) = (self.start(c1), end(self, c1), end(self, c2));
                        if self.start(c2) != b || a == b || b == c {
                            continue;
                        }
                        let (d1, d2) =
                            (self.verts[b] - self.verts[a], self.verts[c] - self.verts[b]);
                        let (l1, l2) = (d1.len(), d2.len());
                        if l1 <= tol || l2 <= tol {
                            continue;
                        }
                        // Back along the same line: antiparallel, and the
                        // far corner of each on the other's line.
                        if d1.dot(d2) >= 0.0
                            || d1.cross(d2).len() > 1e-9 * l1 * l2
                            || (self.verts[c] - self.verts[a]).cross(d1).len() > tol * l1
                        {
                            continue;
                        }
                        let (u1, u2) = (users(self, c1.0), users(self, c2.0));
                        let other = |u: &[usize]| -> Option<usize> {
                            match u {
                                [x, y] if *x == f => Some(*y),
                                [x, y] if *y == f => Some(*x),
                                _ => None,
                            }
                        };
                        if let (Some(g1), Some(g2)) = (other(&u1), other(&u2))
                            && g1 != g2
                        {
                            found = Some((f, c1, c2, a, c, l1, l2));
                            break 'search;
                        }
                    }
                }
            }
            let Some((f, c1, c2, a, c, l1, l2)) = found else {
                break;
            };
            // The needle as two coedges of `f`, `x` out from `p` to `q` and
            // `y` back, each the whole of its edge.
            let (x, y) = if (l1 - l2).abs() <= tol {
                if a != c {
                    // The far corners are one point: join them.
                    for ed in &mut self.edges {
                        if ed.v0 == c {
                            ed.v0 = a;
                        }
                        if ed.v1 == c {
                            ed.v1 = a;
                        }
                    }
                }
                (c1, c2)
            } else if l2 > l1 {
                let e2 = c2.0;
                let ne = self.split_line_at(e2, a);
                gone.push(false);
                // `c2` ran from the tip back past `a`: its part up to `a`
                // is the needle's return.
                (c1, if c2.1 { (e2, true) } else { (ne, false) })
            } else {
                let e1 = c1.0;
                let ne = self.split_line_at(e1, c);
                gone.push(false);
                (if c1.1 { (ne, true) } else { (e1, false) }, c2)
            };
            // `f` drops both; the face beyond `y` takes `x` in its place,
            // running the way `y`'s other use did (the way `x` runs in
            // `f`).
            for lp in &mut self.faces[f].loops {
                lp.retain(|&co| co != x && co != y);
            }
            for g in 0..self.faces.len() {
                if g == f {
                    continue;
                }
                for lp in &mut self.faces[g].loops {
                    for co in lp.iter_mut() {
                        if co.0 == y.0 {
                            *co = x;
                        }
                    }
                }
            }
            self.faces[f].loops.retain(|lp| !lp.is_empty());
            gone[y.0] = true;
            n += 1;
        }
        if n == 0 {
            return 0;
        }
        let mut renum = vec![usize::MAX; self.edges.len()];
        let mut kept = Vec::with_capacity(self.edges.len());
        for (i, ed) in std::mem::take(&mut self.edges).into_iter().enumerate() {
            if !gone[i] {
                renum[i] = kept.len();
                kept.push(ed);
            }
        }
        self.edges = kept;
        for f in &mut self.faces {
            for lp in &mut f.loops {
                for c in lp.iter_mut() {
                    c.0 = renum[c.0];
                }
            }
        }
        n
    }

    /// Collapses open edges shorter than `tol`, joining their two vertices.
    ///
    /// The mesh can leave two corners a hair apart where the exact model
    /// has one: an exact cylinder crossing the edge between two faceted
    /// planes meets it at one point, but the polygonal cylinder crosses
    /// that edge's triangles at two nearby points. Both corners solve onto
    /// the same three surfaces, i.e. the same point, leaving an edge of
    /// length ~1e-13 that readers reject (OCCT: an unorientable face).
    /// Returns how many were collapsed.
    ///
    /// A sliver triangle of the mesh whose three corners solve to two
    /// points (a face flush with a rotated one: the rotation's rounding
    /// leaves Manifold a sliver where the faces meet) has one short edge
    /// between distinct vertices and, once that collapses, another that
    /// starts and ends at the same vertex with no length. Those go too,
    /// in a second pass; the sliver is then a digon ([`Topo::drop_digons`])
    /// or loses its loop altogether.
    pub fn collapse_short_edges(&mut self, tol: f64) -> usize {
        let mut gone = vec![false; self.edges.len()];
        let mut n = 0;
        let length = |ed: &TEdge| -> f64 {
            curve::sample(&ed.curve, ed.range, 8)
                .windows(2)
                .map(|w| (w[1] - w[0]).len())
                .sum()
        };
        for e in 0..self.edges.len() {
            let ed = &self.edges[e];
            if ed.closed() || ed.seam {
                continue;
            }
            let l = length(ed);
            if l >= tol && !ed.retraced(&self.verts, tol) && !ed.stray(&self.verts, l, tol) {
                continue;
            }
            let (keep, drop) = (ed.v0, ed.v1);
            for x in &mut self.edges {
                if x.v0 == drop {
                    x.v0 = keep;
                }
                if x.v1 == drop {
                    x.v1 = keep;
                }
            }
            gone[e] = true;
            n += 1;
        }
        for e in 0..self.edges.len() {
            let ed = &self.edges[e];
            if !gone[e]
                && !ed.seam
                && ed.closed()
                && matches!(ed.curve, Curve::Line { .. })
                && length(ed) < tol
            {
                gone[e] = true;
                n += 1;
            }
        }
        if n == 0 {
            return 0;
        }
        let mut renum = vec![usize::MAX; self.edges.len()];
        let mut kept = Vec::with_capacity(self.edges.len() - n);
        for (i, ed) in std::mem::take(&mut self.edges).into_iter().enumerate() {
            if !gone[i] {
                renum[i] = kept.len();
                kept.push(ed);
            }
        }
        self.edges = kept;
        for f in &mut self.faces {
            for lp in &mut f.loops {
                lp.retain(|&(e, _)| !gone[e]);
                for c in lp.iter_mut() {
                    c.0 = renum[c.0];
                }
            }
            f.loops.retain(|lp| !lp.is_empty());
        }
        n
    }

    /// Removes faces that are two straight edges between the same two
    /// vertices: slivers with no area, which is what a thin mesh triangle
    /// becomes once its short edge collapses. Their two edges are one
    /// line, so the faces either side now share the first, and the second
    /// goes. Returns the removed faces' indices, in order.
    ///
    /// They are common on faceted models (twelve files of the stop-rule
    /// corpora had them, `issue1165.scad` among the eligible ones). OCCT
    /// reads most of them, but the face has no area and no orientation:
    /// left in, it fails our own check for an outer loop. A whole face of
    /// two curved edges that run along each other between the same two
    /// vertices goes the same way (see the comment on `lines` below).
    ///
    /// Two digons on the same two edges are a closed bubble of no volume
    /// (a sliver of the mesh between two flush faces, one of them rotated:
    /// `example017.scad`'s tabs in their slots); both go, with their edges.
    /// And a digon that is one loop of a larger face is a slit of no
    /// width in it (the same tabs, where they come through the face): the
    /// loop goes, and the faces either side share one edge.
    pub fn drop_digons(&mut self, tol: f64) -> Vec<usize> {
        const NONE: usize = usize::MAX;
        let mut replace: Vec<Option<(usize, bool)>> = vec![None; self.edges.len()];
        let mut gone = vec![false; self.edges.len()];
        let mut touched = vec![false; self.edges.len()];
        let mut dropped = Vec::new();
        let mut slits: Vec<(usize, usize)> = Vec::new();
        for (fi, f) in self.faces.iter().enumerate() {
            for li in 0..f.loops.len() {
                if f.loops[li].len() != 2 {
                    continue;
                }
                let whole = f.loops.len() == 1;
                let (a, b) = (f.loops[li][0].0, f.loops[li][1].0);
                if a == b || touched[a] || touched[b] {
                    continue;
                }
                let (ea, eb) = (&self.edges[a], &self.edges[b]);
                let same_ends = ea.v0 != ea.v1
                    && ((ea.v0, ea.v1) == (eb.v0, eb.v1) || (ea.v0, ea.v1) == (eb.v1, eb.v0));
                // Two straight edges between the same corners coincide; so
                // do two curves that run along each other, which is what a
                // face of no area between tangent surfaces is: where a
                // cylinder ends on the equator of a sphere of its own
                // radius, the mesh's cap pokes out between the sphere's
                // polygon and its own in slivers, and exactly the cap and
                // both edges of each sliver are the equator (a BOSL2
                // `stroke()` joint). Only whole faces of curves go, not
                // slits.
                let lines = matches!(ea.curve, Curve::Line { .. })
                    && matches!(eb.curve, Curve::Line { .. })
                    || (whole && same_ends && coincide(ea, eb, 10.0 * tol));
                let other = |e: &TEdge| {
                    if e.faces[0] == fi {
                        e.faces[1]
                    } else {
                        e.faces[0]
                    }
                };
                let (fa, fb) = (other(ea), other(eb));
                if whole && lines && same_ends && fa == fb && fa != NONE && fa != fi {
                    let g = &self.faces[fa];
                    let twin = g.loops.len() == 1 && g.loops[0].len() == 2 && {
                        let mut e2 = [g.loops[0][0].0, g.loops[0][1].0];
                        e2.sort_unstable();
                        let mut e1 = [a, b];
                        e1.sort_unstable();
                        e1 == e2
                    };
                    if twin {
                        gone[a] = true;
                        gone[b] = true;
                        touched[a] = true;
                        touched[b] = true;
                        dropped.push(fi);
                        dropped.push(fa);
                    }
                    continue;
                }
                // The faces either side must be two others, or the edge left
                // would bound one face twice.
                if !lines
                    || !same_ends
                    || fa == fb
                    || [fa, fb].iter().any(|&x| x == NONE || x == fi)
                {
                    continue;
                }
                // A slit whose sides are on one surface, facing opposite
                // ways, runs along a fin of no thickness: two flush walls
                // that the exact model cancels and the mesh's rounding
                // kept. Closing the slit would write a solid OCCT rejects
                // (`example017.scad`). A whole digon there is a sliver
                // between the two, which goes as any other
                // (`issue1165.scad`).
                let (ga, gb) = (&self.faces[fa], &self.faces[fb]);
                // The two sides as oriented surfaces: equal planes (a
                // faceted one and an exact one count) facing opposite ways.
                let opposite = match (&ga.surf, &gb.surf) {
                    (Surf::Plane { n: na, .. }, Surf::Plane { n: nb, .. }) => {
                        let sa = if ga.same_sense { 1.0 } else { -1.0 };
                        let sb = if gb.same_sense { 1.0 } else { -1.0 };
                        ga.surf.same(&gb.surf, tol) && na.dot(*nb) * sa * sb < 0.0
                    }
                    (a, b) => a == b && ga.same_sense != gb.same_sense,
                };
                if !whole && opposite {
                    continue;
                }
                replace[b] = Some((a, eb.v0 != ea.v0));
                touched[a] = true;
                touched[b] = true;
                if whole {
                    dropped.push(fi);
                } else {
                    slits.push((fi, li));
                }
            }
        }
        if dropped.is_empty() && slits.is_empty() {
            return dropped;
        }
        // Slit loops first, while the face numbers are still these.
        for &(fi, li) in slits.iter().rev() {
            self.faces[fi].loops.remove(li);
        }
        dropped.sort_unstable();
        dropped.dedup();
        let mut renum = vec![usize::MAX; self.edges.len()];
        let mut kept = Vec::with_capacity(self.edges.len());
        for (i, ed) in std::mem::take(&mut self.edges).into_iter().enumerate() {
            if replace[i].is_none() && !gone[i] {
                renum[i] = kept.len();
                kept.push(ed);
            }
        }
        self.edges = kept;
        let mut k = 0;
        self.faces.retain(|_| {
            let keep = dropped.binary_search(&k).is_err();
            k += 1;
            keep
        });
        for f in &mut self.faces {
            for lp in &mut f.loops {
                for c in lp.iter_mut() {
                    *c = match replace[c.0] {
                        Some((a, flip)) => (renum[a], c.1 != flip),
                        None => (renum[c.0], c.1),
                    };
                }
            }
        }
        // Which faces each edge separates, again (forward user first).
        for e in &mut self.edges {
            e.faces = [NONE, NONE];
        }
        for (fi, f) in self.faces.iter().enumerate() {
            for lp in &f.loops {
                for &(e, fwd) in lp {
                    self.edges[e].faces[usize::from(!fwd)] = fi;
                }
            }
        }
        dropped
    }

    /// Removes faces left with no loops by [`Topo::collapse_short_edges`]
    /// (a sliver whose every edge had no length), other than the lone
    /// spheres and tori that never had one. Returns the removed faces'
    /// indices, in order.
    pub fn drop_empty_faces(&mut self) -> Vec<usize> {
        const NONE: usize = usize::MAX;
        let dropped: Vec<usize> = (0..self.faces.len())
            .filter(|&f| {
                self.faces[f].loops.is_empty()
                    && !matches!(self.faces[f].surf, Surf::Sphere { .. } | Surf::Torus { .. })
            })
            .collect();
        if dropped.is_empty() {
            return dropped;
        }
        let mut k = 0;
        self.faces.retain(|_| {
            let keep = dropped.binary_search(&k).is_err();
            k += 1;
            keep
        });
        for e in &mut self.edges {
            e.faces = [NONE, NONE];
        }
        for (fi, f) in self.faces.iter().enumerate() {
            for lp in &f.loops {
                for &(e, fwd) in lp {
                    self.edges[e].faces[usize::from(!fwd)] = fi;
                }
            }
        }
        dropped
    }

    /// Moves the vertex of the closed conic edge `e` to the point at curve
    /// parameter `t`, keeping the curve.
    pub fn reseat_closed(&mut self, e: usize, t: f64) {
        let p = curve::eval(&self.edges[e].curve, t);
        let ed = &mut self.edges[e];
        match &mut ed.curve {
            Curve::Circle { center, x_axis, .. } => {
                let c = V::from(*center);
                let n = (p - c).norm();
                *x_axis = n.arr();
                ed.range = [0.0, TAU];
            }
            _ => ed.range = [t, t + TAU],
        }
        let v = ed.v0;
        self.verts[v] = p;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn line(t: &Topo, a: usize, b: usize, chain: Vec<V>) -> TEdge {
        let (p, q) = (t.verts[a], t.verts[b]);
        TEdge {
            v0: a,
            v1: b,
            curve: Curve::Line {
                origin: p.arr(),
                direction: (q - p).norm().arr(),
            },
            range: [0.0, (q - p).len()],
            faces: [usize::MAX; 2],
            chain,
            seam: false,
            dev: 0.0,
            chain_dev: 0.0,
        }
    }

    fn face(n: V, loops: Vec<Vec<(usize, bool)>>) -> TFace {
        TFace {
            surf: Surf::Plane {
                o: v(0.0, 0.0, 0.0),
                n,
            },
            same_sense: true,
            faceted: false,
            loops,
            param: None,
            pcurves: Vec::new(),
            outer: Vec::new(),
            tris: Vec::new(),
            source: Vec::new(),
        }
    }

    /// A wall (face 0) whose loop runs out along the line from corner 1
    /// to the tip, corner 2, and back past corner 1 to corner 0: a needle
    /// between face 1 (beyond the edge out) and face 2 (beyond the edge
    /// back), as a block's foot rounded in two passes left it. The needle
    /// goes, and faces 1 and 2 share its line.
    #[test]
    fn a_needle_between_two_faces_becomes_their_edge() {
        let mut t = Topo::default();
        for p in [
            v(0.0, 0.0, 0.0),
            v(10.0, 0.0, 0.0),
            v(12.0, 0.0, 0.0),
            v(10.0, 0.0, 5.0),
            v(0.0, 0.0, 5.0),
        ] {
            t.add_vertex(p);
        }
        // 0: 1 -> 2 (out), 1: 0 -> 2 (the long edge back), 2: 1 -> 3,
        // 3: 3 -> 4, 4: 4 -> 0.
        let es = [
            line(&t, 1, 2, Vec::new()),
            line(&t, 0, 2, Vec::new()),
            line(&t, 1, 3, Vec::new()),
            line(&t, 3, 4, Vec::new()),
            line(&t, 4, 0, Vec::new()),
        ];
        t.edges.extend(es);
        let y = v(0.0, 1.0, 0.0);
        t.faces = vec![
            face(
                y,
                vec![vec![
                    (2, false),
                    (0, true),
                    (1, false),
                    (4, false),
                    (3, false),
                ]],
            ),
            face(v(0.0, 0.0, 1.0), vec![vec![(0, false)]]),
            face(v(0.0, 0.0, -1.0), vec![vec![(1, true)]]),
        ];
        assert_eq!(t.remove_needles(1e-9), 1);
        // The wall keeps the long edge's part from corner 0 to corner 1.
        let end = |t: &Topo, (e, fwd): (usize, bool)| {
            if fwd { t.edges[e].v1 } else { t.edges[e].v0 }
        };
        let wall: Vec<(usize, usize)> = t.faces[0].loops[0]
            .iter()
            .map(|&c| (t.start(c), end(&t, c)))
            .collect();
        assert_eq!(wall, vec![(3, 1), (1, 0), (0, 4), (4, 3)]);
        // Faces 1 and 2 now share one edge from corner 1 to the tip.
        let e1 = t.faces[1].loops[0][0];
        let e2 = *t.faces[2].loops[0].last().expect("edge");
        assert_eq!(e1.0, e2.0);
        assert_ne!(e1.1, e2.1);
        let ed = &t.edges[e1.0];
        let mut ends = [ed.v0, ed.v1];
        ends.sort_unstable();
        assert_eq!(ends, [1, 2]);
        assert_eq!(t.edges.len(), 5);
    }

    /// A loop that only turns back on itself at an angle is no needle.
    #[test]
    fn a_corner_is_not_a_needle() {
        let mut t = Topo::default();
        for p in [v(0.0, 0.0, 0.0), v(10.0, 0.0, 0.0), v(0.0, 0.0, 5.0)] {
            t.add_vertex(p);
        }
        let es = [
            line(&t, 0, 1, Vec::new()),
            line(&t, 1, 2, Vec::new()),
            line(&t, 2, 0, Vec::new()),
        ];
        t.edges.extend(es);
        t.faces = vec![
            face(
                v(0.0, 1.0, 0.0),
                vec![vec![(0, true), (1, true), (2, true)]],
            ),
            face(
                v(0.0, -1.0, 0.0),
                vec![vec![(2, false), (1, false), (0, false)]],
            ),
        ];
        assert_eq!(t.remove_needles(1e-9), 0);
    }

    /// An edge whose two corners are one point and whose chain runs out
    /// and back along itself is a needle's; one whose chain goes round is
    /// not.
    #[test]
    fn a_retraced_chain_is_a_needle() {
        let mut t = Topo::default();
        for p in [v(0.0, 0.0, 0.0), v(0.0, 0.0, 1e-15)] {
            t.add_vertex(p);
        }
        let mut e = line(&t, 0, 1, vec![t.verts[0], v(0.4, 0.1, 0.0), t.verts[1]]);
        assert!(e.retraced(&t.verts, 1e-9));
        e.chain = vec![t.verts[0], v(0.4, 0.1, 0.0), v(0.5, 0.3, 0.0), t.verts[1]];
        assert!(!e.retraced(&t.verts, 1e-9));
    }

    /// An edge between two corners at one point whose fitted curve is a
    /// wisp while its mesh chain runs on is a contact the model does not
    /// have; one whose chain is as short, or whose corners are apart, is
    /// not.
    #[test]
    fn a_wisp_between_one_point_is_stray() {
        let mut t = Topo::default();
        for p in [v(0.0, 0.0, 0.0), v(1e-8, 0.0, 0.0)] {
            t.add_vertex(p);
        }
        let mut e = line(&t, 0, 1, vec![t.verts[0], v(0.2, 0.1, 0.0), t.verts[1]]);
        assert!(e.stray(&t.verts, 3e-7, 4e-8));
        e.chain = vec![t.verts[0], t.verts[1]];
        assert!(!e.stray(&t.verts, 3e-7, 4e-8));
        t.verts[1] = v(1e-3, 0.0, 0.0);
        e.chain = vec![t.verts[0], v(0.2, 0.1, 0.0), t.verts[1]];
        assert!(!e.stray(&t.verts, 3e-7, 4e-8));
    }
}
