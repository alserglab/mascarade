#ifndef MASCARADE_GEOMETRY_H
#define MASCARADE_GEOMETRY_H

// Shared geometry for the label placer.
//  - BoxFit: a Boost.Geometry R-tree over the cluster polygons, answering the "does a label
//    box overlap free space?" query exactly (replaces the prototype's occupancy grid).
//  - segcross / segbox: the per-pair predicates used by the refine and polish stages, thin
//    wrappers over Boost.Geometry's intersection tests (touching/collinear cases are treated
//    as conflicts -- the extra edge cases do not matter for the conflict scoring).

#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include <cmath>
#include <algorithm>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::d2::point_xy<double>     mpt;
typedef bg::model::polygon<mpt>             mpoly;
typedef bg::model::box<mpt>                 mbox;
typedef bg::model::segment<mpt>             mseg;
typedef std::pair<mbox, unsigned>           mval;

// Axis-aligned rectangle (a label box, bounding box, or the viewport). Bundles the four extents
// that were otherwise passed around as loose doubles, and centralises the overlap / overflow /
// gap tests.
struct Rect {
  double xmin, xmax, ymin, ymax;

  // bounding box of a segment (its two endpoints)
  static Rect ofSegment(double ax, double ay, double bx, double by) {
    return Rect{std::min(ax, bx), std::max(ax, bx), std::min(ay, by), std::max(ay, by)};
  }
  // do the two rectangles overlap? Closed (touching edges count); with plain doubles and no
  // exact-equality/epsilon handling the open/closed boundary case never affects the outcome.
  bool overlaps(const Rect& o) const {
    return xmin <= o.xmax && o.xmin <= xmax && ymin <= o.ymax && o.ymin <= ymax;
  }
  // total distance this box sticks out past `view` (0 when fully inside)
  double overflow(const Rect& view) const {
    return std::max(0.0, view.xmin - xmin) + std::max(0.0, xmax - view.xmax)
         + std::max(0.0, view.ymin - ymin) + std::max(0.0, ymax - view.ymax);
  }
  // Euclidean gap to `o` (0 when overlapping or touching)
  double gap(const Rect& o) const {
    double gx = std::max(std::max(xmin - o.xmax, o.xmin - xmax), 0.0);
    double gy = std::max(std::max(ymin - o.ymax, o.ymin - ymax), 0.0);
    return std::sqrt(gx * gx + gy * gy);
  }
};

// Even-odd ray-cast: is (px, py) inside the ring (vx, vy) of n vertices?
static inline bool pointInPoly(double px, double py,
                               const double* vx, const double* vy, int n) {
  bool inside = false;
  for (int i = 0, j = n - 1; i < n; j = i++) {
    if (((vy[i] > py) != (vy[j] > py))
        && (px < (vx[j] - vx[i]) * (py - vy[i]) / (vy[j] - vy[i]) + vx[i])) {
      inside = !inside;
    }
  }
  return inside;
}

// Length of segment (ax, ay)-(bx, by) that lies inside the simple polygon ring (vx, vy) of n
// vertices. Collects the parameters t where the segment crosses a polygon edge, then sums the
// sub-intervals that are inside (parity toggled from the inside-ness of the A endpoint). This
// is the single-segment case of a linestring/areal intersection length, hand-rolled to avoid
// Boost's general overlay machinery. `ts` is a caller-owned scratch buffer, reused per call.
static inline double segInsidePolyLen(double ax, double ay, double bx, double by,
                                      const double* vx, const double* vy, int n,
                                      std::vector<double>& ts) {
  double dx = bx - ax;
  double dy = by - ay;
  double segLen = std::sqrt(dx * dx + dy * dy);
  if (segLen <= 0.0) {
    return 0.0;
  }
  ts.clear();
  for (int i = 0, j = n - 1; i < n; j = i++) {
    double ex0 = vx[j];
    double ey0 = vy[j];
    double rx = vx[i] - ex0;
    double ry = vy[i] - ey0;
    double den = dx * ry - dy * rx;
    if (std::fabs(den) < 1e-12) {
      continue;                                          // segment parallel to this edge
    }
    double t = ((ex0 - ax) * ry - (ey0 - ay) * rx) / den;   // position along the segment
    double u = ((ex0 - ax) * dy - (ey0 - ay) * dx) / den;   // position along the edge
    if (t > 0.0 && t < 1.0 && u >= 0.0 && u <= 1.0) {
      ts.push_back(t);
    }
  }
  std::sort(ts.begin(), ts.end());
  bool inside = pointInPoly(ax, ay, vx, vy, n);            // inside-ness of the A end (t = 0)
  double prev = 0.0;
  double frac = 0.0;
  for (std::size_t c = 0; c < ts.size(); ++c) {
    if (inside) {
      frac += ts[c] - prev;
    }
    prev = ts[c];
    inside = !inside;
  }
  if (inside) {
    frac += 1.0 - prev;
  }
  return frac * segLen;
}

// First crossing of the directed segment (ax, ay)->(bx, by) with the simple polygon ring
// (vx, vy) of n vertices, returned as the fraction `t` in (0, 1] along the segment (so the hit
// point is A + t*(B - A)). Returns 1.0 when the segment does not cross the ring (the caller
// treats that as "runs to B"). Same crossing test as segInsidePolyLen(), kept hand-rolled for
// the same reason; here we only need the nearest crossing to A, so we track the minimum `t`.
static inline double firstRingHit(double ax, double ay, double bx, double by,
                                  const double* vx, const double* vy, int n) {
  double dx = bx - ax;
  double dy = by - ay;
  double best = 1.0;
  for (int i = 0, j = n - 1; i < n; j = i++) {
    double ex0 = vx[j];
    double ey0 = vy[j];
    double rx = vx[i] - ex0;
    double ry = vy[i] - ey0;
    double den = dx * ry - dy * rx;
    if (std::fabs(den) < 1e-12) {
      continue;                                          // segment parallel to this edge
    }
    double t = ((ex0 - ax) * ry - (ey0 - ay) * rx) / den;   // position along the segment
    double u = ((ex0 - ax) * dy - (ey0 - ay) * dx) / den;   // position along the edge
    if (t > 1e-6 && t <= 1.0 + 1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9 && t < best) {
      best = t;
    }
  }
  return best;
}

// Cached cluster mask rings (raw vertex arrays + bounding boxes) exposing the foreign-cluster
// arc query shared by effectiveLength() and the forcePolish() energy: the total length of a
// leader segment that runs inside any cluster other than its own. Built once per placement from
// the true (undilated) rings; a per-cluster bounding box rejects the common leader-misses-cluster
// case before the segment/polygon length test.
struct ClusterArcs {
  std::vector<std::vector<double> > cvx, cvy;         // per-cluster ring vertices
  std::vector<Rect> bbox;                             // per-cluster bounding box
  int K = 0;
  mutable std::vector<double> scratch;                // crossing-param scratch, reused per call

  void build(const Rcpp::List& polysx, const Rcpp::List& polysy) {
    K = polysx.size();
    cvx.resize(K);
    cvy.resize(K);
    bbox.resize(K);
    for (int k = 0; k < K; ++k) {
      Rcpp::NumericVector vx = polysx[k];
      Rcpp::NumericVector vy = polysy[k];
      cvx[k].assign(vx.begin(), vx.end());
      cvy[k].assign(vy.begin(), vy.end());
      bbox[k] = Rect{*std::min_element(vx.begin(), vx.end()),
                     *std::max_element(vx.begin(), vx.end()),
                     *std::min_element(vy.begin(), vy.end()),
                     *std::max_element(vy.begin(), vy.end())};
    }
  }

  // total length of segment (ax, ay)-(bx, by) inside any cluster other than `own` (0-based index)
  double foreignArc(int own, double ax, double ay, double bx, double by) const {
    Rect leader = Rect::ofSegment(ax, ay, bx, by);
    double arc = 0.0;
    for (int k = 0; k < K; ++k) {
      if (k == own || !leader.overlaps(bbox[k])) {
        continue;
      }
      arc += segInsidePolyLen(ax, ay, bx, by,
                              cvx[k].data(), cvy[k].data(), (int) cvx[k].size(), scratch);
    }
    return arc;
  }
};

// Box-fit acceleration structure: an R-tree of cluster-polygon envelopes plus the polygons.
// Built from the cluster polygons via buildBoxFit() and passed to the candidate / polish
// kernels as an XPtr handle for box-fit queries within a single placement.
struct BoxFit {
  std::vector<mpoly> polys;
  bgi::rtree<mval, bgi::quadratic<16> > tree;

  // does the axis-aligned rectangle `r` intersect ANY cluster polygon?
  bool hit(const Rect& r) const {
    mbox query(mpt(r.xmin, r.ymin), mpt(r.xmax, r.ymax));
    std::vector<mval> candidates;                          // envelopes overlapping the query
    tree.query(bgi::intersects(query), std::back_inserter(candidates));
    for (std::size_t i = 0; i < candidates.size(); ++i) {
      if (bg::intersects(query, polys[candidates[i].second])) {
        return true;
      }
    }
    return false;
  }
};

// Segment-segment intersection of a->b and c->d (Boost.Geometry; touching/collinear count).
static inline bool segcross(double ax, double ay, double bx, double by,
                            double cx, double cy, double dx, double dy) {
  mseg s1(mpt(ax, ay), mpt(bx, by));
  mseg s2(mpt(cx, cy), mpt(dx, dy));
  return bg::intersects(s1, s2);
}

// Segment (ax, ay)-(bx, by) versus the axis-aligned rectangle `r` (Boost.Geometry): a hit if the
// segment intersects the box (interior or boundary).
static inline bool segbox(double ax, double ay, double bx, double by, const Rect& r) {
  mseg s(mpt(ax, ay), mpt(bx, by));
  mbox box(mpt(r.xmin, r.ymin), mpt(r.xmax, r.ymax));
  return bg::intersects(s, box);
}

#endif
