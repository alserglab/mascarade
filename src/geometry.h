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

// Box-fit acceleration structure: an R-tree of cluster-polygon envelopes plus the polygons.
// Built from the cluster polygons via buildBoxFit() and passed to the candidate / polish
// kernels as an XPtr handle for box-fit queries within a single placement.
struct BoxFit {
  std::vector<mpoly> polys;
  bgi::rtree<mval, bgi::quadratic<16> > tree;

  // does the axis-aligned box [xmin, xmax] x [ymin, ymax] intersect ANY cluster polygon?
  bool hit(double xmin, double xmax, double ymin, double ymax) const {
    mbox query(mpt(xmin, ymin), mpt(xmax, ymax));
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

// Segment (ax, ay)-(bx, by) versus the axis-aligned box [x0, x1] x [y0, y1] (Boost.Geometry):
// a hit if the segment intersects the box (interior or boundary).
static inline bool segbox(double ax, double ay, double bx, double by,
                          double x0, double x1, double y0, double y1) {
  mseg s(mpt(ax, ay), mpt(bx, by));
  mbox box(mpt(x0, y0), mpt(x1, y1));
  return bg::intersects(s, box);
}

#endif
