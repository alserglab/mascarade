#ifndef MASCARADE_GEOMETRY_H
#define MASCARADE_GEOMETRY_H

// Shared geometry for the label placer.
//  - BoxFit: a Boost.Geometry R-tree over the cluster polygons, answering the "does a label
//    box overlap free space?" query exactly (replaces the prototype's occupancy grid).
//  - segcross / segbox: the per-pair predicates used by the refine and polish stages, kept
//    hand-rolled so their conflict logic stays identical to the validated prototype.

#include <Rcpp.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include <cmath>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::d2::point_xy<double>  mpt;
typedef bg::model::polygon<mpt>          mpoly;
typedef bg::model::box<mpt>              mbox;
typedef std::pair<mbox, unsigned>        mval;

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

// Proper segment-segment crossing of a->b and c->d (interior only; touching endpoints do not
// count).
static inline bool segcross(double ax, double ay, double bx, double by,
                            double cx, double cy, double dx, double dy) {
  double abx = bx - ax;
  double aby = by - ay;
  double cdx = dx - cx;
  double cdy = dy - cy;
  double den = abx * cdy - aby * cdx;
  if (std::fabs(den) <= 1e-12) {
    return false;
  }
  double t = ((cx - ax) * cdy - (cy - ay) * cdx) / den;   // position along a->b
  double u = ((cx - ax) * aby - (cy - ay) * abx) / den;   // position along c->d
  return t > 1e-9 && t < 1 - 1e-9 && u > 1e-9 && u < 1 - 1e-9;
}

// Segment (ax, ay)-(bx, by) versus the axis-aligned box [x0, x1] x [y0, y1]: a hit if either
// endpoint is strictly inside the box or the segment crosses any box edge.
static inline bool segbox(double ax, double ay, double bx, double by,
                          double x0, double x1, double y0, double y1) {
  if (ax > x0 && ax < x1 && ay > y0 && ay < y1) {
    return true;
  }
  if (bx > x0 && bx < x1 && by > y0 && by < y1) {
    return true;
  }
  if (segcross(ax, ay, bx, by, x0, y0, x1, y0)) {
    return true;
  }
  if (segcross(ax, ay, bx, by, x1, y0, x1, y1)) {
    return true;
  }
  if (segcross(ax, ay, bx, by, x1, y1, x0, y1)) {
    return true;
  }
  if (segcross(ax, ay, bx, by, x0, y1, x0, y0)) {
    return true;
  }
  return false;
}

#endif
