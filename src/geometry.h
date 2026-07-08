#ifndef MASCARADE_GEOMETRY_H
#define MASCARADE_GEOMETRY_H

// Shared geometry for the label placer.
//  - BoxFit: a Boost.Geometry R-tree over the cluster polygons, answering the
//    "does a label box overlap free space?" query exactly (replaces the prototype's
//    occupancy grid + integral image).
//  - segcross / segbox: the per-pair predicates used by the refine and polish stages,
//    kept hand-rolled so their conflict logic stays identical to the validated prototype.

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

// Box-fit acceleration structure: R-tree of cluster-polygon envelopes + the polygons.
// Built from the cluster polygons via buildBoxFit() and passed to the candidate/polish
// kernels as an XPtr handle for box-fit queries within a single placement.
struct BoxFit {
  std::vector<mpoly> polys;
  bgi::rtree<mval, bgi::quadratic<16> > tree;

  // does the axis-aligned box [xmin,xmax] x [ymin,ymax] intersect ANY cluster polygon?
  bool hit(double xmin, double xmax, double ymin, double ymax) const {
    mbox qb(mpt(xmin, ymin), mpt(xmax, ymax));
    std::vector<mval> hits;
    tree.query(bgi::intersects(qb), std::back_inserter(hits));
    for (std::size_t i = 0; i < hits.size(); ++i)
      if (bg::intersects(qb, polys[hits[i].second])) return true;
    return false;
  }
};

// proper segment-segment crossing (interior only; touching endpoints do not count)
static inline bool segcross(double ax, double ay, double bx, double by,
                            double cx, double cy, double dx, double dy) {
  double d1x = bx - ax, d1y = by - ay, d2x = dx - cx, d2y = dy - cy;
  double den = d1x * d2y - d1y * d2x;
  if (std::fabs(den) <= 1e-12) return false;
  double t = ((cx - ax) * d2y - (cy - ay) * d2x) / den;
  double u = ((cx - ax) * d1y - (cy - ay) * d1x) / den;
  return t > 1e-9 && t < 1 - 1e-9 && u > 1e-9 && u < 1 - 1e-9;
}

// segment (ax,ay)-(bx,by) vs axis-aligned box [x0,x1] x [y0,y1] (endpoints inside count)
static inline bool segbox(double ax, double ay, double bx, double by,
                          double x0, double x1, double y0, double y1) {
  if (ax > x0 && ax < x1 && ay > y0 && ay < y1) return true;
  if (bx > x0 && bx < x1 && by > y0 && by < y1) return true;
  if (segcross(ax, ay, bx, by, x0, y0, x1, y0)) return true;
  if (segcross(ax, ay, bx, by, x1, y0, x1, y1)) return true;
  if (segcross(ax, ay, bx, by, x1, y1, x0, y1)) return true;
  if (segcross(ax, ay, bx, by, x0, y1, x0, y0)) return true;
  return false;
}

#endif
