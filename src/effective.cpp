#include "geometry.h"     // Rcpp + Boost.Geometry setup (bg, mpt, mpoly)
#include "constants.h"    // OVERFLOW_WEIGHT
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <vector>
#include <algorithm>
using namespace Rcpp;

typedef bg::model::linestring<mpt>        mline;
typedef bg::model::multi_linestring<mline> mmline;

// Build a Boost polygon (corrected: closed ring, consistent orientation) from parallel x/y.
static mpoly makePolygon(const NumericVector& vx, const NumericVector& vy) {
  mpoly poly;
  for (int e = 0; e < vx.size(); ++e) {
    bg::append(poly.outer(), mpt(vx[e], vy[e]));
  }
  bg::correct(poly);
  return poly;
}

// Effective length the optimizer minimises (after the conflict counts, lexicographically),
// per candidate:
//       base leader length  (`len`)
//     + arc of the leader inside any FOREIGN cluster        -> routes leaders around clusters
//     + OVERFLOW_WEIGHT * box overflow beyond the viewport  -> steers labels in-bounds
// `lab` is the candidate's own 1-indexed cluster; cx/cy min/max are the padded label box;
// xlo/xhi/ylo/yhi are the (already label.buffer-inset) viewport.
// [[Rcpp::export]]
NumericVector effectiveLength(NumericVector len, NumericVector ex, NumericVector ey,
                              NumericVector tx, NumericVector ty, IntegerVector lab,
                              List polysx, List polysy,
                              NumericVector cxmin, NumericVector cxmax,
                              NumericVector cymin, NumericVector cymax,
                              double xlo, double xhi, double ylo, double yhi) {
  int n = ex.size(), K = polysx.size();

  // build each cluster polygon once + cache its bounding box for a cheap broad-phase reject
  std::vector<mpoly> polys(K);
  std::vector<double> bxmin(K), bxmax(K), bymin(K), bymax(K);
  for (int k = 0; k < K; ++k) {
    NumericVector vx = polysx[k], vy = polysy[k];
    polys[k] = makePolygon(vx, vy);
    bxmin[k] = *std::min_element(vx.begin(), vx.end());
    bxmax[k] = *std::max_element(vx.begin(), vx.end());
    bymin[k] = *std::min_element(vy.begin(), vy.end());
    bymax[k] = *std::max_element(vy.begin(), vy.end());
  }

  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double overflow = std::max(0.0, xlo - cxmin[i]) + std::max(0.0, cxmax[i] - xhi)
                    + std::max(0.0, ylo - cymin[i]) + std::max(0.0, cymax[i] - yhi);
    double eff = len[i] + OVERFLOW_WEIGHT * overflow;

    mline leader;
    leader.push_back(mpt(ex[i], ey[i]));
    leader.push_back(mpt(tx[i], ty[i]));
    double sxmin = std::min(ex[i], tx[i]), sxmax = std::max(ex[i], tx[i]);
    double symin = std::min(ey[i], ty[i]), symax = std::max(ey[i], ty[i]);

    for (int k = 0; k < K; ++k) {
      if (lab[i] == k + 1) {
        continue;                                                  // skip own cluster
      }
      if (sxmin > bxmax[k] || sxmax < bxmin[k] || symin > bymax[k] || symax < bymin[k]) {
        continue;                                                  // leader bbox misses cluster
      }
      mmline inside;
      bg::intersection(leader, polys[k], inside);                  // parts of the leader inside k
      eff += bg::length(inside);
    }
    out[i] = eff;
  }
  return out;
}
