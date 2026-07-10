#include "geometry.h"     // Rcpp + Boost.Geometry setup (bg, mpt, mpoly, mline, makePolygon)
#include "constants.h"    // OVERFLOW_WEIGHT
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' Effective length used to rank label candidates
//'
//' Per candidate, the length the optimizer minimises after the conflict counts: the base
//' leader length, plus the arc of the leader inside any FOREIGN cluster (routes leaders around
//' clusters), plus `OVERFLOW_WEIGHT` times how far the box overflows the viewport (steers
//' labels in-bounds).
//'
//' @param len Numeric base leader length per candidate.
//' @param ex,ey Numeric leader start (anchor) per candidate.
//' @param tx,ty Numeric leader target (pole) per candidate.
//' @param lab Integer 1-indexed own-cluster of each candidate.
//' @param polysx,polysy Lists of parallel numeric x/y vectors, one ring per cluster.
//' @param cxmin,cxmax,cymin,cymax Numeric padded box extents per candidate.
//' @param xlo,xhi,ylo,yhi Numeric viewport bounds (already inset by label.buffer).
//' @return Numeric vector of effective lengths, one per candidate.
//' @keywords internal
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
