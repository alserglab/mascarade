#include "geometry.h"     // Rcpp + Boost.Geometry setup; segInsidePolyLen for the arc term
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

  // cache each cluster ring's vertices (raw arrays) + its bounding box for a broad-phase reject
  std::vector<std::vector<double> > cvx(K), cvy(K);
  std::vector<double> bxmin(K), bxmax(K), bymin(K), bymax(K);
  for (int k = 0; k < K; ++k) {
    NumericVector vx = polysx[k], vy = polysy[k];
    cvx[k].assign(vx.begin(), vx.end());
    cvy[k].assign(vy.begin(), vy.end());
    bxmin[k] = *std::min_element(vx.begin(), vx.end());
    bxmax[k] = *std::max_element(vx.begin(), vx.end());
    bymin[k] = *std::min_element(vy.begin(), vy.end());
    bymax[k] = *std::max_element(vy.begin(), vy.end());
  }

  NumericVector out(n);
  std::vector<double> ts;                                          // crossing-param scratch
  for (int i = 0; i < n; ++i) {
    double overflow = std::max(0.0, xlo - cxmin[i]) + std::max(0.0, cxmax[i] - xhi)
                    + std::max(0.0, ylo - cymin[i]) + std::max(0.0, cymax[i] - yhi);
    double eff = len[i] + OVERFLOW_WEIGHT * overflow;

    double ax = ex[i], ay = ey[i], bx = tx[i], by = ty[i];
    double sxmin = std::min(ax, bx), sxmax = std::max(ax, bx);
    double symin = std::min(ay, by), symax = std::max(ay, by);

    for (int k = 0; k < K; ++k) {
      if (lab[i] == k + 1) {
        continue;                                                  // skip own cluster
      }
      if (sxmin > bxmax[k] || sxmax < bxmin[k] || symin > bymax[k] || symax < bymin[k]) {
        continue;                                                  // leader bbox misses cluster
      }
      eff += segInsidePolyLen(ax, ay, bx, by,
                              cvx[k].data(), cvy[k].data(), (int) cvx[k].size(), ts);
    }
    out[i] = eff;
  }
  return out;
}
