#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// How hard a label box sticking out of the viewport is penalised, relative to leader length.
// Must match the same constant in polish.cpp so the discrete stage and the polish agree.
static const double OVERFLOW_WEIGHT = 10.0;

// Point-in-polygon by ray casting. Polygon is parallel x/y arrays, implicitly closed.
static inline bool pointInPolygon(double x, double y,
                                  const std::vector<double>& px, const std::vector<double>& py) {
  int m = (int) px.size();
  bool inside = false;
  for (int i = 0, j = m - 1; i < m; j = i++) {
    bool straddles = (py[i] > y) != (py[j] > y);
    if (straddles && x < (px[j] - px[i]) * (y - py[i]) / (py[j] - py[i]) + px[i])
      inside = !inside;
  }
  return inside;
}

// Length of the leader segment a->b that lies inside polygon (px, py). Split the segment at
// every crossing with a polygon edge; each sub-interval whose midpoint is inside contributes.
static double lengthInsidePolygon(double ax, double ay, double bx, double by,
                                  const std::vector<double>& px, const std::vector<double>& py) {
  double dx = bx - ax, dy = by - ay;
  double seglen = std::sqrt(dx * dx + dy * dy);
  int m = (int) px.size();

  // segment parameters t in [0,1] where the leader crosses a polygon edge, plus the endpoints
  std::vector<double> cuts;
  cuts.push_back(0.0);
  cuts.push_back(1.0);
  for (int e = 0, f = m - 1; e < m; f = e++) {
    double ox = px[f], oy = py[f], sx = px[e] - px[f], sy = py[e] - py[f];
    double den = dx * sy - dy * sx;
    if (std::fabs(den) < 1e-12) continue;                  // leader parallel to this edge
    double t = ((ox - ax) * sy - (oy - ay) * sx) / den;    // position along the leader
    double u = ((ox - ax) * dy - (oy - ay) * dx) / den;    // position along the edge
    if (t >= 0 && t <= 1 && u >= 0 && u <= 1) cuts.push_back(t);
  }
  std::sort(cuts.begin(), cuts.end());

  double inside = 0.0;
  for (std::size_t j = 0; j + 1 < cuts.size(); ++j) {
    double t0 = cuts[j], t1 = cuts[j + 1];
    if (t1 - t0 < 1e-12) continue;
    double tmid = 0.5 * (t0 + t1);
    if (pointInPolygon(ax + tmid * dx, ay + tmid * dy, px, py)) inside += (t1 - t0) * seglen;
  }
  return inside;
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

  // unpack polygons and cache their bounding boxes for a cheap broad-phase reject
  std::vector<std::vector<double> > polyX(K), polyY(K);
  std::vector<double> bxmin(K), bxmax(K), bymin(K), bymax(K);
  for (int k = 0; k < K; ++k) {
    NumericVector vx = polysx[k], vy = polysy[k];
    polyX[k].assign(vx.begin(), vx.end());
    polyY[k].assign(vy.begin(), vy.end());
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

    double ax = ex[i], ay = ey[i], bx = tx[i], by = ty[i];
    double sxmin = std::min(ax, bx), sxmax = std::max(ax, bx);
    double symin = std::min(ay, by), symax = std::max(ay, by);
    for (int k = 0; k < K; ++k) {
      if (lab[i] == k + 1) continue;                              // skip own cluster
      if (sxmin > bxmax[k] || sxmax < bxmin[k]                    // leader bbox misses cluster
          || symin > bymax[k] || symax < bymin[k]) continue;
      eff += lengthInsidePolygon(ax, ay, bx, by, polyX[k], polyY[k]);
    }
    out[i] = eff;
  }
  return out;
}
