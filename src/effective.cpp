#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// point-in-polygon (ray casting), polygon given as parallel x/y (implicitly closed).
static inline bool pip(double x, double y, const std::vector<double>& X, const std::vector<double>& Y) {
  int m = (int) X.size(); bool in = false;
  for (int i = 0, j = m - 1; i < m; j = i++)
    if (((Y[i] > y) != (Y[j] > y)) && (x < (X[j] - X[i]) * (y - Y[i]) / (Y[j] - Y[i]) + X[i])) in = !in;
  return in;
}

// Effective length ranked by the optimizer (after conflicts, lexicographically). Per candidate:
//   base leader length  +  arc of the leader inside any FOREIGN cluster (routes leaders around
//   clusters)  +  how far the label box overflows the viewport on x and y (steers labels
//   in-bounds). lab = candidate's own 1-indexed cluster; cx/cy min/max are the padded box.
// [[Rcpp::export]]
NumericVector effectiveLength(NumericVector len, NumericVector ex, NumericVector ey,
                              NumericVector tx, NumericVector ty, IntegerVector lab,
                              List polysx, List polysy,
                              NumericVector cxmin, NumericVector cxmax,
                              NumericVector cymin, NumericVector cymax,
                              double xlo, double xhi, double ylo, double yhi) {
  int n = ex.size(), K = polysx.size();
  std::vector<std::vector<double> > PX(K), PY(K);
  std::vector<double> bxmin(K), bxmax(K), bymin(K), bymax(K);
  for (int k = 0; k < K; ++k) { NumericVector vx = polysx[k], vy = polysy[k]; int m = vx.size();
    PX[k].assign(vx.begin(), vx.end()); PY[k].assign(vy.begin(), vy.end());
    double xm = 1e300, xM = -1e300, ym = 1e300, yM = -1e300;
    for (int e = 0; e < m; ++e) { xm = std::min(xm, vx[e]); xM = std::max(xM, vx[e]); ym = std::min(ym, vy[e]); yM = std::max(yM, vy[e]); }
    bxmin[k] = xm; bxmax[k] = xM; bymin[k] = ym; bymax[k] = yM; }
  NumericVector out(n);
  for (int i = 0; i < n; ++i) { double ax = ex[i], ay = ey[i], bx = tx[i], by = ty[i];
    // overflow: box extent outside the viewport, summed over both axes
    double over = std::max(0.0, xlo - cxmin[i]) + std::max(0.0, cxmax[i] - xhi)
                + std::max(0.0, ylo - cymin[i]) + std::max(0.0, cymax[i] - yhi);
    double rx = bx - ax, ry = by - ay, seglen = std::sqrt(rx * rx + ry * ry);
    if (seglen <= 0) { out[i] = len[i] + over; continue; }
    double sxmin = std::min(ax, bx), sxmax = std::max(ax, bx), symin = std::min(ay, by), symax = std::max(ay, by);
    double acc = 0;
    for (int k = 0; k < K; ++k) { if (lab[i] == k + 1) continue;
      if (sxmin > bxmax[k] || sxmax < bxmin[k] || symin > bymax[k] || symax < bymin[k]) continue;
      const std::vector<double>& X = PX[k]; const std::vector<double>& Y = PY[k]; int m = (int) X.size();
      std::vector<double> ts; ts.push_back(0.0); ts.push_back(1.0);
      for (int e = 0, f = m - 1; e < m; f = e++) { double p0x = X[f], p0y = Y[f], sx = X[e] - X[f], sy = Y[e] - Y[f];
        double den = rx * sy - ry * sx; if (std::fabs(den) < 1e-12) continue;
        double t = ((p0x - ax) * sy - (p0y - ay) * sx) / den, u = ((p0x - ax) * ry - (p0y - ay) * rx) / den;
        if (t >= 0 && t <= 1 && u >= 0 && u <= 1) ts.push_back(t); }
      std::sort(ts.begin(), ts.end());
      for (std::size_t j = 0; j + 1 < ts.size(); ++j) { double t0 = ts[j], t1 = ts[j + 1]; if (t1 - t0 < 1e-12) continue;
        double tmid = 0.5 * (t0 + t1);
        if (pip(ax + tmid * rx, ay + tmid * ry, X, Y)) acc += (t1 - t0) * seglen; } }
    out[i] = len[i] + acc + over; }
  return out;
}
