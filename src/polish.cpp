#include "geometry.h"
#include <vector>
using namespace Rcpp;

// Force-directed polish: one energy (squared leader length + box-box padding) under a hard
// conflict guard, pattern-search descent. Free-space check (label box vs clusters) via the
// BoxFit R-tree, replacing the prototype's occupancy integral image. Feasibility preserved.
// [[Rcpp::export]]
List forcePolish(SEXP boxfit, NumericVector cx0, NumericVector cy0, NumericVector hw, NumericVector hh,
                 NumericVector tx, NumericVector ty, double pad,
                 double xlo, double xhi, double ylo, double yhi,
                 int iters, double step, double MU, double pad_tgt, double stepmin,
                 bool ll_hard = true, bool sq = true) {
  XPtr<BoxFit> bf(boxfit);
  int n = cx0.size();
  std::vector<double> cx(cx0.begin(), cx0.end()), cy(cy0.begin(), cy0.end());
  std::vector<double> CXm(n), CXM(n), CYm(n), CYM(n), EX(n), EY(n);
  auto lead = [&](int i, double X, double Y, double& ex, double& ey) {
    double dx = tx[i] - X, dy = ty[i] - Y, adx = std::fabs(dx), ady = std::fabs(dy);
    double tt = 1e18; if (adx > 1e-9) tt = std::min(tt, hw[i] / adx); if (ady > 1e-9) tt = std::min(tt, hh[i] / ady);
    if (tt > 1e17) tt = 0; ex = X + dx * tt; ey = Y + dy * tt; };
  auto setg = [&](int i, double X, double Y) { CXm[i] = X - hw[i] - pad; CXM[i] = X + hw[i] + pad; CYm[i] = Y - hh[i] - pad; CYM[i] = Y + hh[i] + pad;
    lead(i, X, Y, EX[i], EY[i]); };
  for (int i = 0; i < n; ++i) setg(i, cx[i], cy[i]);
  // label box (no pad) must not overlap any cluster polygon
  auto forb = [&](int i, double X, double Y)->bool { return bf->hit(X - hw[i], X + hw[i], Y - hh[i], Y + hh[i]); };
  auto valid = [&](int i, double X, double Y)->bool {
    if (X - hw[i] < xlo || X + hw[i] > xhi || Y - hh[i] < ylo || Y + hh[i] > yhi) return false;
    if (forb(i, X, Y)) return false;
    double bxm = X - hw[i] - pad, bxM = X + hw[i] + pad, bym = Y - hh[i] - pad, byM = Y + hh[i] + pad, ex, ey; lead(i, X, Y, ex, ey);
    for (int j = 0; j < n; ++j) { if (j == i) continue;
      if (bxm < CXM[j] && CXm[j] < bxM && bym < CYM[j] && CYm[j] < byM) return false;
      if (segbox(ex, ey, tx[i], ty[i], CXm[j], CXM[j], CYm[j], CYM[j])) return false;
      if (segbox(EX[j], EY[j], tx[j], ty[j], bxm, bxM, bym, byM)) return false;
      if (ll_hard && segcross(ex, ey, tx[i], ty[i], EX[j], EY[j], tx[j], ty[j])) return false; }
    return true; };
  auto energy = [&](int i, double X, double Y)->double {
    double bxm = X - hw[i] - pad, bxM = X + hw[i] + pad, bym = Y - hh[i] - pad, byM = Y + hh[i] + pad;
    double d2 = (X - tx[i]) * (X - tx[i]) + (Y - ty[i]) * (Y - ty[i]);
    double e = sq ? d2 : std::sqrt(d2);
    for (int j = 0; j < n; ++j) { if (j == i) continue;
      double gx = std::max(std::max(bxm - CXM[j], CXm[j] - bxM), 0.0), gy = std::max(std::max(bym - CYM[j], CYm[j] - byM), 0.0);
      double gp = std::sqrt(gx * gx + gy * gy); if (gp < pad_tgt) { double d = pad_tgt - gp; e += MU * d * d; } }
    return e; };
  const int ND = 16; double dvx[ND + 1], dvy[ND + 1];
  for (int k = 0; k < ND; ++k) { double th = 2 * M_PI * k / ND; dvx[k] = std::cos(th); dvy[k] = std::sin(th); }
  for (int it = 0; it < iters; ++it) { bool moved = false;
    for (int i = 0; i < n; ++i) { double X = cx[i], Y = cy[i];
      double tl = std::sqrt((tx[i] - X) * (tx[i] - X) + (ty[i] - Y) * (ty[i] - Y));
      dvx[ND] = tl > 1e-9 ? (tx[i] - X) / tl : 0; dvy[ND] = tl > 1e-9 ? (ty[i] - Y) / tl : 0;
      double bestE = energy(i, X, Y), bX = X, bY = Y;
      for (int k = 0; k <= ND; ++k) { double t = step;
        while (t > stepmin) { double nx2 = X + t * dvx[k], ny2 = Y + t * dvy[k];
          if (valid(i, nx2, ny2)) { double E = energy(i, nx2, ny2); if (E < bestE - 1e-9) { bestE = E; bX = nx2; bY = ny2; } break; } t *= 0.5; } }
      if (bX != X || bY != Y) { cx[i] = bX; cy[i] = bY; setg(i, bX, bY); moved = true; } }
    if (!moved) break; }
  return List::create(_["cx"] = NumericVector(cx.begin(), cx.end()),
                      _["cy"] = NumericVector(cy.begin(), cy.end()));
}
