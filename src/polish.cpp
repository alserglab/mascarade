#include "effective.h"   // effectiveLengthImpl + ClusterArcs; pulls in geometry.h
#include <vector>
#include <algorithm>
using namespace Rcpp;

//' Continuous force-directed label polish
//'
//' Pattern-search descent on the (squared) EFFECTIVE length under a hard conflict guard. The
//' effective length of a label is its centre-to-pole leader length, plus the arc of the leader
//' that runs inside any foreign cluster (routing leaders around clusters), plus a soft viewport
//' overflow penalty -- the same quantity minimised by the upstream `effectiveLength()` ranking,
//' so the continuous polish and the discrete candidate refinement optimise a consistent
//' objective. A box-box spacing penalty is added on top. Starting from a conflict-free layout
//' the search only accepts conflict-free neighbours (free-space check via the BoxFit R-tree),
//' so feasibility is preserved. Overflow is a SOFT term folded into the effective length, not a
//' hard clip, so an off-panel label can be walked back in-bounds and a label leaves the panel
//' only when that lowers total energy. The leader anchor rule mirrors R's `.anchorPoint()`, so
//' both the conflict guard and the foreign-cluster arc match the drawn geometry.
//'
//' @param boxfit External pointer from `buildBoxFit()` (the cluster keep-out).
//' @param cx0,cy0 Numeric starting label-centre coordinates.
//' @param hw,hh Numeric per-label box half-sizes.
//' @param tx,ty Numeric per-label pole (leader target).
//' @param polysx,polysy Lists of parallel numeric x/y vectors, one true mask ring per cluster
//'   (aligned with the labels), used for the foreign-cluster arc term.
//' @param pad Numeric hard box clearance.
//' @param xlo,xhi,ylo,yhi Numeric viewport bounds.
//' @param iters Integer iteration count.
//' @param step Numeric initial pattern-search step.
//' @param MU Numeric weight of the box-box spacing penalty.
//' @param pad_tgt Numeric target inter-box spacing.
//' @param stepmin Numeric smallest step tried before abandoning a direction.
//' @param con_type Integer leader style: 0 = "cl" (corner), otherwise "cm"/"none".
//' @param sq Logical; if `TRUE`, the length term uses the squared distance.
//' @return A list with numeric `cx`, `cy`: the polished label centres.
//' @keywords internal
// [[Rcpp::export]]
List forcePolish(SEXP boxfit, NumericVector cx0, NumericVector cy0,
                 NumericVector hw, NumericVector hh, NumericVector tx, NumericVector ty,
                 List polysx, List polysy,
                 double pad, double xlo, double xhi, double ylo, double yhi,
                 int iters, double step, double MU, double pad_tgt, double stepmin,
                 int con_type, bool sq = true) {
  XPtr<BoxFit> boxFit(boxfit);
  int n = cx0.size();
  std::vector<double> cx(cx0.begin(), cx0.end());
  std::vector<double> cy(cy0.begin(), cy0.end());
  // padded box + leader anchor cached for each label's CURRENT placement; the viewport rectangle
  std::vector<Rect> box(n);
  std::vector<double> anchorX(n), anchorY(n);
  Rect view{xlo, xhi, ylo, yhi};

  // shared foreign-cluster arc field over the true cluster rings (same computation the upstream
  // effectiveLength() uses). Label i's own cluster is i.
  ClusterArcs arcs;
  arcs.build(polysx, polysy);

  // Leader start on label i's box when its centre is (X, Y), aimed at the pole. Mirrors R's
  // .anchorPoint(): con_type 0 = "cl" sign-quadrant corner, else the "cm"/"none" 8-point rule.
  auto updateLeaderAnchor = [&](int i, double X, double Y, double& ax, double& ay) {
    double dx = tx[i] - X;
    double dy = ty[i] - Y;
    double sX = dx >= 0 ? 1.0 : -1.0;
    double sY = dy >= 0 ? 1.0 : -1.0;
    if (con_type == 0) {
      ax = X + sX * hw[i];
      ay = Y + sY * hh[i];
    } else {
      bool xin = std::fabs(dx) < hw[i];
      bool yin = std::fabs(dy) < hh[i];
      ax = (xin && !yin) ? X : X + sX * hw[i];
      ay = (!xin && yin) ? Y : Y + sY * hh[i];
    }
  };
  // padded box of label i centred at (X, Y)
  auto paddedBox = [&](int i, double X, double Y) -> Rect {
    return Rect{X - hw[i] - pad, X + hw[i] + pad, Y - hh[i] - pad, Y + hh[i] + pad};
  };
  // cache label i's padded box + leader anchor for centre (X, Y)
  auto updateGeom = [&](int i, double X, double Y) {
    box[i] = paddedBox(i, X, Y);
    updateLeaderAnchor(i, X, Y, anchorX[i], anchorY[i]);
  };
  for (int i = 0; i < n; ++i) {
    updateGeom(i, cx[i], cy[i]);
  }

  // is centre (X, Y) for label i conflict-free? (the viewport is a soft energy term, not here)
  auto conflictFree = [&](int i, double X, double Y) -> bool {
    // label i's box (no pad) must not overlap any cluster polygon
    if (boxFit->hit(Rect{X - hw[i], X + hw[i], Y - hh[i], Y + hh[i]})) {
      return false;
    }
    Rect me = paddedBox(i, X, Y);
    double ax, ay;
    updateLeaderAnchor(i, X, Y, ax, ay);
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      // box-box overlap
      if (me.overlaps(box[j])) {
        return false;
      }
      // my leader through box j, or leader j through my box
      if (segbox(ax, ay, tx[i], ty[i], box[j])) {
        return false;
      }
      if (segbox(anchorX[j], anchorY[j], tx[j], ty[j], me)) {
        return false;
      }
      // leaders crossing (always a hard constraint)
      if (segcross(ax, ay, tx[i], ty[i], anchorX[j], anchorY[j], tx[j], ty[j])) {
        return false;
      }
    }
    return true;
  };
  // energy of centre (X, Y) for label i: the (squared) EFFECTIVE length -- leader length to the
  // pole, plus the leader's arc inside foreign clusters, plus a SOFT viewport-overflow penalty
  // folded in (no longer a separate additive term) -- then the box-box spacing penalty on top.
  auto energy = [&](int i, double X, double Y) -> double {
    Rect me = paddedBox(i, X, Y);
    double dist = std::sqrt((X - tx[i]) * (X - tx[i]) + (Y - ty[i]) * (Y - ty[i]));
    double ax, ay;
    updateLeaderAnchor(i, X, Y, ax, ay);
    double effLen = effectiveLengthImpl(arcs, i, ax, ay, tx[i], ty[i], dist, me, view);
    double e = sq ? effLen * effLen : effLen;
    for (int j = 0; j < n; ++j) {
      if (j == i) {
        continue;
      }
      double gap = me.gap(box[j]);
      if (gap < pad_tgt) {
        double deficit = pad_tgt - gap;
        e += MU * deficit * deficit;
      }
    }
    return e;
  };

  // search directions: nDir evenly around the circle, plus one straight at the pole (index nDir)
  const int nDir = 16;
  double dirX[nDir + 1];
  double dirY[nDir + 1];
  for (int k = 0; k < nDir; ++k) {
    double theta = 2 * M_PI * k / nDir;
    dirX[k] = std::cos(theta);
    dirY[k] = std::sin(theta);
  }

  for (int iter = 0; iter < iters; ++iter) {
    bool moved = false;
    for (int i = 0; i < n; ++i) {
      double X = cx[i];
      double Y = cy[i];
      double toPole = std::sqrt((tx[i] - X) * (tx[i] - X) + (ty[i] - Y) * (ty[i] - Y));
      dirX[nDir] = toPole > 1e-9 ? (tx[i] - X) / toPole : 0;
      dirY[nDir] = toPole > 1e-9 ? (ty[i] - Y) / toPole : 0;
      double bestEnergy = energy(i, X, Y);
      double bestX = X;
      double bestY = Y;
      for (int k = 0; k <= nDir; ++k) {
        // halve the step until a conflict-free neighbour is found in this direction
        double stepSize = step;
        while (stepSize > stepmin) {
          double nx = X + stepSize * dirX[k];
          double ny = Y + stepSize * dirY[k];
          if (conflictFree(i, nx, ny)) {
            double candEnergy = energy(i, nx, ny);
            if (candEnergy < bestEnergy - 1e-9) {
              bestEnergy = candEnergy;
              bestX = nx;
              bestY = ny;
            }
            break;
          }
          stepSize *= 0.5;
        }
      }
      if (bestX != X || bestY != Y) {
        cx[i] = bestX;
        cy[i] = bestY;
        updateGeom(i, bestX, bestY);
        moved = true;
      }
    }
    if (!moved) {
      break;
    }
  }

  return List::create(_["cx"] = NumericVector(cx.begin(), cx.end()),
                      _["cy"] = NumericVector(cy.begin(), cy.end()));
}
