#include "effective.h"    // effectiveLengthImpl + ClusterArcs (the shared kernel)
#include <vector>
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
  int n = ex.size();
  ClusterArcs arcs;
  arcs.build(polysx, polysy);

  Rect view{xlo, xhi, ylo, yhi};
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    Rect box{cxmin[i], cxmax[i], cymin[i], cymax[i]};
    out[i] = effectiveLengthImpl(arcs, lab[i] - 1,
                                 ex[i], ey[i], tx[i], ty[i], len[i], box, view);
  }
  return out;
}
