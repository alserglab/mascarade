#include "geometry.h"    // firstRingHit
#include <vector>
using namespace Rcpp;

//' Visible leader end (first mask-boundary hit)
//'
//' For each label, the leader runs from its box anchor (`ex`, `ey`) to the cluster pole
//' (`tx`, `ty`), which lies inside the cluster. The drawn leader stops at the mask boundary:
//' the first point where the anchor->pole segment crosses the label's own cluster ring. When
//' the segment does not cross the ring the leader runs all the way to the pole.
//'
//' @param ex,ey Numeric leader-start (anchor) coordinates, one per label.
//' @param tx,ty Numeric pole (leader target) coordinates, one per label.
//' @param polysx,polysy Lists of parallel numeric x/y vectors, one ring per label (its own
//'   cluster).
//' @return A list with numeric `bx`, `by`: the visible leader end, one per label.
//' @keywords internal
// [[Rcpp::export]]
List firstLeaderHit(NumericVector ex, NumericVector ey,
                    NumericVector tx, NumericVector ty,
                    List polysx, List polysy) {
  int n = ex.size();
  NumericVector bx(n);
  NumericVector by(n);
  for (int i = 0; i < n; ++i) {
    NumericVector vx = polysx[i];
    NumericVector vy = polysy[i];
    double t = firstRingHit(ex[i], ey[i], tx[i], ty[i],
                            vx.begin(), vy.begin(), (int) vx.size());
    bx[i] = ex[i] + t * (tx[i] - ex[i]);
    by[i] = ey[i] + t * (ty[i] - ey[i]);
  }
  return List::create(_["bx"] = bx, _["by"] = by);
}
