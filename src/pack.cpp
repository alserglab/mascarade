#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

//' Fixed-order 1-D packing minimising total leader length
//'
//' For labels in a fixed top-to-bottom order, chooses stacked centre y-positions on a fine grid
//' that minimise the total leader length `sqrt(dx^2 + (cy - py)^2)`, subject to a minimum centre
//' separation of `(h_i + h_j)/2 + gap` between neighbours. Solved exactly for the given order by
//' a grid dynamic program (`O(n * gridSize)`): sweeping labels from the bottom up, `g[t]` is the
//' least total length with the current label at some slot at or below `t` and every lower label
//' packed beneath it, via a two-option transition per slot -- either place the label in slot `t`
//' (on top of the best lower label position) or skip the slot and place it lower down. Paired
//' with the Hungarian order (see the seed reorder) this
//' reproduces the assignment optimum -- exactly in the equal-height, single-line case -- because
//' it minimises the true Euclidean leader length rather than a squared-vertical proxy.
//'
//' The grid spans at least the whole viewport `[ylo, yhi]` (and the pole range), so labels may
//' use the full vertical space; when the stacked column is taller than that space the grid is
//' extended symmetrically beyond it so every box still fits.
//'
//' @param dx Numeric horizontal pole-to-column distances, in the fixed top-to-bottom order.
//' @param py Numeric pole y-coordinates, in the same order.
//' @param h Numeric full box heights along the stacking axis, in the same order.
//' @param gap Numeric extra separation added between neighbouring boxes.
//' @param slot Numeric grid resolution (slot height) for the candidate centre positions.
//' @param ylo,yhi Numeric viewport y-bounds the grid must cover.
//' @return Numeric vector of placed centre y-positions, aligned with the inputs.
//' @keywords internal
// [[Rcpp::export]]
NumericVector packLen(NumericVector dx, NumericVector py, NumericVector h,
                      double gap, double slot, double ylo, double yhi) {
  int n = py.size();
  NumericVector out(n);
  if (n == 1) {
    out[0] = py[0];
    return out;
  }
  // grid covers the viewport and the pole range; extend beyond it only if the stack cannot fit
  double lo = ylo;
  double hi = yhi;
  double stack = 0.0;                                   // total column height: sum(h) + gaps
  for (int i = 0; i < n; ++i) {
    stack += h[i];
    if (py[i] < lo) {
      lo = py[i];
    }
    if (py[i] > hi) {
      hi = py[i];
    }
  }
  stack += (n - 1) * gap;
  double extend = stack > (hi - lo) ? stack - (hi - lo) : 0.0;
  double gridLo = lo - extend / 2.0;
  double gridHi = hi + extend / 2.0;
  int slots = (int) std::ceil((gridHi - gridLo) / slot) + 1;
  std::vector<double> gridY(slots);
  for (int t = 0; t < slots; ++t) {
    gridY[t] = gridLo + t * slot;
  }
  // minimum centre separation between consecutive boxes, in whole grid steps
  std::vector<int> lag(n - 1);
  for (int k = 0; k < n - 1; ++k) {
    double sep = (h[k] + h[k + 1]) / 2.0 + gap;
    lag[k] = (int) std::ceil(sep / slot);
  }
  const double INF = std::numeric_limits<double>::infinity();
  // Bottom-up DP. g[k][t] = least total length with label k placed at some slot <= t and every
  // lower label packed below it. One sweep per label, with a two-option transition per slot:
  //   g[k][t] = min( leaderLen(k, t) + g[k+1][t - lag[k]],   // place label k in slot t, or
  //                  g[k][t - 1] )                            // skip slot t (label k sits lower)
  // Only two rows are kept (prev = g[k+1], cur = g[k]); placeHere[k][t] records which option won.
  std::vector<double> prev(slots), cur(slots);
  std::vector<std::vector<char> > placeHere(n, std::vector<char>(slots, 0));
  for (int k = n - 1; k >= 0; --k) {
    for (int t = 0; t < slots; ++t) {
      double d = gridY[t] - py[k];
      double place = std::sqrt(dx[k] * dx[k] + d * d);    // leader length for label k at slot t
      if (k < n - 1) {
        int src = t - lag[k];                             // label k+1 must sit at slot <= src
        place = (src >= 0) ? place + prev[src] : INF;
      }
      double skip = (t > 0) ? cur[t - 1] : INF;           // label k placed lower than slot t
      if (place <= skip) {
        cur[t] = place;
        placeHere[k][t] = (place < INF) ? 1 : 0;
      } else {
        cur[t] = skip;
      }
    }
    prev.swap(cur);                                       // prev now holds g[k]
  }
  // reconstruct: from the top label down, walk down past skipped slots to where each label sits
  int t = slots - 1;
  for (int k = 0; k < n; ++k) {
    while (!placeHere[k][t]) {
      --t;
    }
    out[k] = gridY[t];
    if (k < n - 1) {
      t -= lag[k];                                        // label k+1 sits at least lag[k] below
    }
  }
  return out;
}
