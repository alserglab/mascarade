#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

//' Hungarian (Jonker-Volgenant) assignment
//'
//' O(n^3) minimum-cost assignment on a square cost matrix. Used by the boundary seed (one
//' solve per column) to match labels to stacked slots.
//'
//' @param cost Square numeric cost matrix (`cost[i, j]` = cost of assigning row i to column j).
//' @return Integer vector `res` where `res[i]` is the 0-indexed column assigned to row i,
//'   minimising the total cost.
//' @keywords internal
// [[Rcpp::export]]
IntegerVector hungarian(NumericMatrix cost) {
  // Standard Jonker-Volgenant shortest-augmenting-path solver, kept in its canonical form.
  // Indices are 1-based (column 0 is a sentinel). u/v: row/column potentials; p[j]: the row
  // matched to column j; way[j]: the augmenting-path predecessor column of column j; minv[j]:
  // the minimum reduced cost to reach column j; used[j]: is column j already in the tree.
  int n = cost.nrow();
  const double INF = 1e18;
  std::vector<double> u(n + 1, 0), v(n + 1, 0);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(n + 1, INF);
    std::vector<char> used(n + 1, 0);
    // grow the alternating tree from row i until it reaches an unmatched column
    do {
      used[j0] = 1;
      int i0 = p[j0];
      int j1 = -1;
      double delta = INF;
      for (int j = 1; j <= n; ++j) {
        if (used[j]) {
          continue;
        }
        double cur = cost(i0 - 1, j - 1) - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] != 0);
    // walk the augmenting path back, flipping the matching along it
    do {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }
  IntegerVector res(n);
  for (int j = 1; j <= n; ++j) {
    res[p[j] - 1] = j - 1;
  }
  return res;
}
