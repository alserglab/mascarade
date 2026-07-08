#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// O(n^3) Hungarian / Jonker-Volgenant assignment on a square cost matrix.
// Returns res[i] = 0-indexed column assigned to row i, minimizing total cost.
// Used by the min-cost boundary seed (one Hungarian solve per column).
// [[Rcpp::export]]
IntegerVector hungarian(NumericMatrix cost) {
  int n = cost.nrow();
  const double INF = 1e18;
  std::vector<double> u(n + 1, 0), v(n + 1, 0);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i; int j0 = 0;
    std::vector<double> minv(n + 1, INF);
    std::vector<char> used(n + 1, 0);
    do {
      used[j0] = 1; int i0 = p[j0], j1 = -1; double delta = INF;
      for (int j = 1; j <= n; ++j) if (!used[j]) {
        double cur = cost(i0 - 1, j - 1) - u[i0] - v[j];
        if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
        if (minv[j] < delta) { delta = minv[j]; j1 = j; }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
        else minv[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] != 0);
    do { int j1 = way[j0]; p[j0] = p[j1]; j0 = j1; } while (j0);
  }
  IntegerVector res(n);
  for (int j = 1; j <= n; ++j) res[p[j] - 1] = j - 1;
  return res;
}
