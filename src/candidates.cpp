#include "geometry.h"
#include <vector>
#include <set>
#include <utility>
using namespace Rcpp;

// Radial free-space candidate generation. From each label's pole, march ndir rays; a box
// of the label's size (+pad) fits at (x,y) iff it overlaps no cluster polygon (BoxFit
// R-tree query, replacing the prototype's occupancy integral image). Drop a candidate at
// the near edge of every free interval + every intfill along wide ones.
// [[Rcpp::export]]
DataFrame radialCandidates(SEXP boxfit, NumericMatrix poi, NumericVector hw, NumericVector hh, double pad,
                           double xlo, double xhi, double ylo, double yhi,
                           int ndir, double step, double rstart, double rmax, double intfill, double dedup) {
  XPtr<BoxFit> bf(boxfit);
  int K = poi.nrow();
  std::vector<int> Lo; std::vector<double> Xo, Yo;
  for (int i = 0; i < K; ++i) {
    double px = poi(i, 0), py = poi(i, 1), hwi = hw[i] + pad, hhi = hh[i] + pad;
    std::set<std::pair<long, long> > seen;
    for (int k = 0; k < ndir; ++k) {
      double th = 2 * M_PI * k / ndir, cx = std::cos(th), cy = std::sin(th);
      bool prevfree = false; double lastadd = -1e18, r = rstart;
      while (r <= rmax) {
        double x = px + r * cx, y = py + r * cy;
        bool inb = x - hw[i] >= xlo && x + hw[i] <= xhi && y - hh[i] >= ylo && y + hh[i] <= yhi;
        bool free = inb && !bf->hit(x - hwi, x + hwi, y - hhi, y + hhi);
        if (free && (!prevfree || r - lastadd >= intfill)) {
          long kx = (long) std::lround(x / dedup), ky = (long) std::lround(y / dedup);
          if (seen.insert(std::make_pair(kx, ky)).second) { Lo.push_back(i + 1); Xo.push_back(x); Yo.push_back(y); }
          lastadd = r;
        }
        prevfree = free; r += step;
      }
    }
  }
  return DataFrame::create(Named("label") = Lo, Named("cx") = Xo, Named("cy") = Yo);
}
