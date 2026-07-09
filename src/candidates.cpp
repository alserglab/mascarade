#include "geometry.h"
#include <vector>
#include <set>
#include <utility>
using namespace Rcpp;

//' Radial free-space label candidates
//'
//' From each cluster pole, marches `ndir` rays outward and emits a candidate box centre
//' wherever a box of the label's size (plus `pad`) is cluster-free (a BoxFit R-tree query): at
//' the near edge of every free interval and every `intfill` along wide ones. Candidates may
//' fall partly outside the viewport -- the effective-length overflow term ranks them, so a
//' crowded label can take a minimally-clipped edge slot instead of the far seed.
//'
//' @param boxfit External pointer from `buildBoxFit()`.
//' @param poi K x 2 matrix of cluster poles (the ray origins).
//' @param hw,hh Numeric per-label box half-sizes.
//' @param pad Numeric hard box clearance added around each box.
//' @param xlo,xhi,ylo,yhi Numeric viewport bounds (kept for context; candidates are no longer
//'   clipped to them).
//' @param ndir Integer number of rays per pole.
//' @param step Numeric radial step along each ray.
//' @param rstart,rmax Numeric first and last radius searched.
//' @param intfill Numeric spacing of extra candidates along a wide free interval.
//' @param dedup Numeric grid size for de-duplicating nearby candidates.
//' @return A data.frame with integer `label` (1-indexed) and numeric `cx`, `cy`.
//' @keywords internal
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
        bool free = !bf->hit(x - hwi, x + hwi, y - hhi, y + hhi);   // cluster-free; viewport
                                                                    // overflow is priced later

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
