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
//' @param ndir Integer number of rays per pole.
//' @param step Numeric radial step along each ray.
//' @param rstart,rmax Numeric first and last radius searched.
//' @param intfill Numeric spacing of extra candidates along a wide free interval.
//' @param dedup Numeric grid size for de-duplicating nearby candidates.
//' @return A data.frame with integer `label` (1-indexed) and numeric `cx`, `cy`.
//' @keywords internal
// [[Rcpp::export]]
DataFrame radialCandidates(SEXP boxfit, NumericMatrix poi,
                           NumericVector hw, NumericVector hh, double pad,
                           int ndir, double step, double rstart, double rmax,
                           double intfill, double dedup) {
  XPtr<BoxFit> boxFit(boxfit);
  int nLabels = poi.nrow();
  std::vector<int> outLabel;
  std::vector<double> outX, outY;
  for (int i = 0; i < nLabels; ++i) {
    double poleX = poi(i, 0);
    double poleY = poi(i, 1);
    double halfW = hw[i] + pad;
    double halfH = hh[i] + pad;
    std::set<std::pair<long, long> > seen;                 // dedup grid keys already emitted
    for (int k = 0; k < ndir; ++k) {
      double theta = 2 * M_PI * k / ndir;
      double dirX = std::cos(theta);
      double dirY = std::sin(theta);
      bool prevFree = false;
      double lastAdd = -1e18;
      for (double r = rstart; r <= rmax; r += step) {
        double x = poleX + r * dirX;
        double y = poleY + r * dirY;
        bool free = !boxFit->hit(Rect{x - halfW, x + halfW, y - halfH, y + halfH});
        // emit at the near edge of a free interval, then every `intfill` along a wide one
        // (candidates may fall outside the viewport; the effective-length overflow ranks them)
        if (free && (!prevFree || r - lastAdd >= intfill)) {
          long keyX = (long) std::lround(x / dedup);
          long keyY = (long) std::lround(y / dedup);
          if (seen.insert(std::make_pair(keyX, keyY)).second) {
            outLabel.push_back(i + 1);
            outX.push_back(x);
            outY.push_back(y);
          }
          lastAdd = r;
        }
        prevFree = free;
      }
    }
  }
  return DataFrame::create(Named("label") = outLabel, Named("cx") = outX, Named("cy") = outY);
}
