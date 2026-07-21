#include "geometry.h"
#include <utility>
using namespace Rcpp;

//' Build the box-fit R-tree for the cluster polygons
//'
//' Constructs a Boost.Geometry R-tree over the cluster polygon envelopes (and keeps the
//' polygons) so the candidate and polish kernels can answer "does this label box overlap any
//' cluster?" quickly. Rebuilt once per placement; the mask itself is not.
//'
//' @param polysx,polysy Lists of parallel numeric x/y vectors, one ring per cluster.
//' @return An external-pointer (`XPtr<BoxFit>`) handle to the box-fit structure.
//' @keywords internal
// [[Rcpp::export]]
SEXP buildBoxFit(List polysx, List polysy) {
  int nClusters = polysx.size();
  BoxFit* boxFit = new BoxFit();
  boxFit->polys.resize(nClusters);
  std::vector<mval> entries;                        // (envelope, cluster index) for the R-tree
  entries.reserve(nClusters);
  for (int k = 0; k < nClusters; ++k) {
    NumericVector vx = polysx[k];
    NumericVector vy = polysy[k];
    int m = vx.size();
    mpoly poly;
    for (int i = 0; i < m; ++i) {
      bg::append(poly.outer(), mpt(vx[i], vy[i]));
    }
    bg::correct(poly);                              // close the ring / fix orientation
    boxFit->polys[k] = poly;
    entries.push_back(std::make_pair(bg::return_envelope<mbox>(poly), (unsigned) k));
  }
  boxFit->tree = bgi::rtree<mval, bgi::quadratic<16> >(entries.begin(), entries.end());
  return XPtr<BoxFit>(boxFit, true);
}
