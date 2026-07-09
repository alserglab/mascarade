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
  int np = polysx.size();
  BoxFit* bf = new BoxFit();
  bf->polys.resize(np);
  std::vector<mval> vals;
  vals.reserve(np);
  for (int k = 0; k < np; ++k) {
    NumericVector vx = polysx[k], vy = polysy[k];
    int m = vx.size();
    mpoly p;
    for (int i = 0; i < m; ++i) bg::append(p.outer(), mpt(vx[i], vy[i]));
    bg::correct(p);
    bf->polys[k] = p;
    vals.push_back(std::make_pair(bg::return_envelope<mbox>(p), (unsigned) k));
  }
  bf->tree = bgi::rtree<mval, bgi::quadratic<16> >(vals.begin(), vals.end());
  XPtr<BoxFit> ptr(bf, true);
  return ptr;
}
