// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;

//' Fast vectorised point-in-polygon test (ray casting)
//'
//' @param px Numeric vector of query x-coordinates.
//' @param py Numeric vector of query y-coordinates.
//' @param rx Numeric vector of polygon ring x-coordinates (closed ring: first equals last).
//' @param ry Numeric vector of polygon ring y-coordinates (closed ring: first equals last).
//' @return LogicalVector; TRUE if the corresponding point is inside the ring.
//' @keywords internal
// [[Rcpp::export]]
LogicalVector pip_cpp(NumericVector px, NumericVector py,
                      NumericVector rx, NumericVector ry) {
  int n = px.size();
  int m = rx.size() - 1;   // number of edges (ring closed: first == last)
  LogicalVector inside(n, false);

  for (int i = 0; i < n; i++) {
    double x = px[i], y = py[i];
    bool in_poly = false;
    for (int j = 0, k = m - 1; j < m; k = j++) {
      double xj = rx[j], yj = ry[j];
      double xk = rx[k], yk = ry[k];
      if (((yj > y) != (yk > y)) &&
          (x < (xk - xj) * (y - yj) / (yk - yj) + xj)) {
        in_poly = !in_poly;
      }
    }
    inside[i] = in_poly;
  }
  return inside;
}
