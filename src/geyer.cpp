// src/geyer.cpp
#include <Rcpp.h>
using namespace Rcpp;

// Count saturated neighbors within radius r (capped by 'sat')
inline int neigh_sat(double x, double y,
                     const NumericVector& X,
                     const NumericVector& Y,
                     double r, int sat) {
  const double r2 = r * r;
  int cnt = 0;
  for (int i = 0; i < X.size(); ++i) {
    double dx = X[i] - x;
    double dy = Y[i] - y;
    if (dx*dx + dy*dy <= r2) {
      ++cnt;
      if (cnt >= sat) return sat;
    }
  }
  return cnt;
}

// [[Rcpp::export]]
List rgeyer_bbox_cpp(int n_target,
                     double xmin, double xmax,
                     double ymin, double ymax,
                     double r, double gamma, int sat,
                     int sweeps = 2000, int burnin = 200, int thin = 1) {
  if (n_target <= 0) return List::create(_["x"]=NumericVector(0), _["y"]=NumericVector(0));

  NumericVector X(n_target), Y(n_target);

  // init uniform
  for (int i = 0; i < n_target; ++i) {
    X[i] = R::runif(xmin, xmax);
    Y[i] = R::runif(ymin, ymax);
  }

  // MH updates: random single-point move proposals
  const double area = (xmax - xmin) * (ymax - ymin);
  if (area <= 0) stop("Invalid bbox area.");

  int total = sweeps + burnin;
  double prop_sd = 0.25 * r; // local proposals

  for (int it = 0; it < total; ++it) {
    for (int k = 0; k < n_target; ++k) {
      // current saturated neighbor count around X[k],Y[k] excluding self
      NumericVector Xex = clone(X), Yex = clone(Y);
      // remove k by swapping with last and resizing
      int last = n_target - 1;
      if (k != last) {
        Xex[k] = X[last]; Yex[k] = Y[last];
      }
      Xex = Xex[Range(0, last-1)];
      Yex = Yex[Range(0, last-1)];

      int s_old = neigh_sat(X[k], Y[k], Xex, Yex, r, sat);

      // propose local move (reflect at edges)
      double nx = R::rnorm(X[k], prop_sd);
      double ny = R::rnorm(Y[k], prop_sd);
      if (nx < xmin) nx = xmin + (xmin - nx);
      if (nx > xmax) nx = xmax - (nx - xmax);
      if (ny < ymin) ny = ymin + (ymin - ny);
      if (ny > ymax) ny = ymax - (ny - ymax);

      int s_new = neigh_sat(nx, ny, Xex, Yex, r, sat);

      // Geyer density up to constant: gamma^{sum min(sat, n_i)}
      double logalpha = (double)(s_new - s_old) * std::log(std::max(gamma, 1e-12));
      if (log(R::runif(0.0, 1.0)) < logalpha) {
        X[k] = nx; Y[k] = ny;
      }
    }
    // thinning not strictly necessary here; we return final config
    (void)thin;
  }

  return List::create(_["x"]=X, _["y"]=Y);
}
