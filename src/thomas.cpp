#include <Rcpp.h>
using namespace Rcpp;

// Fast Thomas process generator in a rectangle (domain bbox).
// Parents ~ Poisson(kappa * area), offspring per parent ~ Poisson(mu),
// offspring displacement ~ N(0, sigma^2) in x and y.
// We optionally early-stop once we have generated >= max_points (if max_points > 0).
//
// Returns an n x 2 matrix [x, y] inside the bounding box. Polygon filtering is done in R.
//
// [[Rcpp::export]]
NumericMatrix rthomas_bbox_cpp(double kappa,
                               double mu,
                               double sigma,
                               double xmin,
                               double ymin,
                               double xmax,
                               double ymax,
                               int    max_points = -1) {
  if (xmax <= xmin || ymax <= ymin) {
    stop("Invalid bbox in rthomas_bbox_cpp.");
  }
  if (kappa <= 0.0 || mu <= 0.0 || sigma <= 0.0) {
    // Return empty matrix rather than error; caller can handle
    return NumericMatrix(0, 2);
  }

  const double width  = xmax - xmin;
  const double height = ymax - ymin;
  const double area   = width * height;

  // Number of parents ~ Poisson(kappa * area)
  const double mean_parents = kappa * area;
  int n_par = R::rpois(mean_parents);
  if (n_par <= 0) {
    return NumericMatrix(0, 2);
  }

  // Pre-size a vector of points (we can over-allocate a bit if max_points < 0).
  // If max_points > 0, we will stop once we hit it.
  std::vector<double> xs;
  std::vector<double> ys;
  xs.reserve(max_points > 0 ? max_points : (int)std::ceil(mean_parents * mu * 1.2));
  ys.reserve(xs.capacity());

  for (int i = 0; i < n_par; ++i) {
    // Parent uniform in bbox
    const double px = xmin + unif_rand() * width;
    const double py = ymin + unif_rand() * height;

    // Offspring count ~ Poisson(mu)
    int n_off = R::rpois(mu);
    if (n_off <= 0) continue;

    for (int j = 0; j < n_off; ++j) {
      // Gaussian displacement
      const double dx = R::rnorm(0.0, sigma);
      const double dy = R::rnorm(0.0, sigma);
      const double x  = px + dx;
      const double y  = py + dy;

      // Keep only if in bbox (fast trivial reject)
      if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
        xs.push_back(x);
        ys.push_back(y);
        if (max_points > 0 && (int)xs.size() >= max_points) {
          // Early stop
          NumericMatrix out(xs.size(), 2);
          for (int k = 0; k < (int)xs.size(); ++k) {
            out(k, 0) = xs[k];
            out(k, 1) = ys[k];
          }
          return out;
        }
      }
    }
  }

  // Pack result
  NumericMatrix out(xs.size(), 2);
  for (int k = 0; k < (int)xs.size(); ++k) {
    out(k, 0) = xs[k];
    out(k, 1) = ys[k];
  }
  return out;
}
