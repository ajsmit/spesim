// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
using namespace Rcpp;

// --- Helpers ---------------------------------------------------------------

struct CellGrid {
  double xmin, ymin, cell;
  int nx, ny;
  std::vector< std::vector<int> > buckets; // indices per cell

  CellGrid(double xmin_, double xmax_, double ymin_, double ymax_, double cell_)
    : xmin(xmin_), ymin(ymin_), cell(cell_) {
    nx = std::max(1, (int)std::floor((xmax_ - xmin_) / cell_));
    ny = std::max(1, (int)std::floor((ymax_ - ymin_) / cell_));
    buckets.assign(nx * ny, std::vector<int>());
  }

  inline int idx(int ix, int iy) const {
    return iy * nx + ix;
  }
  inline void clamp_cell(int &ix, int &iy) const {
    if (ix < 0) ix = 0; else if (ix >= nx) ix = nx - 1;
    if (iy < 0) iy = 0; else if (iy >= ny) iy = ny - 1;
  }
  inline void locate(double x, double y, int &ix, int &iy) const {
    ix = (int) std::floor((x - xmin) / cell);
    iy = (int) std::floor((y - ymin) / cell);
    clamp_cell(ix, iy);
  }
  inline void insert(int p, double x, double y) {
    int ix, iy; locate(x,y,ix,iy);
    buckets[idx(ix,iy)].push_back(p);
  }
  inline void remove_from_cell(int p, int ix, int iy) {
    std::vector<int> &v = buckets[idx(ix,iy)];
    for (size_t k=0; k<v.size(); ++k) if (v[k]==p) { v[k] = v.back(); v.pop_back(); break; }
  }
  inline void move(int p, double oldx,double oldy, double newx,double newy) {
    int ox, oy; locate(oldx,oldy,ox,oy);
    int nx_, ny_; locate(newx,newy,nx_,ny_);
    if (ox==nx_ && oy==ny_) return;
    remove_from_cell(p, ox, oy);
    buckets[idx(nx_,ny_)].push_back(p);
  }
};

// Count neighbors of (x0,y0) within r using cell grid; optionally skip index 'skip'
inline int count_neighbors_local(const std::vector<double>& X,
                                 const std::vector<double>& Y,
                                 const CellGrid& grid,
                                 double x0, double y0,
                                 double r2,
                                 int skip = -1) {
  int ix, iy; grid.locate(x0,y0,ix,iy);
  int nbh = 0;
  // loop 3x3 neighborhood of cells
  for (int dy=-1; dy<=1; ++dy) {
    for (int dx=-1; dx<=1; ++dx) {
      int cx = ix + dx, cy = iy + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx,cy) ];
      for (int j : bucket) {
        if (j == skip) continue;
        double dx_ = X[j] - x0;
        double dy_ = Y[j] - y0;
        double d2 = dx_*dx_ + dy_*dy_;
        if (d2 <= r2) ++nbh;
      }
    }
  }
  return nbh;
}

// Update per-point neighbor counts around a point moved from (x_old,y_old) to (x_new,y_new)
// Only points within r of either location are touched.
// Returns deltaPairs = (#pairs_new - #pairs_old) for the moved point.
inline int update_counts_move(std::vector<int>& k,
                              const std::vector<double>& X,
                              const std::vector<double>& Y,
                              const CellGrid& grid,
                              int i,
                              double x_old,double y_old,
                              double x_new,double y_new,
                              double r2) {
  int ix_old, iy_old; grid.locate(x_old,y_old,ix_old,iy_old);
  int ix_new, iy_new; grid.locate(x_new,y_new,ix_new,iy_new);

  // Mark neighbors in union of neighborhoods of old/new cells
  // We will adjust the neighbor counts of affected points only once.
  // Use a small local set (vector<char>) sized by total points? n may be large; instead handle both neighborhoods separately.

  // First: handle neighbors near old position (pairs being removed)
  for (int dy=-1; dy<=1; ++dy) {
    for (int dx=-1; dx<=1; ++dx) {
      int cx = ix_old + dx, cy = iy_old + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx,cy) ];
      for (int j : bucket) {
        if (j == i) continue;
        double dx_ = X[j] - x_old;
        double dy_ = Y[j] - y_old;
        double d2  = dx_*dx_ + dy_*dy_;
        if (d2 <= r2) {
          // was a neighbor at old location; will be removed unless still neighbor at new
          double dxn = X[j] - x_new;
          double dyn = Y[j] - y_new;
          double d2n = dxn*dxn + dyn*dyn;
          if (d2n > r2) {
            // break the pair i-j
            --k[j];
          }
        }
      }
    }
  }

  // Second: handle neighbors near new position (pairs being added)
  for (int dy=-1; dy<=1; ++dy) {
    for (int dx=-1; dx<=1; ++dx) {
      int cx = ix_new + dx, cy = iy_new + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy < 0 || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx,cy) ];
      for (int j : bucket) {
        if (j == i) continue;
        double dxn = X[j] - x_new;
        double dyn = Y[j] - y_new;
        double d2n = dxn*dxn + dyn*dyn;
        if (d2n <= r2) {
          // create the pair i-j unless it already existed at old location
          double dx_ = X[j] - x_old;
          double dy_ = Y[j] - y_old;
          double d2  = dx_*dx_ + dy_*dy_;
          if (d2 > r2) {
            ++k[j];
          }
        }
      }
    }
  }

  // Compute delta for point i itself: new count minus old count
  int old_cnt = count_neighbors_local(X, Y, grid, x_old, y_old, r2, i);
  int new_cnt = count_neighbors_local(X, Y, grid, x_new, y_new, r2, i);
  int delta_i = new_cnt - old_cnt;
  return delta_i;
}

// Reflect to bbox
inline double clamp(double v, double lo, double hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

// Sample a displacement uniformly inside a disk of radius 'step'
inline void draw_uniform_disk(double step, double &dx, double &dy) {
  double u = R::unif_rand();
  double r = step * std::sqrt(u);
  double th = 2.0 * M_PI * R::unif_rand();
  dx = r * std::cos(th);
  dy = r * std::sin(th);
}

// [[Rcpp::export]]
NumericMatrix rstrauss_bbox_cpp(int n,
                                double xmin, double xmax,
                                double ymin, double ymax,
                                double r, double gamma,
                                int sweeps = 2000,
                                int burnin = 200,
                                int thin = 1,
                                double step0 = NA_REAL,
                                double target_acc = 0.35,
                                int tune_every = 200) {
  if (n <= 0) return NumericMatrix(0,2);
  if (!(gamma > 0.0 && gamma <= 1.0))
    stop("gamma must be in (0,1].");
  if (!(xmax > xmin && ymax > ymin))
    stop("Invalid bbox.");
  if (!(r > 0.0))
    stop("r must be > 0.");

  RNGScope scope;

  const double r2 = r * r;
  const double cell = r; // cell width >= r
  CellGrid grid(xmin, xmax, ymin, ymax, cell);

  // --- Initialise points uniformly -----------------------------------
  std::vector<double> X(n), Y(n);
  for (int i=0; i<n; ++i) {
    X[i] = xmin + (xmax - xmin) * R::unif_rand();
    Y[i] = ymin + (ymax - ymin) * R::unif_rand();
    grid.insert(i, X[i], Y[i]);
  }

  // --- Per-point neighbor counts k[i] --------------------------------
  std::vector<int> k(n, 0);
  for (int i=0; i<n; ++i) {
    // count neighbors (excluding i); note each pair will be counted twice in sum(k)
    k[i] = count_neighbors_local(X, Y, grid, X[i], Y[i], r2, i);
  }

  // --- MH parameters --------------------------------------------------
  double step = R_finite(step0) ? step0 : (0.5 * r);
  const double log_gamma = std::log(gamma);

  // number of proposals per sweep: n (one per point on average)
  const long long props_per_sweep = (long long) n;
  long long props_since_tune = 0LL;
  int accepted_since_tune = 0;

  const int total_sweeps = std::max(sweeps, burnin);

  for (int s=0; s<total_sweeps; ++s) {
    for (int p=0; p<n; ++p) {
      // pick random point
      int i = (int) std::floor(n * R::unif_rand());
      if (i == n) i = n - 1;

      // propose local move
      double dx, dy; draw_uniform_disk(step, dx, dy);
      double xn = clamp(X[i] + dx, xmin, xmax);
      double yn = clamp(Y[i] + dy, ymin, ymax);
      if (xn == X[i] && yn == Y[i]) {
        // degenerate (hit boundary exactly); skip
        continue;
      }

      // compute delta for counts (affects k[i] and neighbors)
      int delta_i = update_counts_move(k, X, Y, grid, i, X[i], Y[i], xn, yn, r2);
      // acceptance in log-space: Δlogπ = Δpairs * log(gamma),
      // where Δpairs = new_pairs - old_pairs for point i; all neighbor deltas
      // already applied in update_counts_move
      double dlog = ((double)delta_i) * log_gamma;

      bool accept = false;
      if (dlog >= 0.0) {
        accept = true;
      } else {
        double u = R::unif_rand();
        accept = (std::log(u) < dlog);
      }

      if (accept) {
        // Move accepted: update point and grid, and k[i]
        grid.move(i, X[i], Y[i], xn, yn);
        X[i] = xn; Y[i] = yn;
        k[i] += delta_i;

        ++accepted_since_tune;
      }
      ++props_since_tune;

      // Tuning (only during burn-in sweeps)
      if ((s < burnin) && (props_since_tune >= props_per_sweep * (long long) tune_every)) {
        double acc = (double)accepted_since_tune / (double)props_since_tune;
        if (acc > target_acc + 0.05)      step *= 1.2;
        else if (acc < target_acc - 0.05) step *= 0.8;
        // keep step sensible
        if (step < 1e-6) step = 1e-6;
        if (step > (xmax - xmin)) step = (xmax - xmin);
        accepted_since_tune = 0;
        props_since_tune = 0LL;
      }
    }
  }

  // Return final configuration
  NumericMatrix out(n, 2);
  for (int i=0; i<n; ++i) {
    out(i,0) = X[i];
    out(i,1) = Y[i];
  }
  colnames(out) = CharacterVector::create("x","y");
  return out;
}
