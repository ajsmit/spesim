// src/geyer.cpp
#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Fixed-n Metropolis–Hastings sampler for a (heuristic but internally consistent)
// Geyer saturation process in a rectangular window.
//
// Target (up to constant):  pi(X) ∝ gamma^{ SUM_i min(sat, n_i) }
// where n_i is the number of neighbours of point i within radius r.
//
// IMPORTANT: acceptance ratio must account for changes in neighbours' n_j as
// point i moves. We therefore maintain cached neighbour counts and apply deltas
// only on ACCEPT.
// -----------------------------------------------------------------------------

struct CellGrid {
  double xmin, ymin, cell;
  int nx, ny;
  std::vector< std::vector<int> > buckets;

  CellGrid(double xmin_, double xmax_, double ymin_, double ymax_, double cell_)
    : xmin(xmin_), ymin(ymin_), cell(cell_) {
    nx = std::max(1, (int)std::floor((xmax_ - xmin_) / cell_));
    ny = std::max(1, (int)std::floor((ymax_ - ymin_) / cell_));
    buckets.assign(nx * ny, std::vector<int>());
  }

  inline int idx(int ix, int iy) const { return iy * nx + ix; }

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
    int ix, iy; locate(x, y, ix, iy);
    buckets[idx(ix, iy)].push_back(p);
  }

  inline void remove_from_cell(int p, int ix, int iy) {
    std::vector<int> &v = buckets[idx(ix, iy)];
    for (size_t k = 0; k < v.size(); ++k) {
      if (v[k] == p) {
        v[k] = v.back();
        v.pop_back();
        break;
      }
    }
  }

  inline void move(int p, double oldx, double oldy, double newx, double newy) {
    int ox, oy; locate(oldx, oldy, ox, oy);
    int nx_, ny_; locate(newx, newy, nx_, ny_);
    if (ox == nx_ && oy == ny_) return;
    remove_from_cell(p, ox, oy);
    buckets[idx(nx_, ny_)].push_back(p);
  }
};

inline double clamp(double v, double lo, double hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

// Count neighbours of (x0,y0) within r^2 using cell grid; optionally skip index.
inline int count_neighbors_local(const std::vector<double>& X,
                                 const std::vector<double>& Y,
                                 const CellGrid& grid,
                                 double x0, double y0,
                                 double r2,
                                 int skip = -1) {
  int ix, iy; grid.locate(x0, y0, ix, iy);
  int nbh = 0;
  for (int dy = -1; dy <= 1; ++dy) {
    for (int dx = -1; dx <= 1; ++dx) {
      int cx = ix + dx, cy = iy + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx, cy) ];
      for (int j : bucket) {
        if (j == skip) continue;
        double dx_ = X[j] - x0;
        double dy_ = Y[j] - y0;
        if (dx_*dx_ + dy_*dy_ <= r2) ++nbh;
      }
    }
  }
  return nbh;
}

struct MoveDelta {
  int old_cnt;
  int new_cnt;
  std::vector< std::pair<int,int> > neighbor_delta; // (index, +/-1)
};

inline int sat_contrib(int cnt, int sat) {
  if (sat <= 0) return 0;
  return (cnt >= sat) ? sat : cnt;
}

inline MoveDelta propose_move_deltas(const std::vector<double>& X,
                                    const std::vector<double>& Y,
                                    const CellGrid& grid,
                                    int i,
                                    double x_old, double y_old,
                                    double x_new, double y_new,
                                    double r2) {
  MoveDelta md;
  md.neighbor_delta.clear();

  md.old_cnt = count_neighbors_local(X, Y, grid, x_old, y_old, r2, i);
  md.new_cnt = count_neighbors_local(X, Y, grid, x_new, y_new, r2, i);

  std::unordered_map<int,int> dmap;

  int ix_old, iy_old; grid.locate(x_old, y_old, ix_old, iy_old);
  int ix_new, iy_new; grid.locate(x_new, y_new, ix_new, iy_new);

  // neighbours near old position: pair i-j may be removed
  for (int dy=-1; dy<=1; ++dy) {
    for (int dx=-1; dx<=1; ++dx) {
      int cx = ix_old + dx, cy = iy_old + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx, cy) ];
      for (int j : bucket) {
        if (j == i) continue;
        double dx_ = X[j] - x_old;
        double dy_ = Y[j] - y_old;
        if (dx_*dx_ + dy_*dy_ <= r2) {
          // was neighbour at old; still neighbour at new?
          double dxn = X[j] - x_new;
          double dyn = Y[j] - y_new;
          if (dxn*dxn + dyn*dyn > r2) {
            dmap[j] -= 1;
          }
        }
      }
    }
  }

  // neighbours near new position: pair i-j may be added
  for (int dy=-1; dy<=1; ++dy) {
    for (int dx=-1; dx<=1; ++dx) {
      int cx = ix_new + dx, cy = iy_new + dy;
      if (cx < 0 || cy < 0 || cx >= grid.nx || cy >= grid.ny) continue;
      const std::vector<int> &bucket = grid.buckets[ grid.idx(cx, cy) ];
      for (int j : bucket) {
        if (j == i) continue;
        double dxn = X[j] - x_new;
        double dyn = Y[j] - y_new;
        if (dxn*dxn + dyn*dyn <= r2) {
          // is neighbour at new; was it neighbour at old?
          double dx_ = X[j] - x_old;
          double dy_ = Y[j] - y_old;
          if (dx_*dx_ + dy_*dy_ > r2) {
            dmap[j] += 1;
          }
        }
      }
    }
  }

  md.neighbor_delta.reserve(dmap.size());
  for (auto &kv : dmap) {
    if (kv.second != 0) {
      md.neighbor_delta.emplace_back(kv.first, kv.second);
    }
  }
  return md;
}

// [[Rcpp::export]]
List rgeyer_bbox_cpp(int n_target,
                     double xmin, double xmax,
                     double ymin, double ymax,
                     double r, double gamma, int sat,
                     int sweeps = 2000, int burnin = 200, int thin = 1) {
  if (n_target <= 0) return List::create(_["x"] = NumericVector(0), _["y"] = NumericVector(0));
  if (!(xmax > xmin && ymax > ymin)) stop("Invalid bbox.");
  if (!(r > 0.0)) stop("r must be > 0.");
  if (!(gamma > 0.0)) stop("gamma must be > 0.");
  sat = std::max(0, sat);

  RNGScope scope;

  const double r2 = r * r;
  const double cell = r; // >= r
  CellGrid grid(xmin, xmax, ymin, ymax, cell);

  std::vector<double> X(n_target), Y(n_target);
  for (int i = 0; i < n_target; ++i) {
    X[i] = xmin + (xmax - xmin) * R::unif_rand();
    Y[i] = ymin + (ymax - ymin) * R::unif_rand();
    grid.insert(i, X[i], Y[i]);
  }

  // Cached neighbour counts
  std::vector<int> nbh(n_target, 0);
  for (int i = 0; i < n_target; ++i) {
    nbh[i] = count_neighbors_local(X, Y, grid, X[i], Y[i], r2, i);
  }

  const double log_gamma = std::log(std::max(gamma, 1e-12));

  int total_sweeps = sweeps + burnin;
  double prop_sd = 0.25 * r;

  for (int it = 0; it < total_sweeps; ++it) {
    for (int step = 0; step < n_target; ++step) {
      int i = (int) std::floor(n_target * R::unif_rand());
      if (i == n_target) i = n_target - 1;

      // propose local move (reflect at edges)
      double xn = R::rnorm(X[i], prop_sd);
      double yn = R::rnorm(Y[i], prop_sd);
      if (xn < xmin) xn = xmin + (xmin - xn);
      if (xn > xmax) xn = xmax - (xn - xmax);
      if (yn < ymin) yn = ymin + (ymin - yn);
      if (yn > ymax) yn = ymax - (yn - ymax);

      if (xn == X[i] && yn == Y[i]) continue;

      MoveDelta md = propose_move_deltas(X, Y, grid, i, X[i], Y[i], xn, yn, r2);

      // Compute delta in the saturation potential SUM min(sat, n_j)
      int deltaU = 0;
      deltaU += sat_contrib(md.new_cnt, sat) - sat_contrib(md.old_cnt, sat);
      for (const auto &pr : md.neighbor_delta) {
        int j = pr.first;
        int dj = pr.second;
        int oldj = nbh[j];
        int newj = oldj + dj;
        deltaU += sat_contrib(newj, sat) - sat_contrib(oldj, sat);
      }

      double dlog = ((double)deltaU) * log_gamma;
      bool accept = false;
      if (dlog >= 0.0) {
        accept = true;
      } else {
        accept = (std::log(R::unif_rand()) < dlog);
      }

      if (accept) {
        // apply neighbour deltas
        for (const auto &pr : md.neighbor_delta) {
          nbh[pr.first] += pr.second;
        }
        // update moved point and its neighbour count
        grid.move(i, X[i], Y[i], xn, yn);
        X[i] = xn; Y[i] = yn;
        nbh[i] = md.new_cnt;
      }
    }
    (void)thin;
  }

  NumericVector Xr(n_target), Yr(n_target);
  for (int i = 0; i < n_target; ++i) {
    Xr[i] = X[i];
    Yr[i] = Y[i];
  }
  return List::create(_["x"] = Xr, _["y"] = Yr);
}
