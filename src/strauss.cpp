// src/strauss.cpp
#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <cmath>
using namespace Rcpp;

// ---- Spatial hash for neighbor lookup ---------------------------------

struct HashGrid {
  double xmin, ymin, cell;
  std::unordered_map<long long, std::vector<int> > buckets;

  HashGrid(double xmin_, double ymin_, double cell_) : xmin(xmin_), ymin(ymin_), cell(cell_) {}

  inline long long key_from_xy(double x, double y) const {
    long long ix = (long long) std::floor((x - xmin) / cell);
    long long iy = (long long) std::floor((y - ymin) / cell);
    return (ix << 32) ^ (iy & 0xffffffffLL);
  }

  inline long long key_from_ij(long long ix, long long iy) const {
    return (ix << 32) ^ (iy & 0xffffffffLL);
  }

  inline void ij_from_xy(double x, double y, long long &ix, long long &iy) const {
    ix = (long long) std::floor((x - xmin) / cell);
    iy = (long long) std::floor((y - ymin) / cell);
  }

  void insert(int id, double x, double y) {
    long long k = key_from_xy(x, y);
    buckets[k].push_back(id);
  }

  void erase(int id, double x, double y) {
    long long k = key_from_xy(x, y);
    auto it = buckets.find(k);
    if (it == buckets.end()) return;
    auto &vec = it->second;
    for (size_t i = 0; i < vec.size(); ++i) {
      if (vec[i] == id) {
        vec[i] = vec.back();
        vec.pop_back();
        break;
      }
    }
    if (vec.empty()) buckets.erase(it);
  }

  void move(int id, double oldx, double oldy, double newx, double newy) {
    long long kold = key_from_xy(oldx, oldy);
    long long knew = key_from_xy(newx, newy);
    if (kold == knew) return;
    erase(id, oldx, oldy);
    insert(id, newx, newy);
  }
};

// count neighbors of point id at (x[id], y[id]) within radius r,
// using hash grid; optionally treat (qx,qy) as candidate instead.
// If candidate=true, use (qx,qy) instead of current (x[id],y[id]) and
// exclude id from neighbor lists (we pass its old cell anyway).
static int count_neighbors_for_point(
    int id,
    const std::vector<double> &x,
    const std::vector<double> &y,
    double qx, double qy, bool candidate,
    const HashGrid &grid,
    double r)
{
  const double r2 = r*r;
  long long ix, iy;
  grid.ij_from_xy(candidate ? qx : x[id], candidate ? qy : y[id], ix, iy);

  int cnt = 0;
  for (long long dx = -1; dx <= 1; ++dx) {
    for (long long dy = -1; dy <= 1; ++dy) {
      long long key = grid.key_from_ij(ix + dx, iy + dy);
      auto it = grid.buckets.find(key);
      if (it == grid.buckets.end()) continue;
      const auto &vec = it->second;
      for (int j : vec) {
        if (j == id) continue;
        double xx = candidate ? qx : x[id];
        double yy = candidate ? qy : y[id];
        double dx2 = x[j] - xx;
        double dy2 = y[j] - yy;
        if (dx2*dx2 + dy2*dy2 <= r2) ++cnt;
      }
    }
  }
  return cnt;
}

//' Fast Strauss simulation on a bounding box (fixed-n Metropolis-Hastings)
 //'
 //' Simulate n points in [xmin,xmax] × [ymin,ymax] from a Strauss model
 //' *conditioned on n*, i.e. density ∝ gamma^{s(X)} where s(X) is the
 //' number of unordered point pairs within distance r. This avoids the
 //' need for the beta (intensity) parameter and mixes quickly when n is
 //' moderate. We use single-point relocation proposals with a spatial
 //' hash for O(1) expected neighbor lookups.
 //'
 //' @param n Number of points (fixed).
 //' @param xmin,xmax,ymin,ymax Bounding box (numeric).
 //' @param r Interaction radius (>0).
 //' @param gamma Inhibition strength in (0,1]; smaller = stronger inhibition.
 //' @param sweeps Total MH sweeps (each sweep proposes ~n moves).
 //' @param burnin Number of initial sweeps to discard (still updates state).
 //' @param thin Keep every \code{thin}-th sweep after burnin (for diagnostics; final state is returned).
 //'
 //' @return A numeric matrix with two columns (x,y). The final state after the
 //'         last sweep is returned; thinning is currently not used to average.
 //'
 //' @note This is a *fixed-n* Strauss. If you need variable-n (with \code{beta}),
 //'       use a birth-death sampler; but for your package flow (target counts),
 //'       fixed-n is typically what you want.
 //'
 //' @export
 // [[Rcpp::export]]
 NumericMatrix rstrauss_bbox_cpp(
     int n,
     double xmin, double xmax,
     double ymin, double ymax,
     double r,
     double gamma,
     int sweeps = 2000,
     int burnin = 200,
     int thin   = 1
 ) {
   if (n < 0) stop("n must be >= 0");
   if (n == 0) {
     NumericMatrix out(0, 2);
     colnames(out) = CharacterVector::create("x","y");
     return out;
   }
   if (!(gamma > 0.0 && gamma <= 1.0)) stop("gamma must be in (0, 1]");
   if (!(r > 0.0)) stop("r must be > 0");
   if (xmax <= xmin || ymax <= ymin) stop("invalid bbox");

   // init uniform
   std::vector<double> x(n), y(n);
   for (int i = 0; i < n; ++i) {
     x[i] = R::runif(xmin, xmax);
     y[i] = R::runif(ymin, ymax);
   }

   // hash grid with cell size = r
   HashGrid grid(xmin, ymin, r);
   for (int i = 0; i < n; ++i) grid.insert(i, x[i], y[i]);

   const double loggamma = std::log(gamma);

   auto attempt_move = [&](int i) {
     // count current neighbors of i
     int cur = count_neighbors_for_point(i, x, y, 0.0, 0.0, false, grid, r);

     // propose new uniform position
     double nx = R::runif(xmin, xmax);
     double ny = R::runif(ymin, ymax);
     int nxt = count_neighbors_for_point(i, x, y, nx, ny, true, grid, r);

     int delta = nxt - cur;
     double logacc = delta * loggamma; // acceptance on Strauss energy
     bool accept = false;
     if (logacc >= 0) {
       accept = true;
     } else {
       double u = R::runif(0.0, 1.0);
       accept = (std::log(u) < logacc);
     }
     if (accept) {
       // Update hash and position
       grid.move(i, x[i], y[i], nx, ny);
       x[i] = nx; y[i] = ny;
     }
   };

   // MH sweeps
   for (int s = 0; s < sweeps; ++s) {
     for (int k = 0; k < n; ++k) {
       int i = (int) std::floor(R::runif(0.0, (double) n));
       if (i >= n) i = n - 1;
       attempt_move(i);
     }
     // (Optional) we could collect states every 'thin', but we return final state.
     (void)burnin; (void)thin;
     Rcpp::checkUserInterrupt();
   }

   NumericMatrix out(n, 2);
   for (int i = 0; i < n; ++i) {
     out(i,0) = x[i];
     out(i,1) = y[i];
   }
   colnames(out) = CharacterVector::create("x","y");
   return out;
 }
