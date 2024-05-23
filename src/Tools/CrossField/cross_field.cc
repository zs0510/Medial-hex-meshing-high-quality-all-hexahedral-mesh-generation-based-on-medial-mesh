//
// Created by su on 2022/3/13.
//

#include "cross_field.h"

#include <cfloat>

#include <numeric>
#include <map>
#include <queue>
#include <algorithm>
#include <unordered_map>

#include "basic_types.h"
#include "basic_utils.h"

#include "linear_solver.h"

#ifndef M_PI
#define M_PI 3.14159265358979
#endif // M_PI

using std::vector;
using std::map;
using std::unordered_map;
using std::array;
using std::pair;

namespace {
struct IJV {
  ID2 ij;
  double val;
  bool operator<(const IJV &b) const {
    return ij[0] < b.ij[0] || (ij[0] == b.ij[0] && ij[1] < b.ij[1]);
  }
};

struct IV {
  ID i;
  double val;
  bool operator<(const IV &b) const { return i < b.i; }
};

/* two unknowns (x_2i, x_2i+1) for each edge */
bool stiffness_coefficient(int N,
                           const std::vector<Vec3> &points,
                           const std::vector<ID3> &triangles,
                           ID e,
                           const std::vector<ID2> &uIEdges,
                           const std::vector<ID> &old2IEdge,
                           const std::vector<std::vector<ID>> &uIEdgeToOld,
                           const std::vector<Vec3> &triangle_normal,
                           std::vector<IV> &iv,
                           std::vector<IJV> &ijv) {
  if (uIEdgeToOld[e].size() != 2) {
    error(__FILE__,
          __LINE__,
          "assembly, edge {}: uIEdgeToOld[e].size() = {}",
          e,
          uIEdgeToOld[e].size());
    return false;
  }
  /* edge vertices */
  ID v1 = uIEdges[e][0];
  ID v2 = uIEdges[e][1];
  Vec3 e_x = points[v2] - points[v1];
  double len_r = Length(e_x);
  if (len_r < EPS) {
    error(__FILE__,
          __LINE__,
          "edge too small: v1={}, v2={}, Length = {}",
          v1,
          v2,
          len_r);
    return false;
  }
  e_x = (1. / len_r) * e_x;
  /* four neighbor edges */
  ID b_vars[4] = {NO_ID, NO_ID, NO_ID, NO_ID};
  double alpha[4] = {0., 0., 0., 0.};
  double CR_weight[4] = {-1. / 4, -1. / 4, -1. / 4, -1. / 4};
  Vec3 prevN;
  for (int s = 0; s < 2; ++s) {
    ID oe = uIEdgeToOld[e][s];
    ID t = oe / 3;
    ID le = oe % 3;
    Vec3 tn = triangle_normal[t];
    if (s == 1 && Dot(prevN, tn) < 0.) tn = -1. * tn;
    prevN = tn;
    // tn = {0,0,1};
    Vec3 e_y = Cross(tn, e_x);
    if (MaxAbs(e_y) == 0.) {
      error(__FILE__, __LINE__, "Length(e_y) = {}", Length(e_y));
      return false;
    }
    Normalize(e_y);
    for (int k = 0; k < 2; ++k) { /* two other edges in triangle t */
      ID aoe = 3 * t + (le + 1 + k) % 3;
      ID ae = old2IEdge[aoe];
      ID2 a_edge = uIEdges[ae];
      b_vars[2 * s + k] = ae;
      Vec3 edg = points[a_edge[1]] - points[a_edge[0]];
      double len = Length(edg);
      if (len < EPS) {
        error(__FILE__,
              __LINE__,
              "edge too small: t = {}, k = {}, edge = {}, Length = {}",
              t,
              k,
              a_edge,
              len);
        return false;
      }
      edg = (1. / len) * edg;
      /* 360deg angle (used for the angle rotation) */
      double cx = Dot(edg, e_x);
      double cy = Dot(edg, e_y);
      alpha[2 * s + k] = atan2(cy, cx);
      if (alpha[2 * s + k] < 0.) alpha[2 * s + k] += 2. * M_PI;
      /* 180deg edge-edge angle (used for the Crouzeix Raviart) */
      double agl = 0.;
      if (a_edge[0] == v1) {
        agl = AngleVectorsAlreadyNormalized(edg, e_x);
      } else if (a_edge[1] == v1) {
        agl = AngleVectorsAlreadyNormalized(edg, -1. * e_x);
      } else if (a_edge[0] == v2) {
        agl = AngleVectorsAlreadyNormalized(-1. * edg, e_x);
      } else if (a_edge[1] == v2) {
        agl = AngleVectorsAlreadyNormalized(-1. * edg, -1. * e_x);
      } else {
        error(__FILE__, __LINE__, "should not happen");
        return false;
      }
      CR_weight[2 * s + k] = -2. / tan(agl);
    }
  }
  /* compute coefficients */
  double
      i_sum = -1. * (CR_weight[0] + CR_weight[1] + CR_weight[2] + CR_weight[3]);
  for (auto &v : CR_weight)v /= i_sum;
  iv.clear();
  ijv.clear();
  ID x_i = 2 * e;
  ID y_i = 2 * e + 1;
  iv.push_back({x_i, 1.});
  iv.push_back({y_i, 1.});
  auto Nd = double(N);
  for (int j = 0; j < 4; ++j) {
    ID x_j = 2 * b_vars[j];
    ID y_j = 2 * b_vars[j] + 1;
    ijv.push_back({{x_i, x_j}, CR_weight[j] * cos(Nd * alpha[j])});
    ijv.push_back({{x_i, y_j}, CR_weight[j] * -1. * sin(Nd * alpha[j])});
    ijv.push_back({{y_i, x_j}, CR_weight[j] * sin(Nd * alpha[j])});
    ijv.push_back({{y_i, y_j}, CR_weight[j] * cos(Nd * alpha[j])});
  }
  return true;
}

bool prepare_system(
    const vector<IV> &K_diag, const vector<IJV> &k_coefficients,
    vector<vector<size_t>> &columns, vector<vector<double>> &values) {
  vector<IJV> coefficients = k_coefficients;
  coefficients.resize(coefficients.size() + K_diag.size());
  for (int i = 0; i < (int)K_diag.size(); ++i)
    coefficients[k_coefficients.size() + i] =
        {{ID(K_diag[i].i), ID(K_diag[i].i)}, K_diag[i].val};
  std::sort(coefficients.begin(), coefficients.end());
  size_t cur_i = coefficients[0].ij[0];
  size_t cur_j = coefficients[0].ij[1];
  double acc_val = coefficients[0].val;
  for (size_t k = 1; k < coefficients.size(); ++k) {
    ID i = coefficients[k].ij[0];
    ID j = coefficients[k].ij[1];
    double v = coefficients[k].val;
    if (i != cur_i || j != cur_j) {
      if (std::abs(acc_val) > EPS) {
        columns[cur_i].push_back(cur_j);
        values[cur_i].push_back(acc_val);
      }
      cur_i = i;
      acc_val = v;
      cur_j = j;
    } else {
      acc_val += v;
    }
  }
  if (std::abs(acc_val) > EPS) {
    columns[cur_i].push_back(cur_j);
    values[cur_i].push_back(acc_val);
  }
  return true;
}

int local_frame(const vector<Vec3> &points,
                const vector<ID3> &triangles,
                const vector<Vec3> &triangle_normal,
                const ID v,
                const std::vector<ID> &tris,
                Vec3 &e_x,
                Vec3 &e_y,
                Vec3 &e_z) {
  e_x = Vec3{DBL_MAX, DBL_MAX, DBL_MAX};
  Vec3 avg{0., 0., 0.};
  for (auto t : tris) {
    avg = avg + triangle_normal[t];
    if (e_x[0] == DBL_MAX) {
      for (size_t le = 0; le < 3; ++le) {
        ID v1 = triangles[t][le];
        ID v2 = triangles[t][(le + 1) % 3];
        if (v1 == v) {
          e_x = points[v2] - points[v1];
          break;
        } else if (v2 == v) {
          e_x = points[v1] - points[v2];
          break;
        }
      }
    }
  }
  if (Length(avg) < 1.e-16) {
    error(__FILE__,
          __LINE__,
          "local frame for {} triangles: avg normal = {}, {}, {}",
          tris.size(),
          avg[0],
          avg[1],
          avg[2]);
    return -1;
  }
  if (Length(e_x) < 1.e-16) {
    error(__FILE__,
          __LINE__,
          "local frame for {} triangles: e_x = {}, {}, {}",
          tris.size(),
          e_x[0],
          e_x[1],
          e_x[2]);
    return -1;
  }
  Normalize(e_x);
  Normalize(avg);
  e_z = avg;
  e_y = Cross(e_z, e_x);
  if (Length(e_y) < 1.e-16) {
    error(__FILE__,
          __LINE__,
          "local frame for {} triangles: e_y = {}, {}, {}",
          tris.size(),
          e_y[0],
          e_y[1],
          e_y[2]);
    return -1;
  }
  Normalize(e_y);
  return 0;
}

int ordered_fan_around_vertex(const vector<ID3> &triangles,
                              const ID v,
                              const std::vector<ID> &tris,
                              std::vector<ID> &vPairs_) {
  std::vector<ID2> vPairs;
  vPairs_.clear();
  if (tris.size() < 3) return -1;
  std::vector<bool> done(tris.size(), false);
  { /* Init */
    auto t = tris[0];
    for (int le = 0; le < 3; ++le) {
      ID v1 = triangles[t][le];
      ID v2 = triangles[t][(le + 1) % 3];
      if (v1 != v && v2 != v) {
        vPairs_.push_back(t * 3 + le);
        vPairs.push_back({v1, v2});
        done[0] = true;
        break;
      }
    }
  }
  while (vPairs.back().back() != vPairs[0][0]) {
    ID last = vPairs.back().back();
    /* Look for next edge connected to last */
    bool found = false;
    for (int t = 0; t < (int)tris.size(); ++t) {
      if (done[t]) continue;
      for (int le = 0; le < 3; ++le) {
        ID v1 = triangles[tris[t]][le];
        ID v2 = triangles[tris[t]][(le + 1) % 3];
        if (v1 != v && v2 != v) {
          if (v1 == last) {
            vPairs.push_back({v1, v2});
            vPairs_.push_back(tris[t] * 3 + le);
            done[t] = true;
            found = true;
            break;
          } else if (v2 == last) {
            vPairs.push_back({v2, v1});
            vPairs_.push_back(tris[t] * 3 + le);
            done[t] = true;
            found = true;
            break;
          }
        }
      }
      if (found) break;
    }
    if (!found) {
      error(__FILE__, __LINE__, "failed to order fan vertex pairs");
      return -1;
    }
  }
  return 0;
}

bool buildUniqueEdges(const std::vector<ID3> &triangles,
                      std::vector<SID3> &triangle_neighbors,
                      std::vector<std::vector<ID>> &nm_triangle_neighbors,
                      std::vector<ID2> &uIEdges,
                      std::vector<ID> &old2IEdge,
                      std::vector<vector<ID>> &uIEdgeToOld,
                      std::vector<std::vector<ID>> &vertex2element,
                      int N) {
  triangle_neighbors.clear();
  triangle_neighbors.resize(triangles.size(), {NO_ID, NO_ID, NO_ID});
  nm_triangle_neighbors.clear();
  vertex2element.clear();
  vertex2element.resize(N);
  constexpr size_t NBF = 3;
  /* Store element 'faces', with duplicates for further 'equality test' */
  std::vector<ID2> faces;
  for (size_t i = 0; i < triangles.size(); ++i) {
    for (size_t lf = 0; lf < NBF; ++lf) {
      ID2 face = Sorted(triangles[i][lf], triangles[i][(lf + 1) % NBF]);
      faces.push_back(face);
      vertex2element[triangles[i][lf]].push_back(i);
    }
  }
  /* Reduce 'duplicated faces' to 'unique faces', keeping the 'old2new' mapping */
  size_t nbUniques = sort_unique_with_perm(faces, uIEdges, old2IEdge);
  /* Build the 'unique face to elements' mapping */
  uIEdgeToOld.resize(nbUniques);
  for (int i = 0; i < (int)triangles.size(); ++i) {
    for (int lf = 0; lf < (int)NBF; ++lf) {
      ID facePos = NBF * i + lf;
      uIEdgeToOld[old2IEdge[facePos]].push_back(facePos);
    }
  }
  /* Loop over 'unique face to elements' and set the element adjacent triangles */
  constexpr ID2 NO_FACE = {NO_ID, NO_ID};
  for (int i = 0; i < (int)nbUniques; ++i) {
    if (uIEdges[i] == NO_FACE) continue;
    if (uIEdges[i][0] == NO_ID) return false;
    if (uIEdges[i][1] == NO_ID) return false;
    if (uIEdgeToOld[i].size() == 1) { /* boundary */
      ID eltId = uIEdgeToOld[i][0] / NBF;
      ID lf = uIEdgeToOld[i][0] % NBF;
      triangle_neighbors[eltId][lf] = NO_ID;
    } else if (uIEdgeToOld[i].size() == 2) { /* regular face */
      ID e1 = uIEdgeToOld[i][0] / NBF;
      ID lf1 = uIEdgeToOld[i][0] % NBF;
      ID e2 = uIEdgeToOld[i][1] / NBF;
      ID lf2 = uIEdgeToOld[i][1] % NBF;
      triangle_neighbors[e1][lf1] = (long long)NBF * e2 + lf2;
      triangle_neighbors[e2][lf2] = (long long)NBF * e1 + lf1;
    } else if (uIEdgeToOld[i].size() > 2) { /* non_manifold face */
      for (int j = 0; j < (int)uIEdgeToOld[i].size(); ++j) {
        ID e = uIEdgeToOld[i][j] / NBF;
        ID lf = uIEdgeToOld[i][j] % NBF;
        std::vector<ID> neighs;
        for (int k = 0; k < (int)uIEdgeToOld[i].size(); ++k)
          if (uIEdgeToOld[i][k] != uIEdgeToOld[i][j])
            neighs.push_back(uIEdgeToOld[i][k]);
        neighs.shrink_to_fit();
        ID pos = nm_triangle_neighbors.size();
        nm_triangle_neighbors.push_back(neighs);
        triangle_neighbors[e][lf] = -((sid)pos + 1);
      }
    }
  }
  /* Reduce memory consumption */
  triangle_neighbors.shrink_to_fit();
  nm_triangle_neighbors.shrink_to_fit();
  vertex2element.shrink_to_fit();
  return true;
}

bool SplitElement(const std::vector<ID> &elements,
                  std::vector<std::vector<ID>> &patches,
                  const std::vector<bool> &lines,
                  const vector<ID> &old2IEdge,
                  const std::vector<SID3> &triangle_neighbors) {
  vector<bool> done(elements.size(), false);
  constexpr size_t NBF = 3;
  for (int i = 0; i < (int)elements.size(); ++i)
    if (!done[i]) {
      vector<ID> patch_elements;
      ID curElement = elements[i];
      std::queue<ID> q;
      q.push(curElement);
      done[curElement] = true;
      while (!q.empty()) {
        curElement = q.front();
        q.pop();
        patch_elements.push_back(curElement);
        for (int le = 0; le < (int)NBF; ++le) {
          ID cur_line = old2IEdge[curElement * NBF + le];
          if (lines[cur_line])continue;
          sid nextElement = triangle_neighbors[curElement][le];
          ID nf = nextElement / 3;
          if (!done[nf]) {
            q.push(nf);
            done[nf] = true;
          }
        }
      }
      patches.push_back(patch_elements);
    }
  return true;
}

bool buildNewPointsAndTriangles(const std::vector<Vec3> &points,
                                const std::vector<ID3> &triangles,
                                const std::vector<ID2> &uIEdges,
                                const std::vector<ID> &to_deal_triangle,
                                const std::vector<ID> &to_deal_line,
                                std::vector<Vec3> &d_points,
                                std::vector<ID2> &d_lines,
                                std::vector<ID3> &d_triangles,
                                std::vector<ID> &old2new,
                                std::vector<ID> &new2old) {
  constexpr ID X = std::numeric_limits<ID>::max();
  old2new.resize(points.size(), X);
  new2old.reserve(3 * triangles.size());
  d_points.reserve(3 * triangles.size());
  d_lines.reserve(to_deal_line.size());
  d_triangles.reserve(triangles.size());
  for (auto l : to_deal_line) {
    auto line = uIEdges[l];
    for (int lv = 0; lv < 2; ++lv) {
      auto v = line[lv];
      if (v >= old2new.size()) old2new.resize(v + 1, X);
      if (old2new[v] == X) {
        old2new[v] = d_points.size();
        d_points.push_back(points[v]);
        new2old.push_back(old2new[v]);
      }
      line[lv] = old2new[v];
    }
    std::sort(line.begin(), line.end());
    d_lines.push_back(line);
  }
  for (auto t : to_deal_triangle) {
    auto tri = triangles[t];
    for (int lv = 0; lv < 3; ++lv) {
      auto v = tri[lv];
      if (v >= old2new.size()) old2new.resize(v + 1, X);
      if (old2new[v] == X) {
        old2new[v] = d_points.size();
        d_points.push_back(points[v]);
        new2old.push_back(old2new[v]);
      }
      tri[lv] = old2new[v];
    }
    d_triangles.push_back(tri);
  }
  return true;
}

int computeCrossFieldWithHeatEquation(int N,
                                      const std::vector<Vec3> &points,
                                      const std::vector<ID3> &triangles,
                                      const std::vector<ID2> &lines,
                                      const std::vector<Vec3> &triangle_normal,
                                      std::vector<Vec3> &triEdgeTheta,
                                      int nbDiffusionLevels = 5,
                                      double thresholdNormConvergence = 1.e-3,
                                      int verbosity = 1) {
  if (N != 4 && N != 6) return false;
  /* Build unique edges and association with adjacent triangles,
   * including the non-manifold cases */
  vector<ID2> uIEdges;
  vector<ID> old2IEdge;
  vector<vector<ID>> uIEdgeToOld;
  std::vector<SID3> triangle_neighbors;
  std::vector<std::vector<ID>> nm_triangle_neighbors;
  std::vector<std::vector<ID>> vertex2element;
  bool oka = buildUniqueEdges(triangles,
                              triangle_neighbors,
                              nm_triangle_neighbors,
                              uIEdges,
                              old2IEdge,
                              uIEdgeToOld,
                              vertex2element,
                              (int)points.size());
  if (!oka) {
    error(__FILE__, __LINE__, "failed to compute mesh adjacent elements");
    return false;
  }
  /* Bbox diagonal is used later to specify the diffusion Length */
  double diag = BboxDiag(points);
  if (verbosity > 0)
    info(__FILE__,
         __LINE__,
         "Heat-based Cross field computation, input: {} points, {} lines, {} triangles, {} internal edges, bbox diag = {}",
         points.size(),
         lines.size(),
         triangles.size(),
         uIEdges.size(),
         diag);
  if (uIEdges.empty()) {
    error(__FILE__, __LINE__, "no internal edges");
    return false;
  }
  /* System unknowns: cos(Nt),sin(Nt) for each edge */
  vector<double> x(2 * uIEdges.size(), 0.);
  /* Initial Dirichlet boundary conditions
   * alignment of crosses with edges (relative angle = 0)
   * theta_e = 0 => (cos4t/sin4t) = (1,0) */
  size_t nbc = 0;
  vector<bool> dirichletEdge(uIEdges.size(), false);
  vector<Vec2> dirichletValue(uIEdges.size(), {1., 0.});
  for (int e = 0; e < (int)uIEdges.size(); ++e)
    if (uIEdgeToOld[e].size() != 2) {
      dirichletEdge[e] = true;
      dirichletValue[e] = {1., 0.};
      nbc += 1;
    }
  for (auto l : lines) {
    /* mark the lines as boundary conditions */
    ID2 edge = Sorted(l[0], l[1]);
    auto it = std::find(uIEdges.begin(), uIEdges.end(), edge);
    if (it != uIEdges.end()) {
      ID e = ID(it - uIEdges.begin());
      dirichletEdge[e] = true;
      dirichletValue[e] = {1., 0.};
      nbc += 1;
    }
  }
  if (verbosity >= 2)
    info(__FILE__,
         __LINE__,
         "- boundary conditions: {} crosses fixed on edges",
         nbc);
  if (verbosity >= 2)
    info(__FILE__,
         __LINE__,
         "- compute stiffness matrix coefficients (Crouzeix-Raviart) ...");
  vector<IV> K_diag;
  vector<IJV> k_coefficients;
  vector<double> rhs(2 * uIEdges.size(), 0.);
  vector<double>
      Mass(2 * uIEdges.size(), 1.); /* diagonal for Crouzeix-Raviart */
  for (int e = 0; e < (int)uIEdges.size(); ++e) {
    if (!dirichletEdge[e]) {
      vector<IV> iv;
      vector<IJV> ijv;
      bool oks = stiffness_coefficient(N,
                                       points,
                                       triangles,
                                       e,
                                       uIEdges,
                                       old2IEdge,
                                       uIEdgeToOld,
                                       triangle_normal,
                                       iv,
                                       ijv);
      if (!oks) {
        error(__FILE__,
              __LINE__,
              "failed to compute stiffness matrix coefficients for e = {}",
              e);
        return false;
      }
      for (auto i : iv)K_diag.push_back(i);
      for (auto i : ijv)k_coefficients.push_back(i);
      ID t1 = uIEdgeToOld[e][0] / 3;
      ID t2 = uIEdgeToOld[e][1] / 3;
      double area1 = TriangleArea(points[triangles[t1][0]],
                                  points[triangles[t1][1]],
                                  points[triangles[t1][2]]);
      double area2 = TriangleArea(points[triangles[t2][0]],
                                  points[triangles[t2][1]],
                                  points[triangles[t2][2]]);
      Mass[2 * e + 0] = 1. / 3 * (area1 + area2);
      Mass[2 * e + 1] = 1. / 3 * (area1 + area2);
    } else { /* Dirichlet BC */
      K_diag.push_back({ID(2 * e + 0), 1.});
      rhs[2 * e + 0] = dirichletValue[e][0];
      K_diag.push_back({ID(2 * e + 1), 1.});
      rhs[2 * e + 1] = dirichletValue[e][1];
    }
  }
  double e_avg = 0.;
  auto e_min = DBL_MAX;
  auto e_max = -DBL_MAX;
  for (auto e : uIEdges) {
    double len = Length(points[e[1]] - points[e[0]]);
    e_avg += len;
    e_min = std::min(len, e_min);
    e_max = std::max(len, e_max);
  }
  e_avg /= (double)uIEdges.size();
  if (verbosity >= 2)
    info(__FILE__,
         __LINE__,
         "- edge size: min={}, avg={}, max={} | bbox diag: {}",
         e_min,
         e_avg,
         e_max,
         diag);
  /* prepare system */
  vector<vector<size_t>> K_columns(2 * uIEdges.size());
  vector<vector<double>> K_values(2 * uIEdges.size());
  {
    bool okp = prepare_system(K_diag, k_coefficients, K_columns, K_values);
    if (!okp) {
      error(__FILE__, __LINE__, "failed to prepare system");
      return false;
    }
  }
  double dtInitial = (0.1 * diag) * (0.1 * diag);
  double dtFinal = (3. * e_min) * (3. * e_min);
  if (verbosity >= 1)
    info(__FILE__,
         __LINE__,
         "- diffusion and projection loop ({} levels, {} unknowns, dt = {} .. {})",
         nbDiffusionLevels,
         2 * uIEdges.size(),
         dtInitial,
         dtFinal);
  double wti = Cpu();
  for (int e = 0; e < (int)uIEdges.size(); ++e) {
    x[2 * e + 0] = dirichletValue[e][0];
    x[2 * e + 1] = dirichletValue[e][1];
  }
  vector<double> steps;
  if (nbDiffusionLevels > 1) {
    for (int i = 0; i < nbDiffusionLevels; ++i) {/* resolution transition */
      double dt = dtInitial
          + (dtFinal - dtInitial) * double(i) / double(nbDiffusionLevels - 1);
      steps.push_back(dt);
    }
  } else {
    steps.push_back(dtFinal);
  }
  {
    /* Initialize system (I/dt + M^-1 * L) x_(i+1) = 1/dt * x_i     (Ax = B) */
    vector<vector<size_t>> a_col = K_columns;
    vector<vector<double>> a_val_add = K_values; /* to get sparsity pattern */
    for (auto &i : a_val_add)std::fill(i.begin(), i.end(), 0.);
    vector<double> B = x;
    vector<double> norms(uIEdges.size(), 0.);
    vector<double> prevNorms = norms;
    auto linear_system = new LinearSystem();
    linear_system->Allocate(2 * (int)uIEdges.size());
    vector<double> diag_sum(2 * uIEdges.size(), 0.);
    for (int i = 0; i < (int)a_col.size(); ++i) {
      if (!dirichletEdge[i / 2]) {
        for (int j = 0; j < (int)a_col[i].size(); ++j) {
          a_val_add[i][j] = 1. / Mass[i] * K_values[i][j];
        }
      } else {
        a_val_add[i] = {1.};
      }
    }
    linear_system->AddSparseCoefficients(a_col, a_val_add, true);
    for (auto &i : a_val_add)i.clear();
    for (auto &i : a_col)i.clear();
    bool okp = linear_system->PreprocessSparsityPattern();
    if (!okp) {
      error(__FILE__, __LINE__, "linear solver analysis failed");
      return false;
    }
    /* Loop over the changing time steps */
    auto prev_dt = DBL_MAX;
    for (int iter = 0; iter < (int)steps.size(); ++iter) {
      if (iter > 0 && steps[iter] > prev_dt) continue;
      double dt = steps[iter];
      prev_dt = dt;
      if (verbosity >= 1)
        info(__FILE__,
             __LINE__,
             "  - step {}/{} | dt = {}, linear system loop ...",
             iter + 1,
             steps.size(),
             dt);
      /* Update LHS matrix with the new time step */
      for (int i = 0; i < (int)a_col.size(); ++i)
        if (!dirichletEdge[i / 2]) {
          if (!dirichletEdge[i / 2]) {
            a_col[i] = {(ID)i};
            a_val_add[i] = {-diag_sum[i] + 1. / dt};
            diag_sum[i] = 1. / dt;
          }
        }
      bool oku = linear_system->AddSparseCoefficients(a_col, a_val_add);
      if (!oku) {
        error(__FILE__, __LINE__, "failed to update linear system");
        return false;
      }
      linear_system->Factorize();
      /* Loop at fixed time step */
      constexpr int sub_iter_max = 25;
      for (int sub_iter = 0; sub_iter < sub_iter_max; ++sub_iter) {
        prevNorms = norms;
        /* Update RHS */
        B = x;
        for (int i = 0; i < (int)B.size(); ++i)
          if (!dirichletEdge[i / 2])
            B[i] /= dt;
        linear_system->AddToRightHandSide(B);
        bool oks = linear_system->Solve(x);
        if (!oks) {
          error(__FILE__, __LINE__, "failed to solve linear system");
          return false;
        }
        /* Normalize Cross field and gather norm stats */
        auto nmi = DBL_MAX;
        auto nma = -DBL_MAX;
        for (int e = 0; e < (int)uIEdges.size(); ++e) {
          norms[e] =
              std::sqrt(x[2 * e] * x[2 * e] + x[2 * e + 1] * x[2 * e + 1]);
          nmi = std::min(nmi, norms[e]);
          nma = std::max(nma, norms[e]);
          if (!dirichletEdge[e] && norms[e] > EPS) {
            x[2 * e + 0] /= norms[e];
            x[2 * e + 1] /= norms[e];
          }
        }
        const double EPS_NORM = 1.e-1;
        if (nma > 1. + EPS_NORM) {
          steps[iter] /= 10;
          dt = steps[iter];
          iter -= 1;
          if (verbosity >= 2)
            warn(__FILE__,
                 __LINE__,
                 "    |   max(norm)={} (should be 1.), solve failed, new time step: dt = {}",
                 nma,
                 dt);
          break;
        }
        if (sub_iter > 0 || iter > 0) {
          double lin_f = 0.;
          for (size_t i = 0; i < norms.size(); ++i)
            if (!dirichletEdge[i]) {
              lin_f = std::max(lin_f, norms[i] - prevNorms[i]);
            }
          if (verbosity >= 3)
            info(__FILE__,
                 __LINE__,
                 "   |   system solved, norm diff max: {}, norm range: {} - {}",
                 lin_f,
                 nmi,
                 nma);
          if (lin_f < 1.e-3) break;
        } else {
          if (verbosity >= 3)
            info(__FILE__, __LINE__, "           |   system solved");
        }
      }
    }
    delete linear_system;
  }
  double et = Cpu() - wti;
  if (verbosity >= 2)
    info(__FILE__, __LINE__, "- Cross field elapsed time: {}", et);
  { /* Export solution */
    triEdgeTheta.resize(triangles.size());
    for (int t = 0; t < (int)triangles.size(); ++t) {
      for (int le = 0; le < 3; ++le) {
        ID e = old2IEdge[3 * t + le];
        double
            len = std::sqrt(x[2 * e] * x[2 * e] + x[2 * e + 1] * x[2 * e + 1]);
        double theta = (len > EPS) ? 1. / double(N)
            * atan2(x[2 * e + 1] / len, x[2 * e] / len) : 0.;
        triEdgeTheta[t][le] = theta;
      }
    }
  }
  return true;
}

inline Vec3 crouzeix_raviart_interpolation(Vec3 lambda, Vec3 edge_values[3]) {
  double x = lambda[1];
  double y = lambda[2];
  Vec3 shape = {1.0 - 2.0 * y, -1.0 + 2.0 * (x + y), 1.0 - 2.0 * x};
  return shape[0] * edge_values[0] + shape[1] * edge_values[1]
      + shape[2] * edge_values[2];
}

}

int convertToPerTriangleCrossFieldDirections(
    int Ns,
    const std::vector<Vec3> &points,
    const std::vector<ID2> &uIEdges,
    const std::vector<ID3> &triangles,
    const std::vector<Vec3> &edgeNormal,
    const std::vector<double> &triEdgeTheta,
    std::vector<Vec3> &edgeDirections) {
  edgeDirections.resize(triEdgeTheta.size());
  for (int l = 0; l < (int)triEdgeTheta.size(); ++l) {
    Vec3 N = edgeNormal[l];
    double A = triEdgeTheta[l];
    Vec3 tgt = points[uIEdges[l][1]] - points[uIEdges[l][0]];
    if (Length(tgt) < 1.e-16) {
      warn(__FILE__, __LINE__, "Edge (le={}), Length = {}", l, Length(tgt));
      if (Length(tgt) == 0.) { return 0; }
    }
    Normalize(tgt);
    Vec3 tgt2 = Cross(N, tgt);
    Normalize(tgt2);
    Vec3 cross1 = tgt * cos(A) + tgt2 * sin(A);
    edgeDirections[l] = cross1;
  }
  return 1;
}

int convertToPerTriangleCrossFieldDirections(
    int Ns,
    const std::vector<Vec3> &points,
    const std::vector<ID3> &triangles,
    const std::vector<Vec3> &triangles_normal,
    const std::vector<std::array<double, 3> > &triEdgeTheta,
    std::vector<std::array<double, 9> > &triangleDirections) {

  triangleDirections.resize(triangles.size());

  for (size_t f = 0; f < triangles.size(); ++f) {
    auto t = triangles[f];

    Vec3 N = triangles_normal[f];

    /* Compute one branch at triangle vertices */
    Vec3 refCross = {0., 0., 0.};
    Vec3 avgCross = {0., 0., 0.};
    Vec3 lifted_dirs[3];
    for (size_t le = 0; le < 3; ++le) {
      /* Get cross angle */
      ID2 edge = {t[le], t[(le + 1) % 3]};
      if (edge[0] > edge[1]) std::reverse(edge.begin(), edge.end());

      double A = triEdgeTheta[f][le];

      /* Local reference frame */
      Vec3 tgt = points[edge[1]] - points[edge[0]];

      if (Length(tgt) < 1.e-16) {
        warn(__FILE__, __LINE__, "Edge (tri={},le={}), length = {}", f, le,
             Length(tgt));
        if (Length(tgt) == 0.) { return -1; }
      }
      Normalize(tgt);

      Vec3 tgt2 = Cross(N, tgt);
      Normalize(tgt2);

      Vec3 cross1 = tgt * cos(A) + tgt2 * sin(A);
      Normalize(cross1);

      Vec3 cross2 = Cross(N, cross1);
      Normalize(cross2);

      if (le == 0) {
        refCross = cross1;
        avgCross = avgCross + cross1;
        lifted_dirs[le] = refCross;
      } else {
        Vec3 closest = {0., 0., 0.};
        double dotmax = -DBL_MAX;
        for (int k = 0; k < Ns; ++k) {
          double agl = A + double(k) / double(Ns) * 2. * M_PI;
          Vec3 candidate = tgt * cos(agl) + tgt2 * sin(agl);
          Normalize(candidate);

          if (Dot(candidate, refCross) > dotmax) {
            closest = candidate;
            dotmax = Dot(candidate, refCross);
          }
        }

        lifted_dirs[le] = closest;
        avgCross = avgCross + closest;
      }
    }
    Vec3 vertex_dirs[3];
    std::array<double, 9> tDirs;
    for (size_t lv = 0; lv < 3; ++lv) {
      Vec3 lambda = {0., 0., 0.};
      lambda[lv] = 1.;
      Vec3 dir = crouzeix_raviart_interpolation(lambda, lifted_dirs);
      if (Length2(dir) != 0.) Normalize(dir);
      for (size_t d = 0; d < 3; ++d) {
        tDirs[3 * lv + d] = dir.data()[d];
      }
    }
    triangleDirections[f] = tDirs;
  }

  return 1;
}

int detectCrossFieldSingularities(int Ns,
                                  const std::vector<Vec3> &points,
                                  const std::vector<ID3> &triangles,
                                  const std::vector<ID2> &UIEdges,
                                  const std::vector<std::vector<ID>> &v2elements,
                                  const std::vector<Vec3> &triangle_normal,
                                  const std::vector<bool> &isBorderV,
                                  const std::vector<ID> &old2UIEdge,
                                  const std::vector<double> &triEdgeTheta,
                                  std::vector<ID3> &singularities) {
  singularities.clear();
  if (triEdgeTheta.size() != UIEdges.size()) {
    error(__FILE__, __LINE__, "detect Cross field singularities: wrong inputs");
    return 0;
  }
  /* Compute index at each vertex (looking on the one rings) */
  vector<int> vertexIndex(points.size(), 0);
  for (int i = 0; i < (int)v2elements.size(); ++i)
    if (!isBorderV[i]) {
      ID v = i;
      const vector<ID> &tris = v2elements[i];
      Vec3 e_x, e_y, N;
      int st =
          local_frame(points, triangles, triangle_normal, v, tris, e_x, e_y, N);
      if (st != 0) {
        error(__FILE__, __LINE__, "no local frame for v={}", v);
        continue;
      }
      std::vector<ID> vPairs;
      int st2 = ordered_fan_around_vertex(triangles, v, tris, vPairs);
      if (st2 != 0) {
        error(__FILE__, __LINE__, "no ordered fan around v={}", v);
        continue;
      }
      std::vector<Vec3> rep_vectors(vPairs
                                        .size()); /* Cross field representation vector in local frame (+e_x,+e_y) */
      for (int j = 0; j < (int)vPairs.size(); ++j) {
        ID v1 = UIEdges[old2UIEdge[vPairs[j]]][0];
        ID v2 = UIEdges[old2UIEdge[vPairs[j]]][1];
        Vec3 tgt;
        tgt = points[v2] - points[v1];
        Normalize(tgt);
        Vec3 tgt2 = Cross(N, tgt);
        Normalize(tgt2);
        if (old2UIEdge[vPairs[j]] >= triEdgeTheta.size()) {
          error(__FILE__,
                __LINE__,
                "Edge ({},{}) not found in edgeTheta",
                vPairs[j],
                old2UIEdge[vPairs[j]]);
          return 0;
        }
        double A = triEdgeTheta[old2UIEdge[vPairs[j]]];
        Vec3 branchInQuadrant;
        bool found = false;
        double dp_max = -DBL_MAX;
        auto k_max = (size_t)-1;
        for (int k = 0; k < Ns; ++k) {
          double agl = A + double(k) / double(Ns) * 2. * M_PI;
          Vec3 branch = tgt * cos(agl) + tgt2 * sin(agl);
          if (Dot(branch, e_x) >= 0. && Dot(branch, e_y) >= 0.) {
            branchInQuadrant = branch;
            found = true;
            break;
          }
          double dp = Dot(branch, Vec3{std::sqrt(2), std::sqrt(2), 0.});
          if (dp > dp_max) {
            dp_max = dp;
            k_max = k;
          }
        }
        if (!found && k_max != (size_t)-1) { /* if numerical errors */
          double agl = A + double(k_max) / double(Ns) * 2. * M_PI;
          branchInQuadrant = tgt * cos(agl) + tgt2 * sin(agl);
        }
        double theta = acos(Dot(branchInQuadrant, e_x));
        rep_vectors[j] = {cos(double(Ns) * theta), sin(double(Ns) * theta), 0.};
      }
      double sum_diff = 0.;
      for (int j = 0; j < (int)vPairs.size(); ++j) {
        Vec3 vertical{0., 0., 1.};
        const Vec3 r1 = rep_vectors[j];
        const Vec3 r2 = rep_vectors[(j + 1) % rep_vectors.size()];
        double diff_angle = Angle2Pi(r1, r2, vertical);
        if (diff_angle > M_PI) diff_angle -= 2. * M_PI;
        sum_diff += diff_angle;
      }
      if (std::abs(sum_diff - 2 * M_PI) < 0.05 * M_PI) {
        vertexIndex[v] = -1;
      } else if (std::abs(sum_diff + 2 * M_PI) < 0.05 * M_PI) {
        vertexIndex[v] = 1;
      } else if (std::abs(sum_diff) > 0.05 * M_PI) {
//        debug(__FILE__,
//              __LINE__,
//              "singularity detection, v={}, sum of diff is {} deg, ignored",
//              v,
//              sum_diff);
      }
    }
  /* triangle singularities */
  for (int t = 0; t < (int)triangles.size(); ++t) {
    ID v1 = triangles[t][0];
    ID v2 = triangles[t][1];
    ID v3 = triangles[t][2];
    for (double index : {-1., 1.}) {
      if (vertexIndex[v1] == index
          && vertexIndex[v2] == index
          && vertexIndex[v3] == index) {
        singularities.push_back({NO_ID, NO_ID, (ID)t});
        /* remove from list of available singularities */
        vertexIndex[v1] = 0.;
        vertexIndex[v2] = 0.;
        vertexIndex[v3] = 0.;
      }
    }
  }
  /* edge singularities */
  for (int t = 0; t < (int)triangles.size(); ++t) {
    for (int le = 0; le < 3; ++le) {
      ID v1 = triangles[t][le];
      ID v2 = triangles[t][(le + 1) % 3];
      for (double index : {-1., 1.}) {
        if (vertexIndex[v1] == index
            && vertexIndex[v2] == index) {
          singularities.push_back({NO_ID, (ID)t * 3 + (ID)le, NO_ID});
          /* remove from list of available singularities */
          vertexIndex[v1] = 0.;
          vertexIndex[v2] = 0.;
        }
      }
    }
  }
  /* vertex singularities */
  for (int kv = 0; kv < (int)vertexIndex.size(); ++kv)
    if (vertexIndex[kv] != 0) {
      singularities.push_back({(ID)kv, NO_ID, NO_ID});
    }
  info(__FILE__, __LINE__,
       "detect Cross field singularities: found {} singularities ({} triangles)",
       singularities.size(),
       triangles.size());
  return 1;
}

int BuildBackgroundMeshAndGuidingField(const std::vector<Vec3> &points,
                                       const std::vector<ID3> &triangles,
                                       const std::vector<ID2> &lines,
                                       std::vector<Vec3> &global_edge_directions,
                                       std::vector<std::array<double,9>> &global_triangle_directions,
                                       std::vector<ID3> &global_singularity_list) {
  constexpr int NBF = 3;
  vector<ID2> uIEdges;
  vector<ID> old2IEdge;
  vector<vector<ID>> uIEdgeToOld;
  std::vector<SID3> triangle_neighbors;
  std::vector<std::vector<ID>> nm_triangle_neighbors;
  std::vector<std::vector<ID>> vertex2element;
  bool oka = buildUniqueEdges(triangles,
                              triangle_neighbors,
                              nm_triangle_neighbors,
                              uIEdges,
                              old2IEdge,
                              uIEdgeToOld,
                              vertex2element,
                              (int)points.size());
  if (!oka) error(__FILE__, __LINE__, "build failed");
  vector<bool> done(triangles.size(), false);
  vector<ID> faces(triangles.size());
  std::iota(faces.begin(), faces.end(), 0);
  vector<bool> border_line(uIEdges.size(), false);
  vector<bool> border_v(points.size(), false);
  for (int e = 0; e < (int)uIEdges.size(); ++e)
    if (uIEdgeToOld[e].size() != 2) {
      border_v[uIEdges[e][0]] = true;
      border_v[uIEdges[e][1]] = true;
      border_line[e] = true;
    }
  for (auto &l : lines) {
    ID2 nl = Sorted(l[0], l[1]);
    auto it = std::lower_bound(uIEdges.begin(), uIEdges.end(), nl);
    if (it != uIEdges.end()) {
      ID e = ID(it - uIEdges.begin());
      border_v[uIEdges[e][0]] = true;
      border_v[uIEdges[e][1]] = true;
      border_line[e] = true;
    }
  }
  std::vector<Vec3> triangle_normal(triangles.size());
  for (int t = 0; t < (int)triangles.size(); ++t)
    triangle_normal[t] =
        TriangleNormal(points[triangles[t][0]],
                       points[triangles[t][1]],
                       points[triangles[t][2]]);

  std::vector<Vec3> edgeNormal(uIEdges.size());
  for (int i = 0; i < (int)uIEdges.size(); ++i) {
    edgeNormal[i] = {0., 0., 0.};
    for (int j = 0; j < (int)uIEdgeToOld[i].size(); ++j) {
      ID face = uIEdgeToOld[i][j] / 3;
      edgeNormal[i] = edgeNormal[i] + triangle_normal[face];
    }
    Normalize(edgeNormal[i]);
  }
  vector<double> global_edge_theta(uIEdges.size(), 0);
  vector<Vec3>global_triangle_theta(triangles.size());
  vector<vector<ID>> patches;
  SplitElement(faces, patches, border_line, old2IEdge, triangle_neighbors);
  for (int f = 0; f < (int)patches.size(); ++f) {
    std::vector<Vec3> d_points;
    std::vector<ID2> d_lines;
    std::vector<ID3> d_triangles;
    std::vector<ID> old2new;
    std::vector<ID> new2old;
    vector<ID> to_deal_line;
    to_deal_line.reserve(patches[f].size() * NBF);
    for (auto t : patches[f]) {
      for (int i = 0; i < NBF; ++i)
        if (border_line[old2IEdge[t * NBF + i]])
          to_deal_line.push_back(old2IEdge[t * NBF + i]);
    }
    sort_unique(to_deal_line);
    buildNewPointsAndTriangles(points,
                               triangles,
                               uIEdges,
                               patches[f],
                               to_deal_line,
                               d_points,
                               d_lines,
                               d_triangles,
                               old2new,
                               new2old);
    std::vector<Vec3> patch_normal(d_triangles.size());
    for (int t = 0; t < (int)d_triangles.size(); ++t)
      patch_normal[t] = triangle_normal[patches[f][t]];
    /* Cross field */
    std::vector<Vec3> triEdgeTheta;

    int nbDiffusionLevels = 4;
    double thresholdNormConvergence = 1.e-2;
    int verbosity = 0;
    info(__FILE__,
         __LINE__,
         "- Patch {}/{}: compute Cross field ({} triangles, {} B.C. "
         "edges, {} diffusion levels) ...\n",
         f,
         patches.size(),
         d_triangles.size(),
         d_lines.size(),
         nbDiffusionLevels);
    int scf = computeCrossFieldWithHeatEquation(
        4,
        d_points,
        d_triangles,
        d_lines,
        patch_normal,
        triEdgeTheta,
        nbDiffusionLevels,
        thresholdNormConvergence,
        verbosity);
    if (!scf) {
      warn(__FILE__,
           __LINE__,
           "- Patch {}: failed to compute Cross field\n",
           f);
    }
    for (int t = 0; t < (int)patches[f].size(); ++t) {
      global_triangle_theta[patches[f][t]]=triEdgeTheta[t];
      for (int le = 0; le < 3; ++le) {
        global_edge_theta[old2IEdge[patches[f][t] * 3 + le]] =
            triEdgeTheta[t][le];
      }
    }
  }
  /* Cross field singularities */
  int scsi = detectCrossFieldSingularities(4,
                                           points,
                                           triangles,
                                           uIEdges,
                                           vertex2element,
                                           triangle_normal,
                                           border_v,
                                           old2IEdge,
                                           global_edge_theta,
                                           global_singularity_list);
  if (!scsi)
    error(__FILE__,
          __LINE__,
          "- failed to compute Cross field singularities");

  int sc = convertToPerTriangleCrossFieldDirections(
      4,
      points,
      uIEdges,
      triangles,
      edgeNormal,
      global_edge_theta,
      global_edge_directions);

  sc &= convertToPerTriangleCrossFieldDirections(4,
                                                points,
                                                triangles,
                                                triangle_normal,
                                                global_triangle_theta,
                                                global_triangle_directions);

  //  //the code is present cross field
//  vector<Vec3>ps;
//  vector<ID2>ls;
//  for(int i=0;i<uIEdges.size();++i){
//    Vec3 v1=points[uIEdges[i][0]];
//    Vec3 v2=points[uIEdges[i][1]];
//    double len= Length(v1-v2);
//    Vec3 m=(v1+v2)*0.5;
//    Vec3 v=global_edge_directions[i];
//    Vec3 mv=m+v*len*0.5;
//    ps.emplace_back(m);
//    ps.emplace_back(mv);
//    ls.emplace_back(ID2{2*(ID)i,2*(ID)i+1});
//  }
//  OutputFile("dir.vtk",ps,ls);
  if (!sc)
    error(__FILE__,
          __LINE__,
          "- failed to resample Cross field at triangle corners");

  return 0;
}