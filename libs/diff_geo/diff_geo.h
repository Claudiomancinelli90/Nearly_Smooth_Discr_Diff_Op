#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>
// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

} // namespace yocto
#include <VTP/geodesic_algorithm_exact.h>
#include <VTP/geodesic_mesh.h>
#include <iostream>
#include <stdio.h>
#include <unordered_set>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_shape.h>
using namespace yocto;
enum transport_mode { V2V, V2T, T2V, T2T };
enum simplex { vert, edge };
enum field { exact, graph };

struct shape_geometry {
  vector<vec3i> adjacencies = {};
  vector<vector<int>> v2t = {};
  vector<vector<float>> angles = {};
  vector<float> total_angles = {};
  float avg_edge_length = 0.f;
};

struct shape_op {
  Eigen::SparseMatrix<double> Grad;
  Eigen::SparseMatrix<double> stencil_Grad;
  Eigen::SparseMatrix<double> Lap;
  Eigen::SparseMatrix<double> Hess;
  // vector<Eigen::MatrixXd> CMat;
  vector<Eigen::ColPivHouseholderQR<Eigen::MatrixXd>> CMat;
  vector<Eigen::MatrixXd> Rhs;
  vector<Eigen::MatrixXd> quadrics;
};
void clean_bary(mesh_point &sample);
vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int pid);
vector<vec3f> polyline_pos(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<mesh_point> &poly);

vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p);
float tid_area(const vector<vec3i> &triangles, const vector<vec3f> &positions,
               const int tid);

vector<mesh_point> straightest_path_from_vert_in_triangle(
    const shape_data &data, const shape_geometry &geometry, const int tid,
    const int vid, const float &len, const vec3f &dir);

vector<uint> one_ring(const vector<vec3i> &triangles,
                      const vector<vec3i> &adjacencies,
                      const vector<vector<int>> &v2t, const int vid);
vector<uint> extended_one_ring(const vector<vec3i> &triangles,
                               const vector<vec3i> &adjacencies,
                               const vector<vector<int>> &v2t, const int vid);
vec3f project_vec(const vec3f &v, const vec3f &n);

mesh_point point_from_vert(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int j = 0);
int forced_vert_from_point(const vector<vec3i> &triangles, const mesh_point &p);
vector<float> subdivide_angles(const int number_of_subdivision);

vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle);

mesh_point make_mesh_point(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int tid = -1);

float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vector<int>> &v2t,
                             const vector<vec3f> &normals, const int vid,
                             const vec3f &v);

vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions,
                  const vector<vector<int>> &v2t, const vector<vec3f> &normals,
                  int vid);
vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int tid);

std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions);

std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions,
                       const mesh_point &source);

Eigen::VectorXd wrapper(const vector<float> &f);

float nbr_avg_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, const int vid);
float nbr_min_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, int vid);
float avg_edge_length(const shape_data &mesh, const shape_geometry &topology);

std::tuple<vector<mesh_point>, vector<float>, vector<float>>
weighted_stencil(const vector<vec3i> &triangles, const vector<vec3f> &positions,
                 const geodesic_solver &solver, const vector<vec3f> &normals,
                 const vector<vec3i> &adjacencies,
                 const vector<vector<int>> &v2t,
                 const vector<vector<float>> &angles,
                 const vector<float> &total_angles, const int vid, const int n);
std::tuple<vector<mesh_point>, vector<float>, vector<float>>
uniform_stencil(const shape_data &data, const shape_geometry &geometry,
                const int vid, const int number_of_samples);

vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source);
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const mesh_point &source);
vector<float> heat_distance_field(const shape_data &mesh,
                                  const shape_geometry &topology,
                                  const int sources);

vec3f compute_PCE(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, const vector<float> &f,
                  const mesh_point &p);
Eigen::SparseMatrix<double> PCE_matrix(const vector<vec3i> &triangles,
                                       const vector<vec3f> &positions);

shape_op init_discrete_diff_op_xu(const shape_data &data,
                                  const shape_geometry &geometry,
                                  const dual_geodesic_solver &solver,
                                  const bool grad, const bool lap,
                                  const bool hess);

Eigen::VectorXd compute_hessian(const shape_op &operators,
                                const vector<float> &field);
Eigen::Matrix2d compute_hessian(const Eigen::VectorXd &Hessian, const int vid);

vector<pair<int, float>> cotangent_entries(const shape_data &mesh,
                                           const shape_geometry &topology);
Eigen::SparseMatrix<double>
cot_lap_stiffness_matrix(const shape_data &mesh,
                         const shape_geometry &topology);
Eigen::SparseMatrix<double> mass_matrix(const shape_data &mesh,
                                        const shape_geometry &topology,
                                        const bool lumped = true);
void trace_in_triangles(const vector<vec3f> &positions,
                        const vector<vec3i> &triangles, const vec3f &dir,
                        const vec3f &bary, const int pid, vec3f &sample_pos,
                        vec3f &sample_bary);

int next_tid(const vector<vector<float>> &angles,
             const vector<vec3f> &positions, const vector<vector<int>> &v2t,
             const vector<vec3i> &triangles, const vector<vec3f> &normals,
             const int from, const vec3f &v);
pair<vector<mesh_point>, float> straightest_geodesic(
    const geodesic_solver &solver, const vector<vec3i> &triangles,
    const vector<vec3f> &positions, const vector<vec3f> &normals,
    const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const int max_crossed_tri);
vector<mesh_point> straightest_geodesic(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3f> &normals, const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const bool check_from = true);
void parallel_transp(const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies,
                     const vector<vector<int>> &v2t, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode);

vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const Eigen::VectorXd &f, bool normalized = false);

vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const vector<float> &f, bool normalized = false);

vector<float> compute_laplacian(const shape_op &operators,
                                const vector<float> &distances);

vector<float> cotangent_laplacian(const shape_data &mesh,
                                  const shape_geometry &topology,
                                  const vector<float> &field);

geodesic_solver extended_solver(const shape_data &data,
                                const dual_geodesic_solver &dual_solver,
                                shape_geometry &geometry, const int k);

inline vec3f tri_bary_coords(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid,
                             const vec3f &p) {
  auto px = positions[triangles[tid].x];
  auto py = positions[triangles[tid].y];
  auto pz = positions[triangles[tid].z];
  return tri_bary_coords(px, py, pz, p);
}

inline vec3f tid_normal(const vector<vec3i> &triangles,
                        const vector<vec3f> &positions, const int tid) {
  auto p0 = positions[triangles[tid].x];
  auto p1 = positions[triangles[tid].y];
  auto p2 = positions[triangles[tid].z];

  return normalize(cross(p1 - p0, p2 - p0));
}

inline float tid_area(const shape_data &mesh, const int tid) {
  return tid_area(mesh.triangles, mesh.positions, tid);
}

inline int node_is_adjacent(const geodesic_solver &solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}
template <typename T> inline int find(const T &vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x)
      return i;
  return -1;
}
inline vec3f get_bary(const vec2f &uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}
