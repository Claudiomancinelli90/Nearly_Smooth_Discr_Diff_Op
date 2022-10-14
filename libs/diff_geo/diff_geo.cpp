//
//  karcher.cpp
//  glfw
//
//  Created by Claudio Mancinelli on 21/07/2020.
//

#include "diff_geo.h"
#include <utils/logging.h>
using namespace logging;
using namespace yocto;

#include <deque>

// utility (geometry)
// project a vector v onto a plane having n as normal
vec3f project_vec(const vec3f &v, const vec3f &n) {
  auto proj = n * dot(v, n);

  return v - proj;
}
vec3f rot_vect(const vec3f &p, const vec3f &axis, const float angle) {
  auto M = rotation_frame(axis, angle);
  auto v = p;
  return transform_vector(M, v);
}
vec2f rot_vect(const vec2f &p, const float theta) {
  auto M = mat2f{{yocto::cos(theta), -yocto::sin(theta)},
                 {yocto::sin(theta), yocto::cos(theta)}};
  auto v = M * p;
  return v;
}

vec3f tid_centroid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int tid) {
  vec3f p0 = positions[triangles[tid].x];
  vec3f p1 = positions[triangles[tid].y];
  vec3f p2 = positions[triangles[tid].z];

  return (p0 + p1 + p2) / 3.0;
}
float tid_area(const vector<vec3i> &triangles, const vector<vec3f> &positions,
               const int tid) {
  return triangle_area(positions[triangles[tid].x], positions[triangles[tid].y],
                       positions[triangles[tid].z]);
}

vector<vec3f> polyline_pos(const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<mesh_point> &poly) {
  auto result = vector<vec3f>(poly.size());
  for (auto i = 0; i < poly.size(); ++i) {
    result[i] = eval_position(triangles, positions, poly[i]);
  }
  return result;
}
int forced_vert_from_point(const vector<vec3i> &triangles,
                           const mesh_point &p) {
  auto bary = vector<pair<float, int>>{
      {1 - p.uv.x - p.uv.y, 0}, {p.uv.x, 1}, {p.uv.y, 2}};

  sort(bary.begin(), bary.end());
  return triangles[p.face][bary.back().second];
}

void clean_bary(mesh_point &sample) {

  auto bary = sample.uv;
  auto coords = vector<pair<float, int>>{
      {1 - bary.x - bary.y, 0}, {bary.x, 1}, {bary.y, 2}};
  sort(coords.begin(), coords.end());
  vec3f bary3d = get_bary(bary);
  if (coords[0].first < 0 && coords[1].first > 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[2].second] =
        1 - bary3d[coords[1].second] - bary3d[coords[0].second];
  } else if (coords[0].first < 0) {
    bary3d[coords[0].second] = 0;
    bary3d[coords[1].second] = 0;
    bary3d[coords[2].second] = 1;
  }

  sample = {sample.face, {bary3d.y, bary3d.z}};
}
float nbr_avg_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, const int vid) {
  float len = 0;
  auto s = (int)v2t[vid].size();
  auto vert = positions[vid];
  for (auto &t : v2t[vid]) {
    auto k = find(triangles[t], vid);
    len += length(positions[triangles[t][(k + 1) % 3]] - vert);
  }
  return (s != 0) ? len / s : len;
}
float max_nbr_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, int vid) {
  auto nbr = vector<float>{};
  auto vert = positions[vid];
  for (auto tid : v2t[vid]) {
    auto k = find(triangles[tid], vid);
    nbr.push_back(length(positions[triangles[tid][(k + 1) % 3]] - vert));
  }
  auto lambda = flt_min;
  auto s = (int)nbr.size();
  for (int i = 0; i < s; ++i) {
    lambda = std::max(lambda, nbr[i]);
  }

  return lambda;
}
float nbr_min_edge_length(const vector<vec3i> &triangles,
                          const vector<vec3f> &positions,
                          const vector<vector<int>> &v2t, int vid) {
  auto nbr = vector<float>{};
  auto vert = positions[vid];
  for (auto tid : v2t[vid]) {
    auto k = find(triangles[tid], vid);
    nbr.push_back(length(positions[triangles[tid][(k + 1) % 3]] - vert));
  }
  auto lambda = flt_max;
  auto s = (int)nbr.size();
  for (int i = 0; i < s; ++i) {
    lambda = std::min(lambda, nbr[i]);
  }

  return lambda;
}
float avg_edge_length(const shape_data &data, const shape_geometry &topology) {
  auto avg = 0.f;
  auto count = 0.f;
  for (auto i = 0; i < data.positions.size(); ++i) {
    auto star = topology.v2t[i];
    for (auto tid : star) {
      auto k = find(data.triangles[tid], i);
      auto vid0 = data.triangles[tid][(k + 1) % 3];
      if (vid0 < i)
        continue;

      avg += length(data.positions[vid0] - data.positions[i]);
      ++count;
    }
  }
  return avg / count;
}
vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions,
                  const vector<vector<int>> &v2t, const vector<vec3f> &normals,
                  int vid) {
  auto tid = v2t[vid][0];
  auto k = find(triangles[tid], vid);
  vec3f v = positions[triangles[tid][(k + 1) % 3]] - positions[vid];
  vec3f e = normalize(project_vec(v, normals[vid]));
  return e;
}

vec3f polar_basis(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int tid) {
  auto c = tid_centroid(triangles, positions, tid);
  vec3f v = positions[triangles[tid].x];
  return normalize(v - c);
}

// Compute polar coordinates
float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vector<int>> &v2t,
                             const vector<vec3f> &normals, const int vid,
                             const vec3f &v) {
  auto e = polar_basis(triangles, positions, v2t, normals, vid);

  auto teta = angle(v, e);

  if (dot(cross(e, v), normals[vid]) < 0)
    teta *= -1;

  return teta;
}

float angle_in_tangent_space(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid,
                             const vec3f &v) {
  auto teta = 0.f;

  auto e = polar_basis(triangles, positions, tid);
  auto n = tid_normal(triangles, positions, tid);
  teta = angle(v, e);
  if (dot(cross(e, v), n) < 0)
    teta = 2 * M_PI - teta;

  return teta;
}
mesh_point point_from_vert(const vector<vec3i> &triangles,
                           const vector<vector<int>> &v2t, const int vid,
                           const int j) {
  auto tid = v2t[vid][j];
  auto k = find(triangles[tid], vid);
  auto bary = zero3f;
  bary[k] = 1;
  return {tid, {bary.y, bary.z}};
}
int vert_from_point(const vector<vec3i> &triangles, const mesh_point &p) {
  auto b = vector<pair<float, int>>(3);
  auto bary = get_bary(p.uv);
  for (auto i = 0; i < 3; ++i) {
    b[i] = std::make_pair(bary[i], i);
  }
  sort(b.begin(), b.end());

  return triangles[p.face][b.back().second];
}
vec3f tri_bary_coords(const vec3f &v0, const vec3f &v1, const vec3f &v2,
                      const vec3f &p) {
  vec3f wgts = vec3f{0.0, 0.0, 0.0};
  vec3f u = v1 - v0, v = v2 - v0, w = p - v0;
  float d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
        d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0)
    return zero3f;

  wgts[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(wgts[2]));
  wgts[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(wgts[1]));
  wgts[0] = 1.0 - wgts[1] - wgts[2];
  assert(!isnan(wgts[0]));

  return wgts;
}
vector<int> one_ring(const vector<vec3i> &triangles,
                     const vector<vector<int>> &v2t, const int vid,
                     const bool extended = false) {
  auto ring = vector<int>{};
  for (auto tid : v2t[vid]) {
    auto k = find(triangles[tid], vid);
    ring.push_back(triangles[tid][(k + 1) % 3]);
    if (extended) {
      auto eid =
          vec2i{triangles[tid][(k + 1) % 3], triangles[tid][(k + 2) % 3]};
      ring.push_back(opposite_vertex(triangles[tid], eid));
    }
  }

  return ring;
}
vector<int> k_ring(const vector<vec3i> &triangles,
                   const vector<vector<int>> &v2t, const int vid, const int k) {
  vector<int> ring;
  vector<int> active_set = {vid};
  for (int i = 0; i < k; ++i) {
    vector<int> next_active_set;
    for (int j = 0; j < active_set.size(); ++j) {
      auto nbr = one_ring(triangles, v2t, active_set[j]);
      for (int h = 0; h < nbr.size(); ++h) {
        int curr = nbr[h];
        if (find(ring, curr) == -1 && curr != vid) {
          next_active_set.push_back(curr);
          ring.push_back(curr);
        }
      }
    }
    active_set = next_active_set;
  }

  return ring;
}
int opposite_vert(const vector<vec3i> &triangles,
                  const vector<vec3i> &adjacencies, const int tid,
                  const int vid) {
  auto nei = opposite_face(triangles, adjacencies, tid, vid);
  auto h = find(adjacencies[nei], tid);
  return triangles[nei][(h + 2) % 3];
}
vector<uint> one_ring(const vector<vec3i> &triangles,
                      const vector<vec3i> &adjacencies,
                      const vector<vector<int>> &v2t, const int vid) {
  auto &star = v2t[vid];
  auto nbr = vector<uint>(star.size());
  for (auto i = 0; i < star.size(); ++i) {
    auto tid = star[i];
    auto tr = triangles[tid];
    auto k = find(tr, vid);
    nbr[i] = tr[(k + 1) % 3];
  }
  return nbr;
}
vector<uint> extended_one_ring(const vector<vec3i> &triangles,
                               const vector<vec3i> &adjacencies,
                               const vector<vector<int>> &v2t, const int vid) {
  auto &star = v2t[vid];
  auto nbr = vector<uint>(2 * star.size());
  for (auto i = 0; i < star.size(); ++i) {
    auto tid = star[i];
    auto tr = triangles[tid];
    auto k = find(tr, vid);
    nbr[2 * i] = tr[(k + 1) % 3];
    nbr[2 * i + 1] = opposite_vert(triangles, adjacencies, tid, vid);
  }
  return nbr;
}
// returns the index of the triangles in the star of from which v is pointing
// to
int next_tid(const vector<vector<float>> &angles,
             const vector<vec3f> &positions, const vector<vector<int>> &v2t,
             const vector<vec3i> &triangles, const vector<vec3f> &normals,
             const int from, const vec3f &v) {
  auto teta =
      angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
  if (teta < 0)
    teta += 2 * M_PI;
  auto nbr = angles[from];
  int s = (int)nbr.size();
  if (teta == 0)
    return v2t[from][0];
  for (int i = 0; i < s; ++i) {
    if (nbr[i] < teta)
      continue;

    if (i % 2 == 0) {
      return v2t[from][(i - 2) / 2];
    } else {
      return v2t[from][(i - 1) / 2];
    }
  }
  return v2t[from].back();
}

int next_tid_extended_graph(const vector<float> &total_angles,
                            const vector<vec3f> &positions,
                            const vector<vector<int>> &v2t,
                            const vector<vec3i> &triangles,
                            const vector<vec3f> &normals, const int from,
                            const vec3f &v) {
  auto star = v2t[from];
  auto teta =
      angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
  auto accumulated = 0.f;
  auto scale_factor = 2 * pif / total_angles[from];
  auto vert_pos = positions[from];
  if (teta == 0)
    return star[0];
  if (teta < 0)
    teta += 2 * pif;
  for (int i = 0; i < star.size(); ++i) {
    auto tid = star[i];
    auto offset = find(triangles[tid], from);
    auto vid0 = triangles[tid][(offset + 1) % 3];
    auto vid1 = triangles[tid][(offset + 2) % 3];
    auto curr_angle =
        angle(positions[vid0] - vert_pos, positions[vid1] - vert_pos) *
        scale_factor;
    accumulated += curr_angle;

    if (accumulated > teta)
      return tid;
  }

  return star.back();
}
std::tuple<float, float, int>
next_tid_in_tangent_space(const shape_data &data,
                          const shape_geometry &geometry, const int from,
                          const float &theta) {
  auto nbr = geometry.angles[from];
  int s = nbr.size();
  if (theta == 0)
    return {0.f, nbr[1], geometry.v2t[from][0]};
  for (int i = 0; i < s; ++i) {
    if (nbr[i] < theta)
      continue;
    return {nbr[i - 1], nbr[i], geometry.v2t[from][i - 1]};
  }

  return {nbr.back(), 2 * M_PI, geometry.v2t[from].back()};
}
vec3f trace_segment_vert(const vector<vec3f> &verts, const vec3f n,
                         const vec3f bary, const vec3f baryV, const vec3f &dir,
                         const int offset) {
  auto right = verts[(offset + 1) % 3] - verts[offset];
  auto left = verts[(offset + 2) % 3] - verts[offset];
  auto sample_bary = zero3f;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, left), n) > 0) {
      auto factor = bary[offset] / baryV[offset];
      sample_bary[offset] = 0;
      sample_bary[(offset + 1) % 3] =
          bary[(offset + 1) % 3] - baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] =
          bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;
    } else
      sample_bary[(offset + 2) % 3] = 1;
  } else
    sample_bary[(offset + 1) % 3] = 1;

  return sample_bary;
}

vec3f trace_segment_edge(const vector<vec3f> &verts, const vec3f n,
                         const vec3f bary, const vec3f baryV, const vec3f &dir,
                         const int offset, const vec3f &sample_coords) {
  auto sample_bary = zero3f;
  auto right = verts[(offset + 1) % 3] - sample_coords;
  auto left = verts[offset] - sample_coords;
  auto front = verts[(offset + 2) % 3] - sample_coords;
  if (dot(cross(right, dir), n) > 0) {
    if (dot(cross(dir, front), n) > 0) {
      auto factor = bary[offset] / baryV[offset];
      sample_bary[offset] = 0;
      sample_bary[(offset + 1) % 3] =
          bary[(offset + 1) % 3] - baryV[(offset + 1) % 3] * factor;
      sample_bary[(offset + 2) % 3] =
          bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;

    } else {
      if (dot(cross(dir, left), n) > 0) {
        auto factor = bary[(offset + 1) % 3] / baryV[(offset + 1) % 3];
        sample_bary[(offset + 1) % 3] = 0;
        sample_bary[offset] = bary[offset] - baryV[offset] * factor;
        sample_bary[(offset + 2) % 3] =
            bary[(offset + 2) % 3] - baryV[(offset + 2) % 3] * factor;

      } else {
        if (dot(left, dir) > 0) {
          sample_bary[(offset)] = 1;
        } else {
          sample_bary[(offset + 1) % 3] = 1;
        }
      }
    }
  } else {
    if (dot(right, dir) > 0) {
      sample_bary[(offset + 1) % 3] = 1;
    } else {
      sample_bary[(offset)] = 1;
    }
  }

  return sample_bary;
}
vec3f trace_segment_tri(const vector<vec3f> &verts, const vec3f n,
                        const vec3f bary, const vec3f baryV, const vec3f &dir,
                        const vec3f &sample_coords) {
  auto sample_bary = zero3f;
  vec3f w0 = verts[0] - sample_coords;
  vec3f w1 = verts[1] - sample_coords;
  vec3f w2 = verts[2] - sample_coords;
  if (dot(cross(w0, dir), n) > 0 && dot(cross(dir, w1), n) > 0) {
    sample_bary = vec3f{bary[0] - bary[2] * baryV[0] / baryV[2],
                        bary[1] - bary[2] * baryV[1] / baryV[2], 0};

  } else if (dot(cross(w1, dir), n) > 0 && dot(cross(dir, w2), n) > 0) {
    sample_bary = vec3f{0, bary[1] - bary[0] * baryV[1] / baryV[0],
                        bary[2] - bary[0] * baryV[2] / baryV[0]};
  } else {
    sample_bary = vec3f{bary[0] - bary[1] * baryV[0] / baryV[1], 0,
                        bary[2] - bary[1] * baryV[2] / baryV[1]};
  }

  return sample_bary;
}
// Identify the intersection of the polyline inside triangle pid
// tracing: https://cims.nyu.edu/gcl/papers/campen2016bms.pdf
void trace_in_triangles_old(const vector<vec3f> &positions,
                            const vector<vec3i> &triangles, const vec3f &dir,
                            const vec3f &bary, const int pid, vec3f &sample_pos,
                            vec3f &sample_bary) {
  vec3f baryM = zero3f, baryV = zero3f;
  vec3f v0 = positions[triangles[pid].x];
  vec3f v1 = positions[triangles[pid].y];
  vec3f v2 = positions[triangles[pid].z];
  vector<vec3f> verts = {v0, v1, v2};
  vec3f n = triangle_normal(v0, v1, v2);
  vec3f sample_coords = bary.x * v0 + bary.y * v1 + bary.z * v2;
  vec3f M = sample_coords + dir;

  baryM = tri_bary_coords(v0, v1, v2, M);
  for (int i = 0; i < 3; ++i) {
    baryV[i] = baryM[i] - bary[i];
  }
  auto [is_vertex, k_vert] = bary_is_vert(bary);
  auto [is_on_edge, k_edge] = bary_is_edge(bary);
  if (is_vertex) {
    sample_bary = trace_segment_vert(verts, n, bary, baryV, dir, k_vert);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else if (is_on_edge) {
    sample_bary =
        trace_segment_edge(verts, n, bary, baryV, dir, k_edge, sample_coords);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  } else {
    sample_bary = trace_segment_tri(verts, n, bary, baryV, dir, sample_coords);
    sample_pos = sample_bary.x * verts[0] + sample_bary.y * verts[1] +
                 sample_bary.z * verts[2];
  }
}
void trace_in_triangles(const vector<vec3f> &positions,
                        const vector<vec3i> &triangles, const vec3f &dir,
                        const vec3f &bary, const int pid, vec3f &sample_pos,
                        vec3f &sample_bary) {
  sample_bary = zero3f;
  sample_pos = zero3f;
  auto [is_vertex, k_vert] = bary_is_vert(bary);
  auto [is_on_edge, k_edge] = bary_is_edge(bary);
  if (is_vertex) {
    auto v = positions[triangles[pid][k_vert]];
    auto v0 = positions[triangles[pid][(k_vert + 1) % 3]];
    auto v1 = positions[triangles[pid][(k_vert + 2) % 3]];
    auto theta = angle(v0 - v, dir);
    if (theta <= 1e-8) {
      sample_bary[(k_vert + 1) % 3] = 1;
      sample_pos = sample_bary[k_vert] * v +
                   sample_bary[(k_vert + 1) % 3] * v0 +
                   sample_bary[(k_vert + 2) % 3] * v1;
    } else if (std::abs(theta - angle(v0 - v, v1 - v)) < 1e-8) {
      sample_bary[(k_vert + 2) % 3] = 1;
      sample_pos = sample_bary[k_vert] * v +
                   sample_bary[(k_vert + 1) % 3] * v0 +
                   sample_bary[(k_vert + 2) % 3] * v1;
    } else {
      auto phi = angle(v - v0, v1 - v0);
      auto x = length(v0 - v) * std::sin(theta) /
               std::sin(pif - theta - angle(v - v0, v1 - v0));
      auto l = length(v1 - v0);
      sample_bary[(k_vert + 1) % 3] = 1 - x / l;
      sample_bary[(k_vert + 2) % 3] = x / l;
      sample_pos = sample_bary[k_vert] * v +
                   sample_bary[(k_vert + 1) % 3] * v0 +
                   sample_bary[(k_vert + 2) % 3] * v1;
    }
  } else if (is_on_edge) {
    auto v0 = positions[triangles[pid][k_edge]];
    auto v1 = positions[triangles[pid][(k_edge + 1) % 3]];
    auto v2 = positions[triangles[pid][(k_edge + 2) % 3]];
    auto n = tid_normal(triangles, positions, pid);
    auto p = bary[k_edge] * v0 + bary[(k_edge + 1) % 3] * v1 +
             bary[(k_edge + 2) % 3] * v2;
    auto d = v2 - p;
    if (angle(dir, d) < 1e-10) {
      sample_bary[(k_edge + 2) % 3] = 1;
      sample_pos = v2;
    } else if (dot(cross(d, dir), n) > 0) {
      auto theta = angle(dir, v1 - v0);
      auto k = (1 - bary[k_edge]) * length(v1 - v0);
      auto x =
          k * std::sin(theta) / std::sin(pif - angle(v1 - v0, v2 - v0) - theta);
      auto l = length(v2 - v0);
      sample_bary[k_edge] = 1 - x / l;
      sample_bary[(k_edge + 2) % 3] = x / l;
      sample_pos = sample_bary[k_edge] * v0 +
                   sample_bary[(k_edge + 1) % 3] * v1 +
                   sample_bary[(k_edge + 2) % 3] * v2;
    } else {
      auto theta = angle(dir, v0 - v1);
      auto k = bary[k_edge] * length(v1 - v0);
      auto x =
          k * std::sin(theta) / std::sin(pif - angle(v0 - v1, v2 - v1) - theta);
      auto l = length(v2 - v1);
      sample_bary[(k_edge + 1) % 3] = 1 - x / l;
      sample_bary[(k_edge + 2) % 3] = x / l;
      sample_pos = sample_bary[k_edge] * v0 +
                   sample_bary[(k_edge + 1) % 3] * v1 +
                   sample_bary[(k_edge + 2) % 3] * v2;
    }
  } else {
    std::cout << "tracing in triangles is WIP" << std::endl;
  }
}
void parallel_transp(const vector<vector<float>> &angles,
                     const vector<float> &total_angles,
                     const vector<vec3i> &triangles,
                     const vector<vec3f> &positions,
                     const vector<vec3i> &adjacencies,
                     const vector<vector<int>> &v2t, vec3f &v,
                     const vector<vec3f> &normals, const int &from,
                     const int &to, const int &mode) {
  switch (mode) {
  case V2V: {
    auto mag = length(v);
    float teta =
        angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
    auto nbr_from = extended_one_ring(triangles, adjacencies, v2t, from);
    auto nbr_to = extended_one_ring(triangles, adjacencies, v2t, to);
    float phi_ij = -1;
    float phi_ji = -1;

    for (int i = 0; i < nbr_from.size(); ++i) {
      if (nbr_from[i] == to) {
        phi_ij = angles[from][i];

        break;
      }
    }

    for (int j = 0; j < nbr_to.size(); ++j) {
      if (nbr_to[j] == from) {
        phi_ji = angles[to][j];
        break;
      }
    }
    assert(phi_ij != -1);
    assert(phi_ji != -1);

    vec3f e0 = polar_basis(triangles, positions, v2t, normals, to);
    float rotation = teta + phi_ji + M_PI - phi_ij;

    v = rot_vect(e0, normals[to], rotation);
    v *= mag;

  } break;

  case V2T: {
    float teta =
        angle_in_tangent_space(triangles, positions, v2t, normals, from, v);
    vec3i tri = triangles[to];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f normal = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);

    vec3f coords = positions[from] - centroid;
    float phi_ji = angle(e, coords);
    if (dot(cross(e, coords), normal) < 0)
      phi_ji = 2 * M_PI - phi_ji;
    int offset = find(tri, from);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];
    int vid2 = tri[(offset + 2) % 3];
    float factor = 2 * M_PI / total_angles[from];
    auto nbr_from = extended_one_ring(triangles, adjacencies, v2t, from);
    float phi_ij = -1;
    coords *= -1;
    if (nbr_from[0] == vid2) {
      vec3f edge = positions[vid2] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      curr_angle = 2 * M_PI - curr_angle;
      phi_ij = curr_angle;
    } else {
      for (int i = 0; i < nbr_from.size(); ++i) {
        if (nbr_from[i] == vid1) {
          phi_ij = angles[from][i];
          break;
        }
      }

      vec3f edge = positions[vid1] - positions[from];
      float curr_angle = angle(edge, coords);
      curr_angle *= factor;
      phi_ij += curr_angle;
    }

    float rot = teta + phi_ji + M_PI - phi_ij;

    e *= length(v);
    v = rot_vect(e, normal, rot);

  }

  break;

  case T2V: {
    vec3i tri = triangles[from];
    vec3f p0 = positions[tri.x];
    vec3f p1 = positions[tri.y];
    vec3f p2 = positions[tri.z];
    vec3f n = triangle_normal(p0, p1, p2);
    vec3f centroid = (p0 + p1 + p2) / 3.0;
    vec3f e = normalize(p0 - centroid);
    float teta = angle(e, v);

    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;
    int offset = find(tri, to);
    assert(offset != -1);
    int vid1 = tri[(offset + 1) % 3];

    vec3f vert = positions[tri[offset]];
    vec3f v1 = positions[vid1] - vert;

    vec3f coords = vert - centroid;
    float phi_ij = angle(e, coords);
    if (dot(cross(e, coords), n) < 0)
      phi_ij = 2 * M_PI - phi_ij;

    coords *= -1;
    float phi_ji = angle(v1, coords);
    float factor = 2 * M_PI / total_angles[to];
    phi_ji *= factor;
    auto nbr = extended_one_ring(triangles, adjacencies, v2t, to);
    for (int i = 0; i < nbr.size(); ++i) {
      if (nbr[i] == vid1) {
        float phi = angles[to][i];
        phi_ji += phi;
        break;
      }
    }

    float rot = teta + phi_ji + M_PI - phi_ij;
    vec3f e0 = polar_basis(triangles, positions, v2t, normals, to);
    e0 *= length(v);
    v = rot_vect(e0, normals[to], rot);

  }

  break;

  case T2T: {
    auto flat_from = init_flat_triangle(positions, triangles[from]);
    auto k = find(adjacencies[from], to);
    assert(k != -1);
    auto flat_to = unfold_face(triangles, positions, flat_from, from, to);
    auto bary = vec2f{0.333, 0.333};
    auto c0 =
        interpolate_triangle(flat_from[0], flat_from[1], flat_from[2], bary);
    auto c1 = interpolate_triangle(flat_to[0], flat_to[1], flat_to[2], bary);
    auto e0 = flat_from[0] - c0;
    auto e1 = flat_to[0] - c1;

    auto w = c1 - c0;
    auto phi_ij = angle(e0, w);
    if (cross(e0, w) < 0)
      phi_ij = 2 * M_PI - phi_ij;
    w *= -1;
    auto phi_ji = angle(e1, w);

    if (cross(e1, w) < 0)
      phi_ji = 2 * M_PI - phi_ji;

    auto n = tid_normal(triangles, positions, from);
    auto e = polar_basis(triangles, positions, from);
    float teta = angle(e, v);
    if (dot(cross(e, v), n) < 0)
      teta = 2 * M_PI - teta;

    auto e_to = polar_basis(triangles, positions, to);
    auto n_to = tid_normal(triangles, positions, to);
    float rot = teta + phi_ji + M_PI - phi_ij;
    e_to *= length(v);
    v = rot_vect(e_to, n_to, rot);

  }

  break;
  }
}
// VTP
std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<double> points(3 * V);
  vector<uint> faces(3 * F);
  vector<double> f(V);

  for (int i = 0; i < V; ++i) {
    points[3 * i] = positions[i].x;
    points[3 * i + 1] = positions[i].y;
    points[3 * i + 2] = positions[i].z;
  }
  for (int i = 0; i < F; ++i) {
    faces[3 * i] = triangles[i].x;
    faces[3 * i + 1] = triangles[i].y;
    faces[3 * i + 2] = triangles[i].z;
  }

  return std::make_pair(points, faces);
}
std::pair<vector<double>, vector<unsigned>>
exact_geodesic_wrapper(const vector<vec3i> &triangles,
                       const vector<vec3f> &positions,
                       const mesh_point &source) {
  int V = positions.size() + 1;
  int F = triangles.size();
  auto tid = source.face;
  vector<double> points(3 * V);
  vector<uint> faces(3 * (F + 2));
  vector<float> f(V);
  auto pos = eval_position(triangles, positions, source);

  for (int i = 0; i < V; ++i)
    // if (v2t[i].size() == 0) continue;
    if (i != V - 1) {
      points[3 * i] = positions[i].x;
      points[3 * i + 1] = positions[i].y;
      points[3 * i + 2] = positions[i].z;
    } else {
      points[3 * i] = pos.x;
      points[3 * i + 1] = pos.y;
      points[3 * i + 2] = pos.z;
    }

  for (int i = 0; i < F; ++i) {
    if (i != tid) {
      faces[3 * i] = triangles[i].x;
      faces[3 * i + 1] = triangles[i].y;
      faces[3 * i + 2] = triangles[i].z;
    } else {
      faces[3 * i] = triangles[i].x;
      faces[3 * i + 1] = triangles[i].y;
      faces[3 * i + 2] = V - 1;

      faces[3 * F] = triangles[i].y;
      faces[3 * F + 1] = triangles[i].z;
      faces[3 * F + 2] = V - 1;

      faces[3 * (F + 1)] = triangles[i].z;
      faces[3 * (F + 1) + 1] = triangles[i].x;
      faces[3 * (F + 1) + 2] = V - 1;
    }
  }
  return {points, faces};
}
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const int &source) {
  int V = (int)positions.size();
  int F = (int)triangles.size();
  vector<float> f(V);
  auto [points, faces] = exact_geodesic_wrapper(triangles, positions);
  geodesic_VTP::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
  algorithm.propagate(source);
  vector<geodesic_VTP::Vertex> verts = mesh.vertices();
  for (int j = 0; j < V; ++j) {
    geodesic_VTP::Vertex v = verts[j];
    float value = (float)v.geodesic_distance();
    f[j] = value;
  }

  return f;
}

// utility (gradient matrix)
Eigen::VectorXd wrapper(const vector<float> &f) {
  Eigen::VectorXd F(f.size());
  for (int i = 0; i < f.size(); ++i) {
    F(i) = f[i];
  }
  return F;
}
vector<float> exact_geodesic_distance(const vector<vec3i> &triangles,
                                      const vector<vec3f> &positions,
                                      const mesh_point &source) {
  auto tid = source.face;
  auto [is_vert, offset] = point_is_vert(source);
  if (is_vert)
    return exact_geodesic_distance(triangles, positions,
                                   triangles[tid][offset]);
  else {
    int V = positions.size() + 1;
    int F = triangles.size();

    vector<float> f(V);
    auto [points, faces] = exact_geodesic_wrapper(triangles, positions, source);

    geodesic_VTP::Mesh mesh;
    mesh.initialize_mesh_data(points, faces);
    geodesic_VTP::GeodesicAlgorithmExact algorithm(&mesh);
    algorithm.propagate(V - 1);
    vector<geodesic_VTP::Vertex> verts = mesh.vertices();
    for (int j = 0; j < V; ++j) {
      geodesic_VTP::Vertex v = verts[j];
      float value = (float)v.geodesic_distance();
      f[j] = value;
    }
    f.pop_back();
    return f;
  }
}
Eigen::MatrixXd rhs(int s) {
  Eigen::MatrixXd E(s, s + 1);
  Eigen::MatrixXd X = Eigen::MatrixXd::Constant(s, 1, -1);
  E.topLeftCorner(s, 1) = X;
  Eigen::MatrixXd I(s, s);
  I.setIdentity();
  E.topRightCorner(s, s) = I;
  return E;
}
inline pair<Eigen::MatrixXd, Eigen::MatrixXi>
libigl_wrapper(const vector<vec3f> &positions, const vector<vec3i> &triangles) {
  Eigen::MatrixXd V(positions.size(), 3);
  Eigen::MatrixXi F(triangles.size(), 3);

  for (int i = 0; i < positions.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      V(i, j) = positions[i][j];
    }
  }
  for (int i = 0; i < triangles.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      F(i, j) = triangles[i][j];
    }
  }

  return {V, F};
}
double angle_at_vertex(const shape_data &data, const int tid, const int vid) {
  auto k = find(data.triangles[tid], vid);
  if (k == -1) {
    std::cout << "Error: the vertex must belongs to the triangles. Check the "
                 "construction of the v2t adjacency"
              << std::endl;
    return 0.f;
  }
  return angle(data.positions[data.triangles[tid][(k + 1) % 3]] -
                   data.positions[data.triangles[tid][k]],
               data.positions[data.triangles[tid][(k + 2) % 3]] -
                   data.positions[data.triangles[tid][k]]);
}
inline double cot(const float &theta) {
  return yocto::cos(theta) / yocto::sin(theta);
}
Eigen::VectorXd compute_hessian(const shape_op &operators,
                                const vector<float> &field) {
  auto n = field.size();
  auto f = wrapper(field);
  Eigen::VectorXd Hessian = operators.Hess * f;
  return Hessian;
}
Eigen::Matrix2d compute_hessian(const Eigen::VectorXd &Hessian, const int vid) {
  Eigen::Matrix2d H;
  auto n = Hessian.rows() / 4;
  H(0, 0) = Hessian(vid);
  H(0, 1) = Hessian(n + vid);
  H(1, 0) = Hessian(2 * n + vid);
  H(1, 1) = Hessian(3 * n + vid);

  return H;
}
vector<float> cotangent_laplacian(const shape_data &data,
                                  const shape_geometry &topology,
                                  const vector<float> &field) {
  auto L = cot_lap_stiffness_matrix(data, topology);
  auto M = mass_matrix(data, topology, true);
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  Eigen::SparseMatrix<double> I(data.positions.size(), data.positions.size());
  I.setIdentity();
  Eigen::SparseMatrix<double> M_inv = solver.solve(I);
  Eigen::SparseMatrix<double> A = M_inv * L;
  auto F = wrapper(field);
  Eigen::VectorXd Lap = A * F;
  vector<float> laplacian(field.size());
  for (auto i = 0; i < Lap.size(); ++i) {
    laplacian[i] = (float)Lap(i);
  }
  return laplacian;
}
vector<float> cotangent_laplacian(const Eigen::SparseMatrix<double> &L,
                                  const Eigen::SparseMatrix<double> &M,
                                  const vector<float> &field) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(M);
  Eigen::SparseMatrix<double> I(M.rows(), M.rows());
  I.setIdentity();
  Eigen::SparseMatrix<double> M_inv = solver.solve(I);
  Eigen::SparseMatrix<double> A = M_inv * L;
  auto F = wrapper(field);
  Eigen::VectorXd Lap = A * F;
  vector<float> laplacian(field.size());
  for (auto i = 0; i < Lap.size(); ++i) {
    laplacian[i] = (float)Lap(i);
  }
  return laplacian;
}
vec3f compute_PCE(const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, const vector<float> &f,
                  const mesh_point &p) {
  auto tid = p.face;
  vec3f n = cross(positions[triangles[tid].y] - positions[triangles[tid].x],
                  positions[triangles[tid].z] - positions[triangles[tid].x]);
  double area = length(n);
  n = normalize(n);
  auto vid0 = triangles[tid][0];
  auto vid1 = triangles[tid][1];
  auto vid2 = triangles[tid][2];
  auto v0 = positions[vid0];
  auto v1 = positions[vid1];
  auto v2 = positions[vid2];

  return 1 / area *
         ((f[vid1] - f[vid0]) * cross(n, v0 - v2) +
          (f[vid2] - f[vid0]) * cross(n, v1 - v0));
}
Eigen::SparseMatrix<double> PCE_matrix(const vector<vec3i> &triangles,
                                       const vector<vec3f> &positions) {
  time_function();
  Eigen::SparseMatrix<double> G(triangles.size() * 3, positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T> entries;

  for (int i = 0; i < triangles.size(); ++i) {
    vec3f n = cross(positions[triangles[i].y] - positions[triangles[i].x],
                    positions[triangles[i].z] - positions[triangles[i].x]);
    double area = length(n); // division by 2 is missing because
                             // there is a *2 in the computation of
                             // the gradient
    n = normalize(n);

    for (int off = 0; off < 3; ++off) {
      int prev = triangles[i][off];
      int curr = triangles[i][(off + 1) % 3];
      int next = triangles[i][(off + 2) % 3];
      vec3f u = positions[next] - positions[curr];
      vec3f v = positions[curr] - positions[prev];
      vec3f u_90 = normalize(cross(u, n));
      vec3f v_90 = normalize(cross(v, n));

      vec3f contribute = u_90 * length(u) + v_90 * length(v);
      contribute /= area;

      int row = 3 * i;
      entries.push_back(T(row, curr, contribute.x));
      ++row;
      entries.push_back(T(row, curr, contribute.y));
      ++row;
      entries.push_back(T(row, curr, contribute.z));
    }
  }
  G.setFromTriplets(entries.begin(), entries.end());
  return G;
}
std::pair<float, float> intersect(const vec2f &direction, const vec2f &left,
                                  const vec2f &right) {
  auto v1 = -left;
  auto v2 = right - left;
  auto v3 = vec2f{-direction.y, direction.x};
  auto t0 = cross(v2, v1) / dot(v2, v3);
  auto t1 = -dot(left, v3) / dot(v2, v3);
  return std::make_pair(t0, t1);
};
void trim_path(const vector<vec3i> &triangles, const vector<vec3f> &positions,
               const double &factor, const int tid, vector<mesh_point> &path) {
  auto next = eval_position(triangles, positions, path.back());
  auto prev = eval_position(triangles, positions, path.rbegin()[1]);
  auto w = normalize(prev - next);
  w *= factor;
  w += next;
  auto bary = tri_bary_coords(triangles, positions, tid, w);
  path.pop_back();
  path.push_back({tid, vec2f{bary[1], bary[2]}});
}
std::tuple<vector<float>, vector<uint>>
find_angles_in_tangent_space(const shape_data &data,
                             const shape_geometry &geometry, const int vid,
                             const vector<int> &neighbors) {
  auto result = vector<float>(neighbors.size());
  auto nbr = one_ring(data.triangles, geometry.adjacencies, geometry.v2t, vid);
  for (auto i = 0; i < nbr.size(); ++i) {
    auto it = std::find(neighbors.begin(), neighbors.end(), nbr[i]);
    if (it == neighbors.end())
      continue;
    auto entry = std::distance(neighbors.begin(), it);
    result[entry] = geometry.angles[vid][i];
  }

  return {result, nbr};
}
std::tuple<vec2f, vec2f, vec2f, int, int, int>
dir_left_right(const shape_data &data, const shape_geometry &geometry,
               const int vid, const float &phi) {
  auto theta = (float)std::fmod(phi + M_PI, 2 * M_PI);

  auto [phi0, phi1, new_tid] =
      next_tid_in_tangent_space(data, geometry, vid, theta);
  auto k = find(data.triangles[new_tid], vid);
  auto vid0 = data.triangles[new_tid][(k + 1) % 3];
  auto vid1 = data.triangles[new_tid][(k + 2) % 3];
  auto r0 = length(data.positions[vid0] - data.positions[vid]);
  auto r1 = length(data.positions[vid1] - data.positions[vid]);
  auto dir = vec2f{std::cos(theta), std::sin(theta)};
  auto right = vec2f{r0 * std::cos(phi0), r0 * std::sin(phi0)};
  auto left = vec2f{r1 * std::cos(phi1), r1 * std::sin(phi1)};

  return {dir, left, right, vid0, vid1, new_tid};
}
std::tuple<int, int, unfold_triangle, vec2f>
init_next_flat_tid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const int vid, const int tid,
                   const int vid0, const int vid1, const double &t1,
                   vector<mesh_point> &path) {
  auto k = find(triangles[tid], vid);
  auto bary = vec3f{0, 0, 0};
  bary[(k + 1) % 3] = t1;
  bary[(k + 2) % 3] = 1 - t1;
  auto q = mesh_point{tid, vec2f{bary.y, bary.z}};
  path.push_back(q);
  auto flat_tid = init_flat_triangle(positions, triangles[tid], k);
  auto v = t1 * positions[vid0] + (1 - t1) * positions[vid1] - positions[vid];
  auto x = normalize(positions[vid0] - positions[vid]);
  auto l = dot(v, x);
  auto flat_dir = normalize(vec2f{l, length((v - l * x))});
  return {tid, tid, flat_tid, flat_dir};
}
std::tuple<int, int, unfold_triangle, vec2f>
init_next_flat_tid(const vector<vec3i> &triangles,
                   const vector<vec3f> &positions, const vector<vec3i> &t2t,
                   const int vid, const int tid, const int vid0, const int vid1,
                   const double &t1, vector<mesh_point> &path) {
  auto k = find(triangles[tid], vid);
  auto bary = vec3f{0, 0, 0};
  bary[(k + 1) % 3] = t1;
  bary[(k + 2) % 3] = 1 - t1;
  auto q = mesh_point{tid, vec2f{bary.y, bary.z}};
  path.push_back(q);
  auto flat_tid = init_flat_triangle(positions, triangles[tid], k);
  flat_tid =
      unfold_face(triangles, positions, flat_tid, tid, t2t[tid][(k + 1) % 3]);

  auto v = t1 * positions[vid0] + (1 - t1) * positions[vid1] - positions[vid];
  auto x = normalize(positions[vid0] - positions[vid]);
  auto l = dot(v, x);
  auto flat_dir = normalize(vec2f{l, length(v - l * x)});
  auto prev_tid = tid;
  auto new_tid = t2t[tid][(k + 1) % 3];
  return {new_tid, prev_tid, flat_tid, flat_dir};
}
std::tuple<int, int, unfold_triangle, vec2f>
find_next_triangle(const shape_data &data, const shape_geometry &geometry,
                   const int vid, const float &phi, const double &prev_len,
                   const double &max_len, vector<mesh_point> &path,
                   double &path_len) {
  auto [dir, left, right, vid0, vid1, new_tid] =
      dir_left_right(data, geometry, vid, phi);
  auto [t0, t1] = intersect(dir, left, right);
  auto curr_vid = vid;
  if (path_len + t0 + prev_len >= max_len) {
    path_len += t0;
    return init_next_flat_tid(data.triangles, data.positions, curr_vid, new_tid,
                              vid0, vid1, t1, path);
  }

  while (std::abs(t1) < 1e-4 || std::abs(1 - t1) < 1e-4) {
    path.push_back(point_from_vert(data.triangles, geometry.v2t, curr_vid));
    path_len += t0;
    if (std::abs(t1) < 1e-4) {
      curr_vid = vid1;
      if (geometry.angles[curr_vid].size() == 0)
        return init_next_flat_tid(data.triangles, data.positions, curr_vid,
                                  new_tid, vid0, vid1, t1, path);
      auto [angle, nbr] =
          find_angles_in_tangent_space(data, geometry, curr_vid, {vid1});
      auto theta = (float)std::fmod(angle[0] + M_PI, 2 * M_PI);
      std::tie(dir, left, right, vid0, vid1, new_tid) =
          dir_left_right(data, geometry, curr_vid, theta);
      std::tie(t0, t1) = intersect(dir, left, right);
      if (path_len + t0 + prev_len >= max_len) {
        path_len += t0;
        return init_next_flat_tid(data.triangles, data.positions, curr_vid,
                                  new_tid, vid0, vid1, t1, path);
      }

    } else {
      curr_vid = vid0;
      if (geometry.angles[curr_vid].size() == 0)
        return init_next_flat_tid(data.triangles, data.positions, curr_vid,
                                  new_tid, vid0, vid1, t1, path);
      auto [angle, nbr] =
          find_angles_in_tangent_space(data, geometry, curr_vid, {vid0});
      auto theta = (float)std::fmod(angle[0] + M_PI, 2 * M_PI);
      std::tie(dir, left, right, vid0, vid1, new_tid) =
          dir_left_right(data, geometry, curr_vid, theta);
      std::tie(t0, t1) = intersect(dir, left, right);
      if (path_len + t0 + prev_len >= max_len) {
        path_len += t0;
        return init_next_flat_tid(data.triangles, data.positions, curr_vid,
                                  new_tid, vid0, vid1, t1, path);
      }
    }
  }

  return init_next_flat_tid(data.triangles, data.positions,
                            geometry.adjacencies, curr_vid, new_tid, vid0, vid1,
                            t1, path);
}
std::tuple<int, int, unfold_triangle, vec2f>
handle_vert(const shape_data &data, const shape_geometry &geometry,
            const int vid, const int tid, const mesh_point &from,
            const double &prev_len, const double &max_len,
            vector<mesh_point> &path, double &path_len) {
  auto k = find(data.triangles[tid], vid);
  auto vid0 = data.triangles[tid][(k + 1) % 3];
  auto vid1 = data.triangles[tid][(k + 2) % 3];
  auto [angles, nbr] =
      find_angles_in_tangent_space(data, geometry, vid, {vid0, vid1});
  auto r0 = length(data.positions[vid0] - data.positions[vid]);
  auto r1 = length(data.positions[vid1] - data.positions[vid]);

  auto v0 = vec2f{r0 * std::cos(angles[0]), r0 * std::sin(angles[0])};
  auto v1 = vec2f{r1 * std::cos(angles[1]), r1 * std::sin(angles[1])};
  auto bary = get_bary(from.uv);
  auto p = bary[(k + 1) % 3] * v0 + bary[(k + 2) % 3] * v1;
  auto phi = angles[0] + angle(p, v0);

  return find_next_triangle(data, geometry, vid, phi, prev_len, max_len, path,
                            path_len);
}

vector<mesh_point> straightest_path_from_vert_in_triangle(
    const shape_data &data, const shape_geometry &geometry, const int tid,
    const int vid, const float &len, const vec3f &dir) {
  auto path = vector<mesh_point>{};
  auto k = find(data.triangles[tid], vid);
  auto flat_tid = init_flat_triangle(data.positions, data.triangles[tid], k);
  auto v = normalize(data.positions[data.triangles[tid][(k + 1) % 3]] -
                     data.positions[vid]);
  auto l = dot(dir, v);
  auto flat_dir = normalize(vec2f{l, length(dir - l * v)});
  auto bary = vec3f{0, 0, 0};
  bary[k] = 1;
  path.push_back(mesh_point{tid, vec2f{bary.y, bary.z}});

  auto [t0, t1] =
      intersect(flat_dir, flat_tid[(k + 2) % 3], flat_tid[(k + 1) % 3]);
  bary = vec3f{0, 0, 0};
  bary[(k + 1) % 3] = t1;
  bary[(k + 2) % 3] = 1 - t1;
  auto q = mesh_point{tid, vec2f{bary.y, bary.z}};
  path.push_back(q);

  auto path_len = t0;
  if (t0 > len) {
    trim_path(data.triangles, data.positions, t0 - len, tid, path);
    return path;
  }
  auto prev_tid = tid;
  auto curr_tid = -1;
  auto curr_vid = -1;
  double accumulated = 0;
  if (std::abs(t1) < 1e-4 || std::abs(1 - t1) < 1e-4) {
    accumulated = t0;
    curr_vid = vert_from_point(data.triangles, q);
    if (geometry.angles[curr_vid].size() == 0)
      return path;
    std::tie(curr_tid, prev_tid, flat_tid, flat_dir) =
        handle_vert(data, geometry, curr_vid, tid, path[0], path_len, len, path,
                    accumulated);
    path_len = accumulated;
  } else {
    curr_tid = geometry.adjacencies[tid][(k + 1) % 3];
    flat_tid =
        unfold_face(data.triangles, data.positions, flat_tid, tid, curr_tid);
  }

  while (path_len < len) {
    for (auto h = 0; h < 3; ++h) {
      auto neighbor = geometry.adjacencies[curr_tid][h];
      if (neighbor == prev_tid)
        continue;
      if (neighbor == -1)
        return path;
      auto left = flat_tid[(h + 1) % 3];
      auto right = flat_tid[h];
      std::tie(t0, t1) = intersect(flat_dir, left, right);
      if (t0 <= 0)
        continue;
      if (std::abs(t1) < 1e-6 || std::abs(1 - t1) < 1e-6) {
        accumulated = t0;
        curr_vid = (std::abs(t1) < 1e-6) ? data.triangles[curr_tid][(h + 1) % 3]
                                         : data.triangles[curr_tid][h];
        std::tie(curr_tid, prev_tid, flat_tid, flat_dir) =
            handle_vert(data, geometry, curr_vid, curr_tid, path.back(),
                        path_len, len, path, accumulated);
      } else if (t1 > 0 && t1 < 1) {
        path_len = t0 + accumulated;
        bary = vec3f{0, 0, 0};
        bary[h] = t1;
        bary[(h + 1) % 3] = 1 - t1;
        path.push_back(mesh_point{curr_tid, vec2f{bary.y, bary.z}});
        prev_tid = curr_tid;
        flat_tid = unfold_face(data.triangles, data.positions, flat_tid,
                               curr_tid, neighbor);
        curr_tid = geometry.adjacencies[curr_tid][h];
        break;
      }
    }
  }
  auto factor = path_len - len;
  trim_path(data.triangles, data.positions, factor, prev_tid, path);

  return path;
}
pair<int, vec3f> handle_vert(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions,
                             const vector<vec3f> &normals,
                             const vector<vec3i> &adjacencies,
                             const vector<vector<int>> &v2p_adjacencies,
                             const vector<vector<float>> &angles,
                             const vector<float> &total_angles, const int vid,
                             const int tid, const vec3f dir) {
  auto v = dir;
  parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                  v2p_adjacencies, v, normals, tid, vid, T2V);
  auto next = next_tid_extended_graph(total_angles, positions, v2p_adjacencies,
                                      triangles, normals, vid, v);

  /* auto next = next_tid_extended_graph(
      solver, angles, positions, v2p_adjacencies, triangles, normals, vid,
     v);
   */

  parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                  v2p_adjacencies, v, normals, vid, next, V2T);
  return std::make_pair(next, v);
}
vector<mesh_point> straightest_geodesic(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3f> &normals, const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const bool check_from) {
  auto vid = -1, tid = -1, next_tri = -1;
  auto dir = v;
  float len = 0.0;
  auto next_bary = zero3f, next = zero3f;
  auto prev = eval_position(triangles, positions, from);
  auto samples = vector<mesh_point>{from};
  auto bary = get_bary(from.uv);

  if (check_from) {
    auto [is_vert, kv] = bary_is_vert(bary);
    auto [is_on_edge, ke] = bary_is_edge(bary);
    if (is_vert) {
      vid = triangles[from.face][kv];
      // parallel_transp(angles, total_angles, triangles, positions,
      //                 adjac v2p_adjacencies,encies, dir, normals,
      //                 from.face, vid, T2V);

      tid = next_tid(angles, positions, v2p_adjacencies, triangles, normals,
                     vid, v);
      // next_tid_extended_graph(total_angles, positions, v2p_adjacencies,
      //                               triangles, normals, vid, v);

      /* tid = next_tid_extended_graph(
          solver, angles, positions, v2p_adjacencies, triangles, normals,
         vid, v);
       */
      kv = find(triangles[tid], vid);
      bary = zero3f;
      bary[kv] = 1;
      parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                      v2p_adjacencies, dir, normals, vid, tid, V2T);
    } else if (is_on_edge) {
      auto p0 = triangles[from.face][ke];
      auto p1 = triangles[from.face][(ke + 1) % 3];
      auto p2 = triangles[from.face][(ke + 2) % 3];
      auto n = triangle_normal(positions[p0], positions[p1], positions[p2]);
      auto edge = normalize(positions[p1] - positions[p0]);
      if (dot(cross(edge, v), n) > 0)
        tid = from.face;
      else {
        tid = adjacencies[from.face][ke];
        bary = tri_bary_coords(triangles, positions, tid, prev);
        parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                        v2p_adjacencies, dir, normals, from.face, tid, T2T);
      }
    } else
      tid = from.face;
  } else
    tid = from.face;

  while (len < l) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    samples.push_back({tid, vec2f{next_bary.y, next_bary.z}});
    len += length(next - prev);
    if (len < l) {
      prev = next;
      auto [V, k_v] = bary_is_vert(next_bary);
      auto [E, k_e] = bary_is_edge(next_bary);
      if (V) {
        vid = triangles[tid][k_v];
        if (angles[vid].size() == 0)
          return samples;
        auto out =
            handle_vert(triangles, positions, normals, adjacencies,
                        v2p_adjacencies, angles, total_angles, vid, tid, dir);
        tid = out.first;
        dir = out.second;
        k_v = find(triangles[tid], vid);
        bary = zero3f;
        bary[k_v] = 1;
      } else if (E) {
        next_tri = adjacencies[tid][k_e];
        if (next_tri == -1)
          return samples;
        auto p0 = triangles[tid][k_e];
        auto offset0 = find(triangles[next_tri], p0);
        auto offset1 = (offset0 + 2) % 3;
        bary = zero3f;
        bary[offset0] = next_bary[k_e];
        bary[offset1] = next_bary[(k_e + 1) % 3];

        parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                        v2p_adjacencies, dir, normals, tid, next_tri, T2T);
        tid = next_tri;
      } else
        assert(false);
    }
  }

  auto factor = (len - l);
  auto w = normalize(prev - next);
  w *= factor;
  w += next;
  bary = tri_bary_coords(triangles, positions, tid, w);
  samples.pop_back();
  samples.push_back({tid, vec2f{bary.y, bary.z}});

  return samples;
}
pair<vector<mesh_point>, float> straightest_geodesic(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const vector<vec3f> &normals, const vector<vec3i> &adjacencies,
    const vector<vector<int>> &v2p_adjacencies,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const mesh_point &from, const vec3f &v, const float &l,
    const int max_crossed_tri) {
  auto vid = -1, tid = -1, next_tri = -1;
  auto dir = v;
  float len = 0.0;
  auto next_bary = zero3f, next = zero3f;
  auto prev = eval_position(triangles, positions, from);
  auto samples = vector<mesh_point>{from};
  auto bary = get_bary(from.uv);
  auto [is_vert, kv] = bary_is_vert(bary);
  auto [is_on_edge, ke] = bary_is_edge(bary);
  auto crossed_tri = 0;
  if (is_vert) {
    vid = triangles[from.face][kv];
    tid = next_tid(angles, positions, v2p_adjacencies, triangles, normals, vid,
                   v);
    kv = find(triangles[tid], vid);
    bary = zero3f;
    bary[kv] = 1;
    parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                    v2p_adjacencies, dir, normals, vid, tid, V2T);
  } else if (is_on_edge) {
    auto p0 = triangles[from.face][ke];
    auto p1 = triangles[from.face][(ke + 1) % 3];
    auto p2 = triangles[from.face][(ke + 2) % 3];
    auto n = triangle_normal(positions[p0], positions[p1], positions[p2]);
    auto edge = normalize(positions[p1] - positions[p0]);
    if (dot(cross(edge, v), n) > 0)
      tid = from.face;
    else {
      tid = adjacencies[from.face][ke];
      bary = tri_bary_coords(triangles, positions, tid, prev);
      parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                      v2p_adjacencies, dir, normals, from.face, tid, T2T);
    }

  } else
    tid = from.face;

  while (len < l && crossed_tri < max_crossed_tri) {
    trace_in_triangles(positions, triangles, dir, bary, tid, next, next_bary);
    samples.push_back({tid, vec2f{next_bary.y, next_bary.z}});
    len += length(next - prev);
    ++crossed_tri;
    if (len < l) {
      prev = next;
      auto [V, k_v] = bary_is_vert(next_bary);
      auto [E, k_e] = bary_is_edge(next_bary);
      if (V) {
        vid = triangles[tid][k_v];
        auto out =
            handle_vert(triangles, positions, normals, adjacencies,
                        v2p_adjacencies, angles, total_angles, vid, tid, dir);
        tid = out.first;
        dir = out.second;
        k_v = find(triangles[tid], vid);
        bary = zero3f;
        bary[k_v] = 1;
      } else if (E) {
        next_tri = adjacencies[tid][k_e];
        auto p0 = triangles[tid][k_e];
        auto offset0 = find(triangles[next_tri], p0);
        auto offset1 = (offset0 + 2) % 3;
        bary = zero3f;
        bary[offset0] = next_bary[k_e];
        bary[offset1] = next_bary[(k_e + 1) % 3];

        parallel_transp(angles, total_angles, triangles, positions, adjacencies,
                        v2p_adjacencies, dir, normals, tid, next_tri, T2T);
        tid = next_tri;
      } else
        assert(false);
    }
  }
  if (len > l) {
    auto factor = (len - l);
    auto w = normalize(prev - next);
    w *= factor;
    w += next;
    bary = tri_bary_coords(triangles, positions, tid, w);
    samples.pop_back();
    samples.push_back({tid, vec2f{bary.y, bary.z}});
    len = l;
  } else {
    len = l;
    samples.push_back({tid, {bary.y, bary.z}});
  }

  return {samples, len};
}

float divergence(const shape_data &data, const shape_geometry &topology,
                 const Eigen::VectorXd &pce_grad, const int vid) {
  auto div = 0.f;
  auto star = topology.v2t[vid];
  for (auto i = 0; i < star.size(); ++i) {
    auto grad =
        vec3f{(float)pce_grad(3 * star[i]), (float)pce_grad(3 * star[i] + 1),
              (float)pce_grad(3 * star[i] + 1)};
    auto k = find(data.triangles[star[i]], vid);
    auto e1 = data.positions[data.triangles[star[i]][(k + 1) % 3]] -
              data.positions[data.triangles[star[i]][k]];
    auto e2 = data.positions[data.triangles[star[i]][(k + 2) % 3]] -
              data.positions[data.triangles[star[i]][k]];
    auto theta_1 =
        angle_at_vertex(data, star[i], data.triangles[star[i]][(k + 1) % 3]);
    auto theta_2 =
        angle_at_vertex(data, star[i], data.triangles[star[i]][(k + 2) % 3]);
    auto contribute =
        cot(theta_1) * dot(e1, grad) + cot(theta_2) * dot(e2, grad);
    div += contribute;
  }
  return div / 2.f;
}
Eigen::VectorXd divergence(const shape_data &data,
                           const shape_geometry &topology,
                           const Eigen::VectorXd &pce_grad) {
  Eigen::VectorXd div(data.positions.size());
  for (auto i = 0; i < data.positions.size(); ++i) {
    div(i) = divergence(data, topology, pce_grad, i);
  }
  return div;
}

void fill_riemannian_gradient_entries(
    vector<Eigen::Triplet<double>> &entries,
    const vector<vector<pair<int, float>>> &ring, const Eigen::VectorXd &c,
    const Eigen::VectorXd &a0, const Eigen::VectorXd &a1, const int n) {
  int vid = ring[0][0].first;
  int s = (int)ring.size();
  double c0_squared = pow(c[0], 2);
  double c1_squared = pow(c[1], 2);
  Eigen::Matrix2d g_inv;
  double det = 1 + c0_squared + c1_squared;
  g_inv << 1 + c1_squared, -c[0] * c[1], -c[0] * c[1], 1 + c0_squared;
  g_inv /= det;
  typedef Eigen::Triplet<double> T;
  for (int i = 0; i < s; ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;
      auto w = ring[i][j].second;
      entries.push_back(
          T(vid, entry, w * (g_inv(0, 0) * a0(i) + g_inv(0, 1) * a1(i))));
      entries.push_back(
          T(n + vid, entry, w * (g_inv(1, 0) * a0(i) + g_inv(1, 1) * a1(i))));
    }
  }
}
void fill_riemannian_gradient_entries(
    vector<Eigen::Triplet<double>> &entries,
    const vector<vector<pair<int, float>>> &ring, const Eigen::Matrix2d &g,
    const double &det, const Eigen::VectorXd &c1, const Eigen::VectorXd &c2,
    const vec3f &xu, const vec3f &xv, const int n) {
  int vid = ring[0][0].first;
  int s = (int)ring.size();
  typedef Eigen::Triplet<double> T;
  auto g_nabla_u = 1 / det * ((float)g(1, 1) * xu - (float)g(0, 1) * xv);
  auto g_nabla_v = 1 / det * ((float)g(0, 0) * xv - (float)g(0, 1) * xu);
  for (int i = 0; i < s; ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;
      auto w = ring[i][j].second;

      entries.push_back(
          T(vid, entry, w * (g_nabla_u[0] * c1(i) + g_nabla_v[0] * c2(i))));
      entries.push_back(
          T(n + vid, entry, w * (g_nabla_u[1] * c1(i) + g_nabla_v[1] * c2(i))));
      entries.push_back(T(2 * n + vid, entry,
                          w * (g_nabla_u[2] * c1(i) + g_nabla_v[2] * c2(i))));
    }
  }
}
void fill_hessian_entries(vector<Eigen::Triplet<double>> &entries,
                          const Eigen::Matrix2d &inv_g,
                          const vector<vector<pair<int, float>>> &ring,
                          const Eigen::MatrixXd &C, const vec3f &xu,
                          const vec3f &xv, const vec3f &xuu, const vec3f &xuv,
                          const vec3f &xvv, const int n) {
  typedef Eigen::Triplet<double> T;
  int vid = ring[0][0].first;
  int s = (int)ring.size();
  Eigen::Matrix2d gmj0;
  Eigen::Matrix2d gmj1;
  Eigen::Matrix2d gjm0;
  Eigen::Matrix2d gjm1;
  Eigen::Matrix2d g0jm;
  Eigen::Matrix2d g1jm;
  Eigen::RowVector2d g1m;
  Eigen::RowVector2d g0m;
  Eigen::MatrixXd C0m(2, s);
  Eigen::MatrixXd C1m(2, s);
  Eigen::MatrixXd C1(2, s);
  g0m << inv_g(0, 0), inv_g(0, 1);
  g1m << inv_g(1, 0), inv_g(1, 1);
  gmj0 << dot(xu, xuu), dot(xu, xuv), dot(xv, xuu), dot(xv, xuv);
  gmj1 << dot(xu, xuv), dot(xu, xvv), dot(xv, xuv), dot(xv, xvv);
  gjm0 = gmj0.transpose();
  gjm1 = gmj1.transpose();
  g0jm << dot(xu, xuu), dot(xu, xuv), dot(xu, xuv), dot(xu, xvv);
  g1jm << dot(xv, xuu), dot(xv, xuv), dot(xv, xuv), dot(xv, xvv);
  C1.row(0) = C.row(1);
  C1.row(1) = C.row(2);
  C0m.row(0) = C.row(3);
  C0m.row(1) = C.row(4);
  C1m.row(0) = C.row(4);
  C1m.row(1) = C.row(5);

  Eigen::VectorXd CoeffA;
  Eigen::VectorXd CoeffB;
  Eigen::VectorXd CoeffC;
  Eigen::VectorXd CoeffD;

  CoeffA = -inv_g.row(0) * gjm0 * inv_g.transpose() * C1 +
           inv_g.col(0).transpose() * C0m;
  CoeffB = -inv_g.row(0) * gjm1 * inv_g.transpose() * C1 +
           inv_g.col(0).transpose() * C1m;
  CoeffC = -inv_g.row(1) * gjm0 * inv_g.transpose() * C1 +
           inv_g.col(1).transpose() * C0m;
  CoeffD = -inv_g.row(1) * gjm1 * inv_g.transpose() * C1 +
           inv_g.col(1).transpose() * C1m;

  for (int i = 0; i < s; ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;
      auto w = ring[i][j].second;

      entries.push_back(T(vid, entry, w * CoeffA(i)));
      entries.push_back(T(n + vid, entry, w * CoeffB(i)));
      entries.push_back(T(2 * n + vid, entry, w * CoeffC(i)));
      entries.push_back(T(3 * n + vid, entry, w * CoeffD(i)));
    }
  }

  // typedef Eigen::Triplet<double> T;
  // int vid = ring[0][0].first;
  // int s = (int)ring.size();
  // Eigen::Matrix2d gj0m;
  // Eigen::Matrix2d gj1m;
  // Eigen::Vector2d g0m;
  // Eigen::Vector2d g1m;
  // Eigen::MatrixXd C1m(2, s);
  // Eigen::MatrixXd C2m(2, s);
  // Eigen::MatrixXd C1(2, s);
  // g0m << inv_g(0, 0), inv_g(0, 1);
  // g1m << inv_g(1, 0), inv_g(1, 1);
  // gj0m << dot(xu, xuu), dot(xv, xuu), dot(xu, xuv), dot(xv, xuv);
  // gj1m << dot(xu, xuv), dot(xv, xuv), dot(xu, xvv), dot(xv, xvv);
  // C1.row(0) = C.row(1);
  // C1.row(1) = C.row(2);
  // C1m.row(0) = C.row(3);
  // C1m.row(1) = C.row(4);
  // C2m.row(0) = C.row(4);
  // C2m.row(1) = C.row(5);

  // Eigen::VectorXd CoeffA;
  // Eigen::VectorXd CoeffB;
  // Eigen::VectorXd CoeffC;
  // Eigen::VectorXd CoeffD;

  // CoeffA = g0m.transpose() * (C1m - gj0m * inv_g * C1);
  // CoeffB = g0m.transpose() * (C2m - gj1m * inv_g * C1);
  // CoeffC = g1m.transpose() * (C1m - gj0m * inv_g * C1);
  // CoeffD = g1m.transpose() * (C2m - gj1m * inv_g * C1);

  // for (int i = 0; i < s; ++i) {
  //   for (auto j = 0; j < ring[i].size(); ++j) {
  //     int entry = ring[i][j].first;
  //     auto w = ring[i][j].second;

  //     entries.push_back(T(vid, entry, w * CoeffA(i)));
  //     entries.push_back(T(n + vid, entry, w * CoeffB(i)));
  //     entries.push_back(T(2 * n + vid, entry, w * CoeffC(i)));
  //     entries.push_back(T(3 * n + vid, entry, w * CoeffD(i)));
  //   }
  // }
}
void laplacian_entries(vector<Eigen::Triplet<double>> &entries,
                       const Eigen::Matrix2d &g, const Eigen::Matrix2d &g_inv,
                       const vector<vector<pair<int, float>>> &ring,
                       const Eigen::VectorXd &c, const Eigen::MatrixXd &a) {
  int vid = ring[0][0].first;
  typedef Eigen::Triplet<double> T;
  Eigen::VectorXd b;
  double c0_squared = pow(c[0], 2);
  double c1_squared = pow(c[1], 2);

  double det = 1 + c0_squared + c1_squared;

  double g11u = 2 * c[0] * c[2], g12u = c[0] * c[3] + c[1] * c[2],
         g22_u = 2 * c[1] * c[3], g11_v = 2 * c[0] * c[3],
         g12_v = c[0] * c[4] + c[1] * c[3], g22_v = 2 * c[1] * c[4],
         g_u = g11u + g22_u, g_v = g11_v + g22_v,
         g_invu11 = (det * g22_u - g(1, 1) * g_u) / pow(det, 2),
         g_invu12 = -(det * g12u - g(0, 1) * g_u) / pow(det, 2),
         g_invv21 = -(det * g12_v - g(0, 1) * g_v) / pow(det, 2),
         g_invv22 = (det * g11_v - g(0, 0) * g_v) / pow(det, 2);

  double coeff0 = (g_u * g_inv(0, 0) + g_v * g_inv(1, 0)) / (2 * det) +
                  g_invu11 + g_invv21,
         coeff1 = (g_u * g_inv(0, 1) + g_v * g_inv(1, 1)) / (2 * det) +
                  g_invu12 + g_invv22,
         coeff2 = g_inv(0, 0), coeff3 = 2 * g_inv(0, 1), coeff4 = g_inv(1, 1);

  b = coeff0 * a.row(0) + coeff1 * a.row(1) + coeff2 * a.row(2) +
      coeff3 * a.row(3) + coeff4 * a.row(4);
  for (int i = 0; i < b.size(); ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;

      entries.push_back(T(vid, entry, b(i) * ring[i][j].second));
    }
  }
}
void laplacian_entries(vector<Eigen::Triplet<double>> &entries,
                       const Eigen::Matrix2d &g, const float &det,
                       const vector<vector<pair<int, float>>> &ring,
                       const Eigen::MatrixXd &C, const vec3f &xu,
                       const vec3f &xv, const vec3f &xuu, const vec3f &xuv,
                       const vec3f &xvv) {
  int vid = ring[0][0].first;
  typedef Eigen::Triplet<double> T;
  Eigen::MatrixXd Coeff(5, ring.size());
  Coeff.row(0) = C.row(1);
  Coeff.row(1) = C.row(2);
  Coeff.row(2) = C.row(3);
  Coeff.row(3) = C.row(4);
  Coeff.row(4) = C.row(5);
  auto g_deltau =
      -(g(0, 0) * (g(1, 1) * dot(xu, xvv) - g(0, 1) * dot(xv, xvv)) +
        2 * g(0, 1) * (g(0, 1) * dot(xv, xuv) - g(1, 1) * dot(xu, xuv)) +
        g(1, 1) * (g(1, 1) * dot(xu, xuu) - g(0, 1) * dot(xv, xuu))) /
      pow(det, 2);
  auto g_deltav =
      -(g(0, 0) * (g(0, 0) * dot(xv, xvv) - g(0, 1) * dot(xu, xvv)) +
        2 * g(0, 1) * (g(0, 1) * dot(xu, xuv) - g(0, 0) * dot(xv, xuv)) +
        g(1, 1) * (g(0, 0) * dot(xv, xuu) - g(0, 1) * dot(xu, xuu))) /
      pow(det, 2);
  auto g_deltauu = g(1, 1) / det;
  auto g_deltauv = -2 * g(0, 1) / det;
  auto g_deltavv = g(0, 0) / det;
  Eigen::VectorXd w(5);
  w << g_deltau, g_deltav, g_deltauu, g_deltauv, g_deltavv;
  Eigen::VectorXd b = w.transpose() * Coeff;
  for (int i = 0; i < b.rows(); ++i) {
    for (auto j = 0; j < ring[i].size(); ++j) {
      int entry = ring[i][j].first;

      entries.push_back(T(vid, entry, b(i) * ring[i][j].second));
    }
  }
}
vector<float> subdivide_angles(const int number_of_subdivision) {
  auto tetas = vector<float>(number_of_subdivision);
  auto step = 2 * M_PI / number_of_subdivision;
  auto phi = 0.f;
  for (auto i = 0; i < number_of_subdivision; ++i) {
    tetas[i] = phi;
    phi += step;
  }
  return tetas;
}

vector<float> angle_at_vertices(const vector<vec3i> &triangles,
                                const vector<vec3f> &positions,
                                const vector<vector<int>> &v2t, const int vid) {
  auto star = v2t[vid];
  auto result = vector<float>(star.size());
  auto vert_pos = positions[vid];
  auto accumulated = 0.f;
  for (auto i = 0; i < star.size(); ++i) {
    auto tid = star[i];
    auto k = find(triangles[tid], vid);
    auto vid0 = triangles[tid][(k + 1) % 3];
    auto vid1 = triangles[tid][(k + 2) % 3];
    accumulated +=
        angle(positions[vid0] - vert_pos, positions[vid1] - vert_pos);
    result[i] = accumulated;
  }

  return result;
}
std::tuple<vector<mesh_point>, vector<float>, vector<float>>
uniform_stencil(const shape_data &data, const shape_geometry &geometry,
                const int vid, const int number_of_samples) {
  vector<vec2f> directions(number_of_samples);
  vector<int> tris(number_of_samples);
  vector<mesh_point> samples(number_of_samples);
  auto len =
      nbr_avg_edge_length(data.triangles, data.positions, geometry.v2t, vid);

  auto lens = vector<float>(number_of_samples, len);
  auto tetas = subdivide_angles(number_of_samples);
  auto step = geometry.total_angles[vid] / number_of_samples;
  auto steps =
      angle_at_vertices(data.triangles, data.positions, geometry.v2t, vid);
  auto accumulated = 0.f;
  auto alpha = 0.f;
  auto start = point_from_vert(data.triangles, geometry.v2t, vid);
  auto curr_tid = geometry.v2t[vid][0];
  auto vert_pos = data.positions[vid];
  auto tid_count = 0;
  for (auto i = 0; i < number_of_samples; ++i) {
    auto start = point_from_vert(data.triangles, geometry.v2t, vid, tid_count);
    auto k = find(data.triangles[curr_tid], vid);
    auto vid0 = data.triangles[curr_tid][(k + 1) % 3];
    auto vid1 = data.triangles[curr_tid][(k + 2) % 3];
    auto curr_angle =
        angle(data.positions[vid0] - vert_pos, data.positions[vid1] - vert_pos);
    auto offset = std::abs(curr_angle - alpha);
    if (offset < 0.02)
      alpha = curr_angle;

    auto v =
        rot_vect(data.positions[vid0] - vert_pos,
                 tid_normal(data.triangles, data.positions, curr_tid), alpha);
    auto path = straightest_path_from_vert_in_triangle(data, geometry, curr_tid,
                                                       vid, len, normalize(v));
    // straightest_geodesic(data.triangles, data.positions, data.normals,
    //                      geometry.adjacencies, geometry.v2t,
    //                      geometry.angles, geometry.total_angles, start,
    //                      normalize(v), len, false);

    samples[i] = path.back();
    accumulated += step;
    alpha += step;

    while (alpha > curr_angle) {
      alpha = accumulated - steps[tid_count];
      curr_tid = geometry.adjacencies[curr_tid][(k + 2) % 3];
      k = find(data.triangles[curr_tid], vid);
      vid0 = data.triangles[curr_tid][(k + 1) % 3];
      vid1 = data.triangles[curr_tid][(k + 2) % 3];
      curr_angle = angle(data.positions[vid0] - vert_pos,
                         data.positions[vid1] - vert_pos);
      ++tid_count;
    }
    // if (steps[tid_count] < accumulated) {
    //   while (steps[tid_count] < accumulated) {
    //     alpha = accumulated - steps[tid_count];
    //     curr_tid = geometry.adjacencies[curr_tid][(k + 2) % 3];
    //     ++tid_count;
    //   }
    // } else {
    //   alpha += step;
    // }
  }
  return {samples, lens, tetas};
}
vec3f flip_bary_to_adjacent_tri(const vector<vec3i> &adjacencies,
                                const int tid0, const int tid1,
                                const vec3f &bary) {
  if (tid0 == tid1)
    return bary;
  auto new_bary = zero3f;
  auto k1 = find(adjacencies[tid1], tid0);
  auto k0 = find(adjacencies[tid0], tid1);
  if (k1 == -1) {
    std::cout << "Error, faces are not adjacent" << std::endl;
    return zero3f;
  }
  new_bary[k1] = bary[(k0 + 1) % 3];
  new_bary[(k1 + 1) % 3] = bary[k0];
  new_bary[(k1 + 2) % 3] = bary[(k0 + 2) % 3];

  return new_bary;
}
std::tuple<vector<int>, mesh_point, mesh_point>
handle_short_strips(const vector<vec3i> &triangles,
                    const vector<vec3f> &positions, const vector<int> &strip,
                    const mesh_point &start, const mesh_point &end) {
  if (strip.size() == 1) {
    return {strip, start, end};
  } else if (strip.size() == 2) {
    auto [inside, b2f] =
        point_in_triangle(triangles, positions, start.face,
                          eval_position(triangles, positions, end));
    if (inside) {
      auto new_end = mesh_point{start.face, b2f};
      return {{start.face}, start, new_end};
    }
    std::tie(inside, b2f) =
        point_in_triangle(triangles, positions, end.face,
                          eval_position(triangles, positions, start));

    if (inside) {
      auto new_start = mesh_point{end.face, b2f};
      return {{end.face}, new_start, end};
    }

    return {strip, start, end};
  }
  return {{-1}, {}, {}};
}
std::tuple<vector<int>, mesh_point, mesh_point>
cleaned_strip(const vector<vec3i> &triangles, const vector<vec3f> &positions,
              const vector<vec3i> &adjacencies, const vector<int> &strip,
              const mesh_point &start, const mesh_point &end) {
  vector<int> cleaned = strip;

  auto start_entry = 0, end_entry = (int)strip.size() - 1;
  auto b3f = zero3f;
  auto new_start = start;
  auto new_end = end;
  auto [is_vert, kv] = point_is_vert(end);
  auto [is_edge, ke] = point_is_edge(end);
  if (strip.size() <= 2)
    return handle_short_strips(triangles, positions, strip, start, end);
  // Erasing from the bottom
  if (is_vert) {
    auto vid = triangles[end.face][kv];
    auto curr_tid = strip[end_entry - 1];
    kv = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.pop_back();
      --end_entry;
      if (end_entry == 1)
        break;
      // see comment below
      auto curr_tid = strip[end_entry - 1];
      kv = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned.back()], vid);
    assert(kv != -1);
    b3f[kv] = 1;
    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}}; // updating end
  } else if (is_edge) {
    if (end.face != strip.back()) {
      assert(adjacencies[end.face][ke] == strip.back());

      if (end.face == strip[end_entry - 1])
        cleaned.pop_back();
    } else if (adjacencies[end.face][ke] == strip[end_entry - 1])
      cleaned.pop_back();

    b3f = flip_bary_to_adjacent_tri(adjacencies, end.face, cleaned.back(),
                                    get_bary(end.uv));

    new_end = mesh_point{cleaned.back(), vec2f{b3f.y, b3f.z}}; // updating end
  }
  std::tie(is_vert, kv) = point_is_vert(start);
  std::tie(is_edge, ke) = point_is_vert(start);

  if (is_vert) {
    auto vid = triangles[start.face][kv];
    auto curr_tid = strip[start_entry + 1];
    kv = find(triangles[curr_tid], vid);
    while (kv != -1) {
      cleaned.erase(cleaned.begin());
      ++start_entry;
      if (start_entry > end_entry - 1)
        break;
      auto curr_tid = strip[start_entry + 1];
      kv = find(triangles[curr_tid], vid);
    }
    kv = find(triangles[cleaned[0]], vid);
    assert(kv != -1);
    b3f = zero3f;
    b3f[kv] = 1;
    new_start = mesh_point{cleaned[0], vec2f{b3f.y, b3f.z}}; // udpdating start

  } else if (is_edge) {
    if (start.face != strip[0]) {
      assert(adjacencies[start.face][ke] == strip[0]);
      if (start.face == strip[1])
        cleaned.erase(cleaned.begin());
    } else if (adjacencies[start.face][ke] == strip[1]) {
      cleaned.erase(cleaned.begin());
    }
    b3f = flip_bary_to_adjacent_tri(adjacencies, start.face, cleaned[0],
                                    get_bary(start.uv));
    new_start = {cleaned[0], vec2f{b3f.y, b3f.z}}; // updating start
  }
  return {cleaned, new_start, new_end};
}
template <typename Update, typename Stop, typename Exit>
void search_strip(vector<float> &field, vector<bool> &in_queue,
                  const dual_geodesic_solver &solver,
                  const vector<vec3i> &triangles,
                  const vector<vec3f> &positions, int start, int end,
                  Update &&update, Stop &&stop, Exit &&exit) {
  auto destination_pos =
      eval_position(triangles, positions, {end, {1.0f / 3, 1.0f / 3}});

  auto estimate_dist = [&](int face) {
    auto p = eval_position(triangles, positions, {face, {1.0f / 3, 1.0f / 3}});
    return length(p - destination_pos);
  };
  field[start] = estimate_dist(start);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

  while (!queue.empty()) {
    auto node = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight)
        break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node))
      break;
    if (stop(node))
      continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1)
        continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance)
        continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      if (update(node, neighbor, new_distance))
        return;
    }
  }
}
vector<int> compute_strip_tlv(const shape_data &data,
                              const dual_geodesic_solver &solver, int start,
                              int end) {
  if (start == end)
    return {start};

  thread_local static auto parents = vector<int>{};
  thread_local static auto field = vector<float>{};
  thread_local static auto in_queue = vector<bool>{};

  if (parents.size() != solver.graph.size()) {
    parents.assign(solver.graph.size(), -1);
    field.assign(solver.graph.size(), flt_max);
    in_queue.assign(solver.graph.size(), false);
  }

  // initialize once for all and sparsely cleanup at the end of every solve
  auto visited = vector<int>{start};
  auto sources = vector<int>{start};
  auto update = [&visited, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(field, in_queue, solver, data.triangles, data.positions, start,
               end, update, stop, exit);

  // extract_strip
  auto strip = vector<int>{};
  auto node = end;
  strip.reserve((int)sqrt(parents.size()));
  while (node != -1) {
    assert(find(strip, node) != 1);
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  for (auto &v : visited) {
    parents[v] = -1;
    field[v] = flt_max;
    in_queue[v] = false;
  }
  // assert(check_strip(geometry.adjacencies, strip));
  return strip;
}

geodesic_path compute_geodesic_path(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const dual_geodesic_solver &solver,
                                    const mesh_point &start,
                                    const mesh_point &end) {
  // profile_function();

  vector<int> parents;
  auto strip = compute_strip_tlv(data, solver, end.face, start.face);
  // get_strip(mesh.solver, data.triangles, data.positions,
  // geometry.adjacencies,
  //           geometry.v2t, data.angles, end, start);
  auto path = geodesic_path{};
  auto [cleaned, new_start, new_end] = cleaned_strip(
      data.triangles, data.positions, geometry.adjacencies, strip, start, end);

  if (new_start.face == new_end.face) {
    path.start = new_start;
    path.end = new_end;
    path.strip = {new_start.face};
    return path;
  }

  path = shortest_path(data.triangles, data.positions, geometry.adjacencies,
                       new_start, new_end, cleaned);
  return path;
}
std::tuple<vector<mesh_point>, vector<float>, vector<float>>
ring_stencil(const shape_data &data, const shape_geometry &geometry,
             const dual_geodesic_solver &solver, const int vid) {

  // vector<mesh_point> start(number_of_samples);
  bool use_the_1_ring = true;
  auto k = 1;
  auto nbr = k_ring(data.triangles, geometry.v2t, vid, k);
  auto lens = vector<float>{};
  auto tetas = vector<float>{};
  auto samples = vector<mesh_point>{};
  while (nbr.size() < 5) {
    use_the_1_ring = false;
    ++k;
    nbr = k_ring(data.triangles, geometry.v2t, vid, k);
  }

  lens.resize(nbr.size());
  tetas.resize(nbr.size());
  samples.resize(nbr.size());
  auto ring = one_ring(data.triangles, geometry.v2t, vid);
  auto point_vid = point_from_vert(data.triangles, geometry.v2t, vid);
  // auto e = polar_basis(data.triangles, geometry.v2t, positions, normals,
  // vid);
  for (auto i = 0; i < nbr.size(); ++i) {
    if (use_the_1_ring) {
      samples[i] = point_from_vert(data.triangles, geometry.v2t, nbr[i]);
      tetas[i] = geometry.angles[vid][i];
      lens[i] = length(data.positions[vid] - data.positions[nbr[i]]);
    } else if (i < ring.size()) {
      samples[i] = point_from_vert(data.triangles, geometry.v2t, nbr[i]);
      tetas[i] = geometry.angles[vid][i];
      lens[i] = length(data.positions[vid] - data.positions[nbr[i]]);
    } else {
      samples[i] = point_from_vert(data.triangles, geometry.v2t, nbr[i]);
      auto path =
          compute_geodesic_path(data, geometry, solver, point_vid, samples[i]);
      auto pos = path_positions(path, data.triangles, data.positions,
                                geometry.adjacencies);
      auto dir = pos[1] - data.positions[vid];
      auto k = find(data.triangles[path.start.face], vid);
      auto theta = angle(
          dir, data.positions[data.triangles[path.start.face][(k + 1) % 3]] -
                   data.positions[vid]);
      auto it = find(nbr.begin(), nbr.end(),
                     data.triangles[path.start.face][(k + 1) % 3]);
      auto entry = std::distance(nbr.begin(), it);
      theta = geometry.angles[vid][entry] +
              theta * 2 * pif / geometry.total_angles[vid];
      tetas[i] = theta;
      lens[i] = path_length(pos);
    }
  }
  return {samples, lens, tetas};
}
// given a vertex v and an angle theta, returns the direction in 2D of the
// line that forms an angle theta at vid
std::pair<mesh_point, float>
path_length_for_stencil(const vector<vec3i> &triangles,
                        const vector<vec3f> &positions,
                        const vector<vec3i> &adjacencies, const int tid,
                        const int vid, const int it, const int n) {
  auto p = mesh_point{};
  auto k = find(triangles[tid], vid);
  auto flat_tid = init_flat_triangle(positions, triangles[tid]);

  auto v0 = flat_tid[k];
  auto v1 = flat_tid[(k + 1) % 3];
  auto v2 = flat_tid[(k + 2) % 3];
  auto phi12 = angle(v1 - v0, v2 - v0);
  auto phi02 = angle(v0 - v1, v2 - v1);
  auto phi01 = angle(v1 - v2, v0 - v2);
  auto l01 = length(v1 - v0);
  auto l02 = length(v0 - v2);
  auto l12 = length(v1 - v2);
  auto theta = phi12 * it / n;

  auto y0 = l12 * std::sin(phi01) * std::sin(theta) /
            (std::sin(phi02) * std::sin(phi12 - theta) +
             std::sin(phi01) * std::sin(theta));
  auto x0 = std::sin(phi02) * y0 / std::sin(theta);
  auto next = adjacencies[tid][(k + 1) % 3];
  auto flat_n = unfold_face(triangles, positions, flat_tid, tid, next);
  auto h = find(adjacencies[next], tid);
  auto v3 = flat_n[(h + 2) % 3];
  auto phi23 = angle(v2 - v1, v3 - v1);
  auto l13 = length(v3 - v1);
  auto xi = pif - theta - phi02 - phi23;
  auto x1 = y0 * std::sin(phi23) / std::sin(xi);

  auto y1 = l01 * std::sin(theta) / std::sin(xi);
  if (y1 <= l13) {
    auto bary = zero3f;
    bary[(h + 1) % 3] = 1 - y1 / l13;
    bary[(h + 2) % 3] = y1 / l13;

    p = {next, vec2f{bary.y, bary.z}};
  } else {
    auto phi13 = angle(v1 - v2, v3 - v2);
    xi = pif - (phi12 - theta + phi01 + phi13);
    x1 = (l12 - y0) * std::sin(phi13) / std::sin(xi);
    y1 = l02 * std::sin(phi12 - theta) / std::sin(xi);
    auto bary = zero3f;
    bary[h] = 1 - y1 / l13;
    bary[(h + 2) % 3] = y1 / l13;
    p = {next, vec2f{bary.y, bary.z}};
  }

  return {p, x0 + x1};
}
std::tuple<vector<mesh_point>, vector<float>, vector<float>> weighted_stencil(
    const vector<vec3i> &triangles, const vector<vec3f> &positions,
    const geodesic_solver &solver, const vector<vec3f> &normals,
    const vector<vec3i> &adjacencies, const vector<vector<int>> &v2t,
    const vector<vector<float>> &angles, const vector<float> &total_angles,
    const int vid, const int n) {

  auto start = point_from_vert(triangles, v2t, vid);
  auto e = polar_basis(triangles, positions, v2t, normals, vid);
  auto nbr = solver.graph[vid];
  auto s = (int)nbr.size() / 2 * (n + 1);
  auto tetas = vector<float>(s);
  vector<vec2f> directions(s);
  vector<float> lens(s);
  vector<mesh_point> samples(s);
  for (auto j = 0; j < nbr.size(); ++j) {
    if (j % 2)
      continue;
    auto step = (j != nbr.size() - 2)
                    ? (angles[vid][j + 2] - angles[vid][j]) / n
                    : (2 * pif - angles[vid][j]) / n;
    for (auto i = 0; i < n; ++i) {
      auto entry = n * j / 2 + i;
      tetas[entry] = angles[vid][j] + i * step;

      auto curr_dir = rot_vect(e, normals[vid], tetas[entry]);
      if (i == 0) {
        samples[entry] = point_from_vert(triangles, v2t, nbr[j].node);
        lens[entry] = nbr[j].length;
      } else {

        auto tid = v2t[vid][j / 2];
        std::tie(samples[entry], lens[entry]) = path_length_for_stencil(
            triangles, positions, adjacencies, tid, vid, i, n);
      }
    }
  }
  return {samples, lens, tetas};
}
// we follow the order 00 01 10 11
vector<vec2f> Christoffel_symbol(const Eigen::Matrix2d &g_inv,
                                 const Eigen::VectorXd &c) {
  vector<vec2f> Christoffel(4);
  double g11u = 2 * c[0] * c[2], g12u = c[0] * c[3] + c[1] * c[2], g21u = g12u,
         g22u = 2 * c[1] * c[3], g11v = 2 * c[0] * c[3],
         g12v = c[0] * c[4] + c[1] * c[3], g21v = g12v, g22v = 2 * c[1] * c[4];

  Christoffel[0].x = 0.5 * (g_inv(0, 0) * g11u - g_inv(0, 1) * g11v);
  Christoffel[0].y =
      0.5 * (g_inv(1, 0) * g11u + g_inv(1, 1) * (2 * g12u - g11v));

  Christoffel[1].x = 0.5 * (g_inv(0, 0) * g11v + g_inv(0, 1) * g22u);
  Christoffel[1].y = 0.5 * (g_inv(1, 0) * g11v + g_inv(1, 1) * g22u);

  Christoffel[2] = Christoffel[1];

  Christoffel[3].x = 0.5 * (g_inv(0, 0) * g22u + g_inv(0, 1) * g22v);
  Christoffel[3].y =
      0.5 * (g_inv(1, 0) * (2 * g21v - g22u) + g_inv(1, 1) * g22v);

  return Christoffel;
}
Eigen::MatrixXd Moore_Penrose_inverse(const Eigen::MatrixXd &C,
                                      const double &epsilon = 1e-8) {

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(C, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  double cond = svd.singularValues()(0) /
                svd.singularValues()(svd.singularValues().size() - 1);
  if (cond > 50) {
    double tol = epsilon * std::max(C.cols(), C.rows()) *
                 svd.singularValues().array().abs()(0);
    return svd.matrixV() *
           (svd.singularValues().array().abs() > tol)
               .select(svd.singularValues().array().inverse(), 0)
               .matrix()
               .asDiagonal() *
           svd.matrixU().adjoint();
  } else
    return C.inverse();
}
float greatest_bary_coord(const vec3f &bary) {
  auto max_b = flt_min;
  for (auto i = 0; i < 3; ++i) {
    max_b = std::max(max_b, bary[i]);
  }
  return max_b;
}
shape_op init_discrete_diff_op_xu(const shape_data &data,
                                  const shape_geometry &geometry,
                                  const dual_geodesic_solver &solver,
                                  const bool grad, const bool lap,
                                  const bool hess) {
  time_function();
  shape_op op = {};
  typedef Eigen::Triplet<double> T;
  vector<T> L_entries;
  vector<T> G_entries;
  vector<T> H_entries;
  int V = (int)data.positions.size();
  op.quadrics.resize(V);
  op.CMat.resize(V);
  op.Rhs.resize(V);

  Eigen::SparseMatrix<double> Grad;
  Eigen::SparseMatrix<double> Lap;
  vector<Eigen::MatrixXd> A_coeff(V);
  Eigen::SparseMatrix<double> Hess;
  vector<Eigen::Matrix2d> g_invs(V);

  vector<vector<vector<pair<int, float>>>> contributes_map(V);
  vector<float> gaussian_curv(data.positions.size());
  float min_curv = flt_min;
  int invertible = 0;
  auto [Verts, Tris] = libigl_wrapper(data.positions, data.triangles);
  for (int i = 0; i < V; ++i) {
    if (geometry.angles[i].size() == 0)
      continue;
    // auto [nbr, lens, tetas] = ring_stencil(data, geometry, solver, i);
    auto [nbr, lens, tetas] = uniform_stencil(data, geometry, i, 36);
    vec3f vert = data.positions[i];
    vec3f n = data.normals[i];
    int s = (int)nbr.size();
    Eigen::MatrixXd A(s + 1, 6);
    Eigen::MatrixXd P(s + 1, 3);
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> W(s + 1);
    W.diagonal()[0] = 1;
    for (auto h = 0; h < 3; ++h) {
      P(0, h) = vert[h];
    }
    auto teta = 0.f;
    auto pos = zero2f;
    auto d = 0.f;
    auto coords = zero3f;
    vector<vector<pair<int, float>>> a_map(s + 1);
    a_map[0].push_back({i, 1});
    A(0, 0) = 1;
    A(0, 1) = 0;
    A(0, 2) = 0;
    A(0, 3) = 0;
    A(0, 4) = 0;
    A(0, 5) = 0;

    for (int j = 0; j < s; ++j) {
      pos =
          vec2f{lens[j] * yocto::cos(tetas[j]), lens[j] * yocto::sin(tetas[j])};
      coords = eval_position(data.triangles, data.positions, nbr[j]);
      auto bary_end = get_bary(nbr[j].uv);
      for (auto h = 0; h < 3; ++h) {
        a_map[j + 1].push_back({data.triangles[nbr[j].face][h], bary_end[h]});
        P(j + 1, h) = coords[h];
      }
      A(j + 1, 0) = 1;
      A(j + 1, 1) = pos[0];
      A(j + 1, 2) = pos[1];
      A(j + 1, 3) = pow(pos[0], 2);
      A(j + 1, 4) = pos[0] * pos[1];
      A(j + 1, 5) = pow(pos[1], 2);
      W.diagonal()[j + 1] = greatest_bary_coord(bary_end);
    }

    Eigen::MatrixXd At = Eigen::Transpose<Eigen::MatrixXd>(A);
    Eigen::MatrixXd B = At * W * A;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(B);
    Eigen::MatrixXd c(3, 6);
    // Eigen::MatrixXd C;
    // Eigen::JacobiSVD<Eigen::MatrixXd> svd(B, Eigen::ComputeFullU |
    //                                              Eigen::ComputeFullV);
    // double cond = svd.singularValues()(0) /
    //               svd.singularValues()(svd.singularValues().size() - 1);
    // if (cond > 1e-3) {

    //   double tol = 1e-8;
    //   C = svd.matrixV() *
    //       (svd.singularValues().array().abs() > tol)
    //           .select(svd.singularValues().array().inverse(), 0)
    //           .matrix()
    //           .asDiagonal() *
    //       svd.matrixU().adjoint();
    // } else {
    //   ++invertible;
    //   C = B.inverse();
    // }
    // C = C * At;

    for (auto h = 0; h < 3; ++h) {
      c.row(h) = dec.solve(At * W * P.col(h));
    }

    // if (dec.isInvertible()) {
    //   Eigen::MatrixXd inv = B.inverse();
    //   C = inv * At;
    //   op.CMat[i] = C;
    //   for (auto h = 0; h < 3; ++h) {
    //     c.row(h) = C * P.col(h);
    //   }
    //   ++invertible;
    // } else {
    //   Eigen::MatrixXd Bi = Moore_Penrose_inverse(B);
    //   C = Bi * At;
    //   op.CMat[i] = C;
    //   for (auto h = 0; h < 3; ++h) {
    //     c.row(h) = dec.solve(At * P.col(h));
    //   }
    // }
    op.quadrics[i] = c;
    op.CMat[i] = dec;
    op.Rhs[i] = At * W;
    Eigen::Matrix2d g;
    Eigen::Matrix2d g_inv;
    auto xu = vec3f{(float)c(0, 1), (float)c(1, 1), (float)c(2, 1)};
    auto xv = vec3f{(float)c(0, 2), (float)c(1, 2), (float)c(2, 2)};
    auto xuu =
        vec3f{2 * (float)c(0, 3), 2 * (float)c(1, 3), 2 * (float)c(2, 3)};
    auto xuv = vec3f{(float)c(0, 4), (float)c(1, 4), (float)c(2, 4)};
    auto xvv =
        vec3f{2 * (float)c(0, 5), 2 * (float)c(1, 5), 2 * (float)c(2, 5)};
    g << dot(xu, xu), dot(xu, xv), dot(xu, xv), dot(xv, xv);
    auto det = g.determinant();
    g_inv = g.inverse();

    // if (lap)
    //   laplacian_entries(L_entries, g, det, a_map, C, xu, xv, xuu, xuv, xvv);

    // if (grad)
    //   fill_riemannian_gradient_entries(G_entries, a_map, g, det, C.row(1),
    //                                    C.row(2), xu, xv, V);
    // if (hess)
    //   fill_hessian_entries(H_entries, g_inv, a_map, C, xu, xv, xuu, xuv, xvv,
    //                        V);
  }

  if (lap) {
    op.Lap.resize(V, V);
    op.Lap.setFromTriplets(L_entries.begin(), L_entries.end());
  }
  if (grad) {
    op.Grad.resize(3 * V, V);
    op.Grad.setFromTriplets(G_entries.begin(), G_entries.end());
  }

  if (hess) {
    op.Hess.resize(4 * V, V);
    op.Hess.setFromTriplets(H_entries.begin(), H_entries.end());
  }

  std::cout << "\n";
  std::printf("%f %% of the vertices had invertible matrices",
              (float)invertible / V * 100);
  return op;
}

vector<pair<int, double>> cotangent_weights(const shape_data &data,
                                            const shape_geometry &topology,
                                            const int vid) {
  auto star = topology.v2t[vid];
  auto s = (int)star.size();
  auto weights = vector<pair<int, double>>(s);
  for (auto i = 0; i < s; ++i) {
    auto sum = 0.f;
    auto tid0 = star[(s - 1 + i) % s];
    auto tid1 = star[i];
    auto k = find(data.triangles[tid1], vid);
    auto neighbor = data.triangles[tid1][(k + 1) % 3];
    auto opp = opposite_vertex(data.triangles[tid1],
                               vec2i{vid, data.triangles[tid1][(k + 2) % 3]});
    auto alpha0 =
        angle_at_vertex(data, tid1, data.triangles[tid1][(k + 2) % 3]);

    auto alpha1 = angle_at_vertex(data, tid0, opp);
    auto c0 = cot(alpha0);
    auto c1 = cot(alpha1);
    auto count = 0.f;
    if (!std::isnan(c0)) {
      sum += std::max(1e-10, c0);
      count += 1.0;
    }
    if (!std::isnan(c1)) {
      sum += std::max(1e-10, c1);
      count += 1.0;
    }
    weights[i] = (count == 0.f) ? std::make_pair(neighbor, 0.f)
                                : std::make_pair(neighbor, sum / count);
  }
  return weights;
}
vector<std::pair<vec2i, double>> vert_mass(const shape_data &data,
                                           const shape_geometry &topology,
                                           const int vid, const bool lumped) {
  auto star = topology.v2t[vid];
  auto contributes = vector<std::pair<vec2i, double>>{};
  if (lumped) {
    auto mass = 0.f;
    for (auto i = 0; i < star.size(); ++i) {
      mass += tid_area(data, star[i]);
    }
    contributes.push_back(std::make_pair(vec2i{vid, vid}, mass / 3.f));
  } else {
    auto mass = 0.f;
    auto s = (int)star.size();
    for (auto i = 0; i < s; ++i) {
      auto curr_area = tid_area(data, star[i]);
      mass += curr_area / 6.f;
      auto k = find(data.triangles[star[i]], vid);
      auto vid0 = data.triangles[star[i]][(k + 1) % 3];
      if (vid0 < vid)
        continue;
      auto tid = star[(s - 1 + i) % s];
      auto tot_contribute = (tid_area(data, tid) + curr_area) / 12.f;
      contributes.push_back(std::make_pair(vec2i{vid, vid0}, tot_contribute));
    }
    contributes.push_back(std::make_pair(vec2i{vid, vid}, mass));
  }

  return contributes;
}
Eigen::SparseMatrix<double>
cot_lap_stiffness_matrix(const shape_data &data,
                         const shape_geometry &topology) {
  Eigen::SparseMatrix<double> L(data.positions.size(), data.positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T> entries;
  for (auto i = 0; i < data.positions.size(); ++i) {
    auto wgts = cotangent_weights(data, topology, i);
    auto sum = 0.f;
    for (auto weight : wgts) {
      entries.push_back(T(i, weight.first, weight.second));
      sum -= weight.second;
    }
    if (sum == 0) {
      std::cerr << "Error: null row in the matrix! (disconnected vertex? I "
                   "put 1 in the diagonal)"
                << std::endl;
      sum = 1.0;
    }
    entries.push_back(T(i, i, sum));
  }

  L.setFromTriplets(entries.begin(), entries.end());

  return L;
}

Eigen::SparseMatrix<double> mass_matrix(const shape_data &data,
                                        const shape_geometry &topology,
                                        const bool lumped) {
  Eigen::SparseMatrix<double> M(data.positions.size(), data.positions.size());
  typedef Eigen::Triplet<double> T;
  vector<T> entries;
  if (lumped) {
    for (auto i = 0; i < data.positions.size(); ++i) {
      auto mass = vert_mass(data, topology, i, lumped);
      entries.push_back(T(i, i, mass[0].second));
    }
  } else {
    for (auto i = 0; i < data.positions.size(); ++i) {
      auto mass = vert_mass(data, topology, i, lumped);
      for (auto contribute : mass) {
        entries.push_back(
            T(contribute.first.x, contribute.first.y, contribute.second));
        if (contribute.first.x != contribute.first.y)
          entries.push_back(
              T(contribute.first.y, contribute.first.x, contribute.second));
      }
    }
  }

  M.setFromTriplets(entries.begin(), entries.end());

  return M;
}

vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const Eigen::VectorXd &f, bool normalized) {
  auto V = (int)positions.size();

  auto F = (int)triangles.size();

  if (G.rows() == 3 * V) {
    vector<vec3f> g(V);

    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < V; ++i) {
      g[i] = vec3f{(float)Grad(i), (float)Grad(V + i), (float)Grad(2 * V + i)};
      // g[i] = polar_to_cartesian(solver, positions, normals, Grad(i),
      //                           Grad(V + i), i);
      // if (normalized)
      //   g[i] = normalize(g[i]);
    }
    return g;
  } else {
    std::printf(
        "Wong size of Gradient matrix (size is %dx%d while it should be "
        "%dx%d)",
        G.rows(), G.cols(), 2 * V, V);
    vector<vec3f> g(F);
    Eigen::VectorXd Grad = G * f;

    for (int i = 0; i < F; ++i) {
      g[i].x = Grad(3 * i);
      g[i].y = Grad(3 * i + 1);
      g[i].z = Grad(3 * i + 2);

      if (normalized)
        g[i] = normalize(g[i]);
    }
    return g;
  }
}
vector<vec3f> compute_grad(const geodesic_solver &solver,
                           const vector<vec3i> &triangles,
                           const vector<vec3f> &positions,
                           const vector<vec3f> &normals,
                           const Eigen::SparseMatrix<double, 1> &G,
                           const vector<float> &f, bool normalized) {
  auto F = wrapper(f);
  return compute_grad(solver, triangles, positions, normals, G, F, normalized);
}
vector<float> compute_laplacian(const shape_op &operators,
                                const vector<float> &distances) {

  auto DIST = wrapper(distances);
  vector<float> laplacian(distances.size());
  Eigen::VectorXd Lap = operators.Lap * DIST;
  for (auto i = 0; i < Lap.size(); ++i) {
    laplacian[i] = Lap(i);
  }
  return laplacian;
}

pair<float, int> optimal_sample(const vector<vec3i> &triangles,
                                const geodesic_path &path, const int vid) {
  auto lerp = path.lerps[0];
  auto tid = path.strip[0];
  auto entry = 0;
  while (lerp == 0 || lerp == 1) {
    ++entry;
    if (path.lerps.size() == entry)
      return {lerp, tid};
    if (find(triangles[path.strip[entry]], vid) != -1) {
      lerp = path.lerps[entry];
      tid = path.strip[entry];
    } else
      return {lerp, tid};
  }

  return {lerp, tid};
}
geodesic_solver extended_solver(const shape_data &data,
                                const dual_geodesic_solver &dual_solver,
                                shape_geometry &geometry, const int k) {
  auto old_solver = make_geodesic_solver(data.triangles, data.positions,
                                         geometry.adjacencies, geometry.v2t);
  geometry.angles = compute_angles_wo_opposite(
      data.triangles, data.positions, geometry.adjacencies, geometry.v2t,
      geometry.total_angles);
  auto avg_valence = 0;
  geodesic_solver solver;
  solver.graph.resize(data.positions.size());
  // data.angles.resize(mesh.positions.size());

  for (auto i = 0; i < data.positions.size(); ++i) {

    auto neighborhood = k_ring(data.triangles, geometry.v2t, i, k);
    auto valence = neighborhood.size();
    // auto new_angles = vector<pair<float, int>>(neighborhood.size());
    auto lengths = unordered_map<int, float>{};
    auto source = point_from_vert(data.triangles, geometry.v2t, i);
    for (auto j = 0; j < neighborhood.size(); ++j) {
      auto curr_point =
          point_from_vert(data.triangles, geometry.v2t, neighborhood[j]);
      auto path = compute_geodesic_path(data, geometry, dual_solver, source,
                                        curr_point);
      if (path.lerps.size() > 0) {
        auto [lerp, tid] = optimal_sample(data.triangles, path, i);
        // new_angles[j] =
        //     std::make_pair(angle_in_tangent_space(
        //                        old_solver, old_angles, mesh.total_angles,
        //                        data.triangles, data.positions, tid, i, lerp),
        //                    neighborhood[j]);
        lengths[neighborhood[j]] = path_length(
            path, data.triangles, data.positions, geometry.adjacencies);
      } else {
        auto entry = node_is_adjacent(old_solver, i, neighborhood[j]);
        if (entry == -1) {
          lengths[neighborhood[j]] = flt_max;
          continue;
          //          std::cout << "Error, the vertices should be adjacent" <<
          //          std::endl; return solver;
        }

        // new_angles[j] = std::make_pair(old_angles[i][entry],
        // neighborhood[j]);
        lengths[neighborhood[j]] = old_solver.graph[i][entry].length;
      }
    }
    avg_valence += valence;
    // sort(new_angles.begin(), new_angles.end());
    // data.angles[i].resize(new_angles.size());
    solver.graph[i].resize(valence);
    for (auto j = 0; j < valence; ++j) {
      solver.graph[i][j] = {neighborhood[j], lengths.at(neighborhood[j])};
    }
  }
  return solver;
  // std::cout << "avg valence is" << std::endl;
  // std::cout << avg_valence / data.positions.size() << std::endl;
}
