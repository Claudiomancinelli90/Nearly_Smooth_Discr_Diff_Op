#include "utilities.h"
std::tuple<vector<float>, geodesic::Mesh, geodesic::GeodesicAlgorithmExact>
exact_geodesic_distance_surazhsky(const vector<vec3i> &triangles,
                                  const vector<vec3f> &positions,
                                  const mesh_point &source) {
  auto [points, faces] = exact_geodesic_wrapper(triangles, positions, source);
  geodesic::Mesh mesh;
  mesh.initialize_mesh_data(points, faces);
  geodesic::GeodesicAlgorithmExact algorithm(&mesh);
  geodesic::SurfacePoint vert_source(&mesh.vertices()[(int)positions.size()]);
  std::vector<geodesic::SurfacePoint> all_sources(1, vert_source);
  double const distance_limit = geodesic::GEODESIC_INF;
  algorithm.propagate(all_sources, distance_limit);
  auto f = vector<float>(positions.size());
  for (unsigned i = 0; i < mesh.vertices().size(); ++i) {
    geodesic::SurfacePoint p(&mesh.vertices()[i]);

    double distance;
    unsigned best_source = algorithm.best_source(
        p, distance); // for a given surface point, find closets source
                      // and distance to this source

    f[i] = distance;
  }

  return {f, mesh, algorithm};
}

vector<vec3f> get_basis_at_vertex(const vector<vec3i> &triangles,
                                  const vector<vec3f> &positions,
                                  const vector<vec3f> &normals,
                                  const vector<vector<int>> &v2t,
                                  const int vid) {

  auto e = polar_basis(triangles, positions, v2t, normals, vid);
  auto n = normals[vid];
  auto e_perp = rot_vect(e, n, pif / 2);
  if (length(cross(e, e_perp) - n) > 1e-10)
    e_perp *= -1;
  return {e, e_perp, n};
}
vector<vec3f> get_basis_at_p(const vector<vec3i> &triangles,
                             const vector<vec3f> &positions, const int tid) {
  auto v0 = positions[triangles[tid][0]];
  auto v1 = positions[triangles[tid][1]];
  auto v2 = positions[triangles[tid][2]];
  auto e = normalize(v1 - v0);
  auto n = triangle_normal(v0, v1, v2);
  auto e_perp = cross(n, e);
  return {e, e_perp, n};
}
vector<vec3f> get_basis_at_vid(const vector<vec3i> &triangles,
                               const vector<vec3f> &positions,
                               const vector<vec3f> &normals,
                               const vector<vector<int>> &v2t, const int vid) {
  auto e0 = polar_basis(triangles, positions, v2t, normals, vid);
  auto n = normals[vid];
  auto e1 = cross(n, e0);
  return {e0, e1, n};
}
// Transformation matrix s.t [x,y,1]=T[x1,y1,1]
Eigen::Matrix3d tranformation_matrix(const vec2f &x1, const vec2f &y1,
                                     const vec2f &O1) {
  Eigen::Matrix3d T = Eigen::Matrix3d::Zero();
  T << x1.x, y1.x, O1.x, x1.y, y1.y, O1.y, 0, 0, 1;

  return T;
}

Eigen::Matrix4d tranformation_matrix(const vec3f &x1, const vec3f &y1,
                                     const vec3f &z1, const vec3f &O1) {
  Eigen::Matrix4d T;

  T << x1.x, x1.y, x1.z, -dot(x1, O1), y1.x, y1.y, y1.z, -dot(y1, O1), z1.x,
      z1.y, z1.z, -dot(z1, O1), 0, 0, 0, 1;

  return T;
}
Eigen::Matrix3d switch_reference_frame(const shape_data &data,
                                       const shape_geometry &geometry,
                                       const int tid0, const int tid1) {
  auto flat_t0 = init_flat_triangle(data.positions, data.triangles[tid0], 0);
  auto flat_t1 =
      unfold_face(data.triangles, data.positions, flat_t0, tid0, tid1);
  auto e0 = normalize(flat_t1[1] - flat_t1[0]);
  auto e1 = vec2f{-e0.y, e0.x};
  return tranformation_matrix(e0, e1, flat_t1[0]);
}
vec3f switch_reference_frame(const Eigen::Matrix4d &T, const vec3f &p) {
  Eigen::Vector4d V;
  V << p.x, p.y, p.z, 1;
  Eigen::Vector4d t = T * V;
  return vec3f{(float)t(0), (float)t(1), (float)t(2)};
}

vec3f switch_reference_frame_vector(const Eigen::Matrix4d &T, const vec3f &p) {
  Eigen::Vector4d V;
  V << p.x, p.y, p.z, 0;
  Eigen::Vector4d t = T * V;
  return vec3f{(float)t(0), (float)t(1), (float)t(2)};
}

vec2f switch_reference_frame(const Eigen::Matrix3d &T, const vec2f &p) {
  Eigen::Vector3d V;
  V << p.x, p.y, 1;
  Eigen::Vector3d t = T * V;
  return vec2f{(float)t(0), (float)t(1)};
}
float angle_in_tg_space(const vector<vec3i> &triangles,
                        const vector<vec3f> &positions,
                        const vector<vec3i> &adjacencies,
                        const vector<vector<float>> &angles,
                        const vector<vector<int>> &v2t, const int vid,
                        const int vid0) {

  auto nbr = one_ring(triangles, adjacencies, v2t, vid);

  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i] == vid0) {
      return angles[vid][i];
    }
  }

  std::cout << "error in computing local coordinates" << std::endl;

  return 0.f;
}
vector<vec2f> get_local_coordinates(const vector<vec3i> &triangles,
                                    const vector<vec3f> &positions,
                                    const vector<vec3i> &adjacencies,
                                    const vector<vector<float>> &angles,
                                    const vector<vector<int>> &v2t,
                                    const int vid, const int vid0,
                                    const int vid1) {
  auto result = vector<vec2f>(2, zero2f);
  auto nbr = one_ring(triangles, adjacencies, v2t, vid);

  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i] == vid0) {
      auto r = length(positions[vid] - positions[vid0]);
      auto theta = angles[vid][i];
      result[0] = vec2f{r * std::cos(theta), r * std::sin(theta)};
    }

    if (nbr[i] == vid1) {
      auto r = length(positions[vid] - positions[vid1]);
      auto theta = angles[vid][i];
      result[1] = vec2f{r * std::cos(theta), r * std::sin(theta)};
    }
  }

  if (length(result[0]) == 0)
    std::cout << "error in computing local coordinates" << std::endl;

  if (length(result[1]) == 0)
    std::cout << "error in computing local coordinates" << std::endl;

  return result;
}
std::tuple<vector<vec2f>, vector<vec2f>>
get_coordinates_for_quadric_fitting(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const int vid, const int tid) {
  auto pos_in_t = vector<vec2f>(6);
  auto pos_tg = vector<vec2f>(6);
  auto k = find(data.triangles[tid], vid);
  auto vid1 = data.triangles[tid][(k + 1) % 3];
  auto vid2 = data.triangles[tid][(k + 2) % 3];

  auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  auto Tp = tranformation_matrix(ep[0], ep[1], ep[2],
                                 data.positions[data.triangles[tid].x]);
  auto p0p = switch_reference_frame(Tp, data.positions[vid]);
  auto p1p = switch_reference_frame(Tp, data.positions[vid1]);
  auto p2p = switch_reference_frame(Tp, data.positions[vid2]);
  auto tri = get_local_coordinates(data.triangles, data.positions,
                                   geometry.adjacencies, geometry.angles,
                                   geometry.v2t, vid, vid1, vid2);
  pos_tg[0] = zero2f;
  pos_tg[1] = tri[0];
  pos_tg[2] = tri[1];

  pos_in_t[0] = {p0p.x, p0p.y};
  pos_in_t[1] = {p1p.x, p1p.y};
  pos_in_t[2] = {p2p.x, p2p.y};

  for (auto i = 0; i < 3; ++i) {
    pos_in_t[3 + i] = 0.5 * (pos_in_t[i] + pos_in_t[(i + 1) % 3]);
    pos_tg[3 + i] = 0.5 * (pos_tg[i] + pos_tg[(i + 1) % 3]);
  }

  return {pos_in_t, pos_tg};
}
float evaluate_quadric(const Eigen::VectorXd &Q, const vec2f &uv) {

  return Q(0) + uv.x * Q(1) + uv.y * Q(2) + Q(3) * std::pow(uv.x, 2) +
         Q(4) * uv.x * uv.y + Q(5) * std::pow(uv.y, 2);
}
vec3f evaluate_quadric(const Eigen::MatrixXd &Q, const vec2f &uv) {

  auto pos = zero3f;
  for (auto j = 0; j < 3; ++j) {
    pos[j] = Q(j, 0) + uv.x * Q(j, 1) + uv.y * Q(j, 2) +
             Q(j, 3) * std::pow(uv.x, 2) + Q(j, 4) * uv.x * uv.y +
             Q(j, 5) * std::pow(uv.y, 2);
  }

  return pos;
}
inline vec3f evaluate_quadric(const shape_op &op, const int vid,
                              const vec2f &uv) {
  auto Q = op.quadrics[vid];
  return evaluate_quadric(Q, uv);
}
float evaluate_quadric_du(const Eigen::VectorXd &Q, const vec2f &uv) {
  return Q(1) + 2 * uv.x * Q(3) + uv.y * Q(4);
}
float evaluate_quadric_dv(const Eigen::VectorXd &Q, const vec2f &uv) {

  return Q(2) + 2 * uv.y * Q(5) + uv.x * Q(4);
}
float evaluate_quadric_duu(const Eigen::VectorXd &Q, const vec2f &uv) {

  return 2 * Q(3);
}
float evaluate_quadric_dvv(const Eigen::VectorXd &Q, const vec2f &uv) {
  return 2 * Q(5);
}
float evaluate_quadric_duv(const Eigen::VectorXd &Q, const vec2f &uv) {

  return Q(4);
}
vec3f evaluate_quadric_du(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xu = zero3f;
  for (auto j = 0; j < 3; ++j) {
    xu[j] = Q(j, 1) + 2 * uv.x * Q(j, 3) + uv.y * Q(j, 4);
  }

  return xu;
}
vec3f evaluate_quadric_dv(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xv = zero3f;
  for (auto j = 0; j < 3; ++j) {
    xv[j] = Q(j, 2) + 2 * uv.y * Q(j, 5) + uv.x * Q(j, 4);
  }

  return xv;
}
vec3f evaluate_quadric_duu(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xuu = zero3f;
  for (auto j = 0; j < 3; ++j) {
    xuu[j] = 2 * Q(j, 3);
  }

  return xuu;
}
vec3f evaluate_quadric_dvv(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xvv = zero3f;
  for (auto j = 0; j < 3; ++j) {
    xvv[j] = 2 * Q(j, 5);
  }

  return xvv;
}
vec3f evaluate_quadric_duv(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xuv = zero3f;
  for (auto j = 0; j < 3; ++j) {
    xuv[j] = Q(j, 4);
  }

  return xuv;
}
vec3f evaluate_quadric_du(const shape_op &op, const int vid, const vec2f &uv) {
  auto Q = op.quadrics[vid];
  return evaluate_quadric_du(Q, uv);
}
vec3f evaluate_quadric_dv(const shape_op &op, const int vid, const vec2f &uv) {
  auto Q = op.quadrics[vid];
  return evaluate_quadric_dv(Q, uv);
}
vec3f evaluate_quadric_n(const Eigen::MatrixXd &Q, const vec2f &uv) {
  auto xu = evaluate_quadric_du(Q, uv);
  auto xv = evaluate_quadric_dv(Q, uv);
  return normalize(cross(xu, xv));
}
vec2f coordinates_in_tangent_space_of_p(const shape_data &data,
                                        const shape_geometry &geometry,
                                        const int vid, const mesh_point &p) {
  auto tid = p.face;
  auto k = find(data.triangles[tid], vid);
  auto vid1 = data.triangles[tid][(k + 1) % 3];
  auto vid2 = data.triangles[tid][(k + 2) % 3];
  auto tri = get_local_coordinates(data.triangles, data.positions,
                                   geometry.adjacencies, geometry.angles,
                                   geometry.v2t, vid, vid1, vid2);
  auto bary = get_bary(p.uv);

  return bary[(k + 1) % 3] * tri[0] + bary[(k + 2) % 3] * tri[1];
}
Eigen::VectorXd reparametrized_quadric_full(const Eigen::VectorXd &Q,
                                            const Eigen::Matrix3d T) {
  Eigen::VectorXd result(Q.rows());

  result(0) = Q(0) + Q(1) * T(0, 2) + Q(2) * T(1, 2) +
              Q(3) * T(0, 2) * T(0, 2) + Q(4) * T(0, 2) * T(1, 2) +
              Q(5) * T(1, 2) * T(1, 2);
  result(1) = 2 * Q(3) * T(0, 0) * T(0, 2) +
              Q(4) * (T(0, 0) * T(1, 2) + T(0, 2) * T(1, 0)) +
              2 * Q(5) * T(1, 0) * T(1, 2) + Q(1) * T(0, 0) + Q(2) * T(1, 0);
  result(2) = 2 * Q(3) * T(0, 1) * T(0, 2) +
              Q(4) * (T(0, 1) * T(1, 2) + T(0, 2) * T(1, 1)) +
              2 * Q(5) * T(1, 1) * T(1, 2) + Q(1) * T(0, 1) + Q(2) * T(1, 1);
  result(3) = Q(3) * T(0, 0) * T(0, 0) + Q(4) * T(0, 0) * T(1, 0) +
              Q(5) * T(1, 0) * T(1, 0);
  result(4) = 2 * Q(3) * T(0, 0) * T(0, 1) +
              Q(4) * (T(0, 0) * T(1, 1) + T(0, 1) * T(1, 0)) +
              2 * Q(5) * T(1, 0) * T(1, 1);
  result(5) = Q(3) * T(0, 1) * T(0, 1) + Q(4) * T(0, 1) * T(1, 1) +
              Q(5) * T(1, 1) * T(1, 1);

  return result;
}
Eigen::MatrixXd reparametrized_quadric_full(const Eigen::MatrixXd &Q,
                                            const Eigen::Matrix3d T) {
  Eigen::MatrixXd result(3, 6);
  for (auto i = 0; i < 3; ++i) {

    result(i, 0) = Q(i, 0) + Q(i, 1) * T(0, 2) + Q(i, 2) * T(1, 2) +
                   Q(i, 3) * T(0, 2) * T(0, 2) + Q(i, 4) * T(0, 2) * T(1, 2) +
                   Q(i, 5) * T(1, 2) * T(1, 2);
    result(i, 1) = 2 * Q(i, 3) * T(0, 0) * T(0, 2) +
                   Q(i, 4) * (T(0, 0) * T(1, 2) + T(0, 2) * T(1, 0)) +
                   2 * Q(i, 5) * T(1, 0) * T(1, 2) + Q(i, 1) * T(0, 0) +
                   Q(i, 2) * T(1, 0);
    result(i, 2) = 2 * Q(i, 3) * T(0, 1) * T(0, 2) +
                   Q(i, 4) * (T(0, 1) * T(1, 2) + T(0, 2) * T(1, 1)) +
                   2 * Q(i, 5) * T(1, 1) * T(1, 2) + Q(i, 1) * T(0, 1) +
                   Q(i, 2) * T(1, 1);
    result(i, 3) = Q(i, 3) * T(0, 0) * T(0, 0) + Q(i, 4) * T(0, 0) * T(1, 0) +
                   Q(i, 5) * T(1, 0) * T(1, 0);
    result(i, 4) = 2 * Q(i, 3) * T(0, 0) * T(0, 1) +
                   Q(i, 4) * (T(0, 0) * T(1, 1) + T(0, 1) * T(1, 0)) +
                   2 * Q(i, 5) * T(1, 0) * T(1, 1);
    result(i, 5) = Q(i, 3) * T(0, 1) * T(0, 1) + Q(i, 4) * T(0, 1) * T(1, 1) +
                   Q(i, 5) * T(1, 1) * T(1, 1);
  }

  return result;
}
Eigen::MatrixXd Jacobian_at_p(const Eigen::MatrixXd &Q, const vec2f &p_uv) {
  Eigen::MatrixXd J(3, 2);
  auto xu = evaluate_quadric_du(Q, p_uv);
  auto xv = evaluate_quadric_dv(Q, p_uv);
  J(0, 0) = xu.x;
  J(1, 0) = xu.y;
  J(2, 0) = xu.z;

  J(0, 1) = xv.x;
  J(1, 1) = xv.y;
  J(2, 1) = xv.z;

  return J;
}
Eigen::MatrixXd Jacobian_inverse(const Eigen::MatrixXd &J) {
  auto Jt = J.transpose();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> dec(Jt * J);
  return dec.inverse() * Jt;
}
std::pair<float, float> sin_and_cos_between_vectors(const vec3f &a,
                                                    const vec3f &b) {

  auto cos = dot(a, b);
  auto sin = length(b - cos * a);
  return {cos, sin};
}
vector<vec3f> compute_basis_tg_plane(const Eigen::MatrixXd &Q,
                                     const vec2f &p_uv) {
  Eigen::MatrixXd J = Jacobian_at_p(Q, p_uv);
  Eigen::Vector2d p;
  p << p_uv.x, p_uv.y;
  Eigen::Vector3d e0 = J * p;

  auto xu = evaluate_quadric_du(Q, p_uv);
  auto xv = evaluate_quadric_dv(Q, p_uv);
  auto n = normalize(cross(xu, xv));
  auto result = vector<vec3f>(2);
  result[0] = vec3f{(float)e0(0), (float)e0(1), (float)e0(2)};
  result[1] = rot_vect(result[0], n, pif / 2);
  return result;
}
vec2f basis_in_vert_tg_space(const Eigen::MatrixXd &Q, const vec3f &e,
                             const vec2f &p_uv) {
  auto J = Jacobian_at_p(Q, p_uv);
  auto Jinv = Jacobian_inverse(J);
  Eigen::Vector3d e0;
  e0 << e.x, e.y, e.z;
  Eigen::Vector2d result = Jinv * e0;
  return vec2f{(float)result(0), (float)result(1)};
}
Eigen::MatrixXd blended_quadric(const vector<Eigen::MatrixXd> &Q,
                                const vec2f &bary) {
  Eigen::MatrixXd QBlen(3, 6);
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 6; ++j) {
      QBlen(i, j) = (1 - bary.x - bary.y) * Q[0](i, j) + bary.x * Q[1](i, j) +
                    bary.y * Q[2](i, j);
    }
  }

  return QBlen;
}
Eigen::Matrix3d reparametrization_matrix_from_vert_to_point(
    const shape_data &data, const shape_geometry &geometry,
    const Eigen::MatrixXd &Q, const int vid, const mesh_point &p) {
  auto p_uv = coordinates_in_tangent_space_of_p(data, geometry, vid, p);
  if (length(p_uv) == 0)
    return Eigen::Matrix3d::Identity();

  auto bases = compute_basis_tg_plane(Q, p_uv);
  auto e0 = normalize(p_uv);
  auto e1 = normalize(basis_in_vert_tg_space(Q, bases[1], p_uv));
  return tranformation_matrix(e0, e1, p_uv);
}
std::tuple<Eigen::MatrixXd, vector<Eigen::Matrix3d>>
blended_quadric(const shape_data &data, const shape_geometry &geometry,
                const shape_op &op, const mesh_point &p) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto v0 = data.positions[vid0];
  auto v1 = data.positions[vid1];
  auto v2 = data.positions[vid2];
  auto pos_p = interpolate_triangle(v0, v1, v2, p.uv);
  auto verts = vector<int>{vid0, vid1, vid2};
  auto Q = vector<Eigen::MatrixXd>{op.quadrics[vid0], op.quadrics[vid1],
                                   op.quadrics[vid2]};
  auto T = vector<Eigen::Matrix3d>(3);
  for (auto i = 0; i < 3; ++i) {
    T[i] = reparametrization_matrix_from_vert_to_point(data, geometry, Q[i],
                                                       verts[i], p);
  }

  auto ep = get_basis_at_p(data.triangles, data.positions, p.face);

  auto theta0p = angle(pos_p - v0, ep[0]);
  if (dot(cross(ep[0], pos_p - v0), ep[2]) > 0)
    theta0p = 2 * pif - theta0p;
  auto theta1p = angle(pos_p - v1, ep[0]);
  if (dot(cross(ep[0], pos_p - v1), ep[2]) > 0)
    theta1p = 2 * pif - theta1p;
  auto theta2p = angle(pos_p - v2, ep[0]);
  if (dot(cross(ep[0], pos_p - v2), ep[2]) > 0)
    theta0p = 2 * pif - theta2p;

  Eigen::Matrix3d R0p;
  Eigen::Matrix3d R1p;
  Eigen::Matrix3d R2p;

  R0p << std::cos(theta0p), -std::sin(theta0p), 0, std::sin(theta0p),
      std::cos(theta0p), 0, 0, 0, 1;
  R1p << std::cos(theta1p), -std::sin(theta1p), 0, std::sin(theta1p),
      std::cos(theta1p), 0, 0, 0, 1;
  R2p << std::cos(theta2p), -std::sin(theta2p), 0, std::sin(theta2p),
      std::cos(theta2p), 0, 0, 0, 1;

  auto Tr = vector<Eigen::Matrix3d>{T[0] * R0p, T[1] * R1p, T[2] * R2p};
  Eigen::MatrixXd Qt0 = reparametrized_quadric_full(Q[0], Tr[0]);
  Eigen::MatrixXd Qt1 = reparametrized_quadric_full(Q[1], Tr[1]);
  Eigen::MatrixXd Qt2 = reparametrized_quadric_full(Q[2], Tr[2]);
  return {blended_quadric(vector<Eigen::MatrixXd>{Qt0, Qt1, Qt2}, p.uv), Tr};
}
Eigen::Matrix3d reparametrization(const vec2f &a, const vec2f &b,
                                  const vec2f &c, const vec2f &a1,
                                  const vec2f &b1, const vec2f &c1) {
  Eigen::MatrixXf A(6, 6);
  A << a.x, a.y, 1, 0, 0, 0, 0, 0, 0, a.x, a.y, 1, b.x, b.y, 1, 0, 0, 0, 0, 0,
      0, b.x, b.y, 1, c.x, c.y, 1, 0, 0, 0, 0, 0, 0, c.x, c.y, 1;
  Eigen::VectorXf B(6);
  B << a1.x, a1.y, b1.x, b1.y, c1.x, c1.y;
  Eigen::VectorXf sol = A.colPivHouseholderQr().solve(B);
  Eigen::Matrix3d T;
  T << sol(0), sol(1), sol(2), sol(3), sol(4), sol(5), 0, 0, 1;

  return T;
}

std::vector<Eigen::Matrix3d>
reparametrization_matrices(const shape_data &data,
                           const shape_geometry &geometry,
                           const mesh_point &p) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];

  auto ep = get_basis_at_p(data.triangles, data.positions, p.face);

  auto Tp = tranformation_matrix(
      ep[0], ep[1], ep[2], eval_position(data.triangles, data.positions, p));

  auto p0p = switch_reference_frame(Tp, data.positions[vid0]);
  auto p1p = switch_reference_frame(Tp, data.positions[vid1]);
  auto p2p = switch_reference_frame(Tp, data.positions[vid2]);

  auto tri0 = get_local_coordinates(data.triangles, data.positions,
                                    geometry.adjacencies, geometry.angles,
                                    geometry.v2t, vid0, vid1, vid2);

  auto T0 = reparametrization(vec2f{p0p.x, p0p.y}, vec2f{p1p.x, p1p.y},
                              vec2f{p2p.x, p2p.y}, zero2f, tri0[0], tri0[1]);
  auto tri1 = get_local_coordinates(data.triangles, data.positions,
                                    geometry.adjacencies, geometry.angles,
                                    geometry.v2t, vid1, vid0, vid2);

  auto T1 = reparametrization(vec2f{p0p.x, p0p.y}, vec2f{p1p.x, p1p.y},
                              vec2f{p2p.x, p2p.y}, tri1[0], zero2f, tri1[1]);
  auto tri2 = get_local_coordinates(data.triangles, data.positions,
                                    geometry.adjacencies, geometry.angles,
                                    geometry.v2t, vid2, vid0, vid1);

  auto T2 = reparametrization(vec2f{p0p.x, p0p.y}, vec2f{p1p.x, p1p.y},
                              vec2f{p2p.x, p2p.y}, tri2[0], tri2[1], zero2f);
  return {T0, T1, T2};
}
Eigen::VectorXd blended_quadric(const vector<Eigen::VectorXd> &Q,
                                const vec2f &bary) {
  Eigen::VectorXd QBlen(6);
  for (auto i = 0; i < 6; ++i) {

    QBlen(i) =
        (1 - bary.x - bary.y) * Q[0](i) + bary.x * Q[1](i) + bary.y * Q[2](i);
  }

  return QBlen;
}

vec3f normal_at_p_in_tangent_space(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const int vid,
                                   const mesh_point &p) {
  auto uv = coordinates_in_tangent_space_of_p(data, geometry, vid, p);
  auto xu = evaluate_quadric_du(op, vid, uv);
  auto xv = evaluate_quadric_dv(op, vid, uv);
  return normalize(cross(xu, xv));
}

vector<Eigen::MatrixXd>
quadric_at_p_on_edge(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const mesh_point &p,
                     const mesh_point &q, const int vid) {
  auto k0 = find(data.triangles[p.face], vid);
  auto k1 = find(data.triangles[q.face], vid);
  auto T0 = reparametrization_matrices(data, geometry, p)[k0];
  auto T1 = reparametrization_matrices(data, geometry, q)[k1];

  Eigen::Matrix3d T = T1.inverse() * T0;

  auto Q = op.quadrics[vid];

  auto Q0 = reparametrized_quadric_full(Q, T0);
  auto Q1 = reparametrized_quadric_full(Q, T1);
  Q1 = reparametrized_quadric_full(Q1, T);

  return {Q0, Q1};
}
Eigen::MatrixXd blend_quadric_on_edge(const Eigen::MatrixXd &Q0,
                                      const Eigen::MatrixXd &Q1,
                                      const float &t) {
  Eigen::MatrixXd QBlended(3, 6);
  for (auto i = 0; i < 3; ++i) {
    for (auto j = 0; j < 6; ++j) {
      QBlended(i, j) = (1 - t) * Q0(i, j) + t * Q1(i, j);
    }
  }

  return QBlended;
}

vector<Eigen::MatrixXd> quadric_at_p_on_edge_w_repara(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const mesh_point &p, const mesh_point &q) {
  auto k = find(geometry.adjacencies[p.face], q.face);

  assert(k != -1);

  auto Tp = reparametrization_matrices(data, geometry, p);
  auto T0p = Tp[k];
  auto T1p = Tp[(k + 1) % 3];
  auto vid0 = data.triangles[p.face][k];
  auto vid1 = data.triangles[p.face][(k + 1) % 3];
  auto Q0 = op.quadrics[vid0];
  auto Q1 = op.quadrics[vid1];
  auto Q0p = reparametrized_quadric_full(Q0, T0p);
  auto Q1p = reparametrized_quadric_full(Q1, T1p);
  auto bary = get_bary(p.uv);
  auto Qp = blend_quadric_on_edge(Q0p, Q1p, 1 - bary[k]);
  auto h = find(data.triangles[q.face], vid0);
  auto Tq = reparametrization_matrices(data, geometry, q);
  auto T0q = Tq[h];
  auto T1q = Tq[(h + 1) % 3];
  Eigen::Matrix3d TC0 = T0q.inverse() * T0p;
  Eigen::Matrix3d TC1 = T1q.inverse() * T1p;

  auto Q0q = reparametrized_quadric_full(Q0, T0q);
  auto Q1q = reparametrized_quadric_full(Q1, T1q);
  Q0q = reparametrized_quadric_full(Q0q, TC0);
  Q1q = reparametrized_quadric_full(Q1q, TC1);
  auto Rp = blend_quadric_on_edge(Q0q, Q1q, 1 - bary[k]);

  // auto Q0 = quadric_at_p_on_edge(data, geometry, op, p, q, vid0);
  // auto Q1 = quadric_at_p_on_edge(data, geometry, op, p, q, vid1);

  return {Qp, Rp};
}
std::tuple<Eigen::Matrix2f, vector<vec3f>>
metric_tensor(const shape_data &data, const shape_geometry &geometry,
              const shape_op &op, const int &vid) {
  auto c = op.quadrics[vid];

  auto xu = vec3f{(float)c(0, 1), (float)c(1, 1), (float)c(2, 1)};
  auto xv = vec3f{(float)c(0, 2), (float)c(1, 2), (float)c(2, 2)};
  auto xuu = vec3f{2 * (float)c(0, 3), 2 * (float)c(1, 3), 2 * (float)c(2, 3)};
  auto xuv = vec3f{(float)c(0, 4), (float)c(1, 4), (float)c(2, 4)};
  auto xvv = vec3f{2 * (float)c(0, 5), 2 * (float)c(1, 5), 2 * (float)c(2, 5)};

  Eigen::Matrix2f g;
  g << dot(xu, xu), dot(xu, xv), dot(xv, xu), dot(xv, xv);

  return {g, {xu, xv, xuu, xvv, xuv}};
}
vec3f quadric_at_p(const shape_data &data, const shape_geometry &geometry,
                   const shape_op &op, const mesh_point &p) {
  auto [Qblen, T] = blended_quadric(data, geometry, op, p);
  return evaluate_quadric(Qblen, zero2f);
}
std::tuple<Eigen::Matrix2f, vector<vec3f>, vector<Eigen::Matrix3d>>
metric_tensor(const shape_data &data, const shape_geometry &geometry,
              const shape_op &op, const mesh_point &p) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];

  auto Q0 = op.quadrics[vid0];
  auto Q1 = op.quadrics[vid1];
  auto Q2 = op.quadrics[vid2];

  vector<Eigen::MatrixXd> Q = {Q0, Q1, Q2};
  vector<Eigen::MatrixXd> QT(3);
  // for (auto i = 0; i < 3; ++i) {
  //   auto [pos_in_tri, pos_tg] = get_coordinates_for_quadric_fitting(
  //       data, geometry, data.triangles[tid][i], tid);
  //   QT[i] = quadric_fitting(Q[i], pos_in_tri, pos_tg);
  // }

  auto [Qblen, T] = blended_quadric(data, geometry, op, p);
  // auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  // auto Tp = tranformation_matrix(ep[0], ep[1], ep[2], data.positions[vid0]);
  // auto pos_p_3d = switch_reference_frame(
  //     Tp, eval_position(data.triangles, data.positions, p));
  // auto pos_p = vec2f{pos_p_3d.x, pos_p_3d.y};
  auto xu = evaluate_quadric_du(Qblen, zero2f);
  auto xv = evaluate_quadric_dv(Qblen, zero2f);
  auto xuu = evaluate_quadric_duu(Qblen, zero2f);
  auto xvv = evaluate_quadric_dvv(Qblen, zero2f);
  auto xuv = evaluate_quadric_duv(Qblen, zero2f);

  Eigen::Matrix2f g;
  g << dot(xu, xu), dot(xu, xv), dot(xv, xu), dot(xv, xv);

  return {g, {xu, xv, xuu, xvv, xuv}, T};
}
vec3f normal_at_vid(const shape_data &data, const shape_geometry &geometry,
                    const shape_op &op, const int vid) {

  auto [g, X] = metric_tensor(data, geometry, op, vid);

  return normalize(cross(X[0], X[1]));
}
vec3f normal_at_p(const shape_data &data, const shape_geometry &geometry,
                  const shape_op &op, const mesh_point &p) {

  auto [is_vert, kv] = point_is_vert(p);
  auto tid = p.face;
  if (is_vert)
    return normal_at_vid(data, geometry, op, data.triangles[tid][kv]);

  auto [Qblen, T] = blended_quadric(data, geometry, op, p);

  auto xu = evaluate_quadric_du(Qblen, zero2f);
  auto xv = evaluate_quadric_dv(Qblen, zero2f);

  return normalize(cross(xu, xv));
}
vec3f normal_at_p_on_edge(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const mesh_point &p, const vec2i &edge) {

  auto tid = p.face;
  // auto vid0 = edge.x;
  // auto vid1 = edge.y;
  // auto bary = get_bary(p.uv);
  // auto k = find(data.triangles[tid], vid0);
  // auto h = find(data.triangles[tid], vid1);
  // assert(k != -1);
  // assert(h != -1);
  // auto Q0 = op.quadrics[vid0];
  // auto Q1 = op.quadrics[vid1];
  // vector<Eigen::MatrixXd> Q_vid = {Q0, Q1};
  // vector<Eigen::MatrixXd> Q(2);

  // for (auto i = 0; i < 2; ++i) {
  //   auto [pos_in_tri, pos_tg] =
  //       get_coordinates_for_quadric_fitting(data, geometry, edge[i], tid);
  //   Q[i] = quadric_fitting(Q_vid[i], pos_in_tri, pos_tg);
  // }

  auto [Qblen, T] = blended_quadric(
      data, geometry, op, p); // blend_quadric_on_edge(Q[0], Q[1], 1 - bary[k]);

  // auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  // auto Tp = tranformation_matrix(ep[0], ep[1], ep[2],
  //                                data.positions[data.triangles[tid].x]);
  // auto pos_p = switch_reference_frame(
  //     Tp, eval_position(data.triangles, data.positions, p));
  auto xu = evaluate_quadric_du(Qblen, zero2f);
  auto xv = evaluate_quadric_dv(Qblen, zero2f);

  return normalize(cross(xu, xv));
}
float curvature_at_p_on_edge(const shape_data &data,
                             const shape_geometry &geometry, const shape_op &op,
                             const mesh_point &p, const vec2i &edge) {

  auto tid = p.face;
  // auto vid0 = edge.x;
  // auto vid1 = edge.y;
  // auto bary = get_bary(p.uv);
  // auto k = find(data.triangles[tid], vid0);
  // auto h = find(data.triangles[tid], vid1);
  // assert(k != -1);
  // assert(h != -1);
  // auto Q0 = op.quadrics[vid0];
  // auto Q1 = op.quadrics[vid1];
  // string Q0_name = "quadric_at_" + std::to_string(vid0) + "_" +
  //                  std::to_string(tid) + "_new_method";
  // string Q1_name = "quadric_at_" + std::to_string(vid1) + "_" +
  //                  std::to_string(tid) + "_new_method";

  // write_quadric(Q0, Q0_name);
  // write_quadric(Q1, Q1_name);
  // vector<Eigen::MatrixXd> Q_vid = {Q0, Q1};
  // vector<Eigen::MatrixXd> Q(2);

  // for (auto i = 0; i < 2; ++i) {
  //   auto [pos_in_tri, pos_tg] =
  //       get_coordinates_for_quadric_fitting(data, geometry, edge[i], tid);
  //   Q[i] = quadric_fitting(Q_vid[i], pos_in_tri, pos_tg);
  // }

  auto [Qblen, T] = blended_quadric(data, geometry, op, p);
  // blend_quadric_on_edge(Q[0], Q[1], 1 - bary[k]);
  // string blended_name =
  //     "Blended_quadric_in_" + std::to_string(tid) + "_new_method";
  // write_quadric(Qblen, blended_name);
  // auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  // auto Tp = tranformation_matrix(ep[0], ep[1], ep[2],
  //                                data.positions[data.triangles[tid].x]);
  // auto pos_p = switch_reference_frame(
  //     Tp, eval_position(data.triangles, data.positions, p));
  auto xu = evaluate_quadric_du(Qblen, zero2f);
  auto xv = evaluate_quadric_dv(Qblen, zero2f);
  auto xuu = evaluate_quadric_duu(Qblen, zero2f);
  auto xvv = evaluate_quadric_dvv(Qblen, zero2f);
  auto xuv = evaluate_quadric_duv(Qblen, zero2f);
  auto n = normalize(cross(xu, xv));
  auto E = dot(xu, xu);
  auto F = dot(xu, xv);
  auto G = dot(xv, xv);
  auto L = dot(xuu, n);
  auto M = dot(xuv, n);
  auto N = dot(xvv, n);
  return (L * N - M * M) / (E * G - F * F);
}

float gaussian_curvature_at_vid(const shape_data &data,
                                const shape_geometry &geometry,
                                const shape_op &op, const int &vid) {
  auto [g, Xuv] = metric_tensor(data, geometry, op, vid);
  auto xu = Xuv[0];
  auto xv = Xuv[1];
  auto xuu = Xuv[2];
  auto xvv = Xuv[3];
  auto xuv = Xuv[4];
  auto n = normalize(cross(xu, xv));
  auto L = dot(xuu, n);
  auto M = dot(xuv, n);
  auto N = dot(xvv, n);
  auto num = L * N - pow(M, 2);
  auto den = g.determinant();

  return (float)num / den;
}
float gaussian_curvature_at_p(const shape_data &data,
                              const shape_geometry &geometry,
                              const shape_op &op, const mesh_point &p) {
  auto [g, Xuv, T] = metric_tensor(data, geometry, op, p);
  auto xu = Xuv[0];
  auto xv = Xuv[1];
  auto xuu = Xuv[2];
  auto xvv = Xuv[3];
  auto xuv = Xuv[4];
  auto n = normalize(cross(xu, xv));
  auto L = dot(xuu, n);
  auto M = dot(xuv, n);
  auto N = dot(xvv, n);
  auto num = L * N - pow(M, 2);
  auto den = g.determinant();
  return (float)num / den;
}
float gaussian_curvature_at_p_lerp(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const mesh_point &p) {
  auto tid = p.face;

  auto c0 =
      gaussian_curvature_at_vid(data, geometry, op, data.triangles[tid][0]);
  auto c1 =
      gaussian_curvature_at_vid(data, geometry, op, data.triangles[tid][1]);
  auto c2 =
      gaussian_curvature_at_vid(data, geometry, op, data.triangles[tid][2]);

  return interpolate_triangle(c0, c1, c2, p.uv);
}

Eigen::VectorXd sample_field(const shape_data &data,
                             const shape_geometry &geometry,
                             const vector<float> &f, const int vid) {
  auto [nbr, thetas, lens] = uniform_stencil(data, geometry, vid, 36);
  Eigen::VectorXd result(37);
  result(0) = f[vid];
  for (auto i = 0; i < 36; ++i) {
    auto bary = get_bary(nbr[i].uv);
    result(i + 1) = bary.x * f[data.triangles[nbr[i].face][0]] +
                    bary.y * f[data.triangles[nbr[i].face][1]] +
                    bary.z * f[data.triangles[nbr[i].face][2]];
  }
  return result;
}

vec3f gradient_at_vid(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f,
                      const int &vid) {
  auto c = op.quadrics[vid];
  auto xu = vec3f{(float)c(0, 1), (float)c(1, 1), (float)c(2, 1)};
  auto xv = vec3f{(float)c(0, 2), (float)c(1, 2), (float)c(2, 2)};
  Eigen::Matrix2f g;
  g << dot(xu, xu), dot(xu, xv), dot(xu, xv), dot(xv, xv);

  auto field_at_neighbors = sample_field(data, geometry, f, vid);
  // Eigen::VectorXd Q = op.CMat[vid] * field_at_neighbors;
  Eigen::VectorXd Q = op.CMat[vid].solve(op.Rhs[vid] * field_at_neighbors);
  Eigen::Vector2f diff_f;
  diff_f << Q(1), Q(2);
  Eigen::Vector2f components = g.inverse() * diff_f;
  auto grad = xu * components(0) + xv * components(1);
  return grad / (2 * std::sqrt((float)Q(0)));
}
std::tuple<vec3f, float>
gradient_at_vid_w_magnitude(const shape_data &data,
                            const shape_geometry &geometry, const shape_op &op,
                            const vector<float> &f, const int &vid) {
  auto c = op.quadrics[vid];
  auto xu = vec3f{(float)c(0, 1), (float)c(1, 1), (float)c(2, 1)};
  auto xv = vec3f{(float)c(0, 2), (float)c(1, 2), (float)c(2, 2)};
  Eigen::Matrix2f g;
  g << dot(xu, xu), dot(xu, xv), dot(xu, xv), dot(xv, xv);

  auto field_at_neighbors = sample_field(data, geometry, f, vid);
  // Eigen::VectorXd Q = op.CMat[vid] * field_at_neighbors;
  Eigen::VectorXd Q = op.CMat[vid].solve(op.Rhs[vid] * field_at_neighbors);
  Eigen::Vector2f diff_f;
  diff_f << Q(1), Q(2);
  Eigen::Vector2f components = g.inverse() * diff_f;
  auto grad = xu * components(0) + xv * components(1);
  auto mag = g(0, 0) * components(0) * components(0) +
             2 * g(0, 1) * components(0) * components(1) +
             g(1, 1) * components(1) * components(1);
  return {grad, mag};
  // auto e0 = polar_basis(data.triangles, data.positions, geometry.v2t,
  //                       data.normals, vid);
  // auto e1 = rot_vect(e0, data.normals[vid], pif / 2);
  // auto grad = e0 * components(0) + e1 * components(1);
  // return grad / (2 * std::sqrt((float)Q(0)));
}
float field_at_p(const shape_data &data, const shape_geometry &geometry,
                 const shape_op &op, const vector<float> &f,
                 const mesh_point &p) {
  auto [g, xuv, T] = metric_tensor(data, geometry, op, p);
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto verts = vector<int>{vid0, vid1, vid2};
  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);

  // Eigen::VectorXd Q0 = op.CMat[vid0] * f0;
  // Eigen::VectorXd Q1 = op.CMat[vid1] * f1;
  // Eigen::VectorXd Q2 = op.CMat[vid2] * f2;
  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f fpi;
  for (auto i = 0; i < 3; ++i) {
    auto p_uv = coordinates_in_tangent_space_of_p(data, geometry, verts[i], p);
    fpi[i] = evaluate_quadric(Q[i], p_uv);
  }
  return (1 - p.uv.x - p.uv.y) * fpi[0] + p.uv.x * fpi[1] + p.uv.y * fpi[2];
}
float field_at_p_lerp(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f,
                      const mesh_point &p) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto verts = vector<int>{vid0, vid1, vid2};
  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);

  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f fpi;
  for (auto i = 0; i < 3; ++i) {
    fpi[i] = evaluate_quadric(Q[i], zero2f);
  }
  return (1 - p.uv.x - p.uv.y) * fpi[0] + p.uv.x * fpi[1] + p.uv.y * fpi[2];
}

vec3f gradient_at_p(const shape_data &data, const shape_geometry &geometry,
                    const shape_op &op, const vector<float> &f,
                    const mesh_point &p, const bool &squared = true) {

  auto [g, xuv, T] = metric_tensor(data, geometry, op, p);
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto verts = vector<int>{vid0, vid1, vid2};
  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);

  // Eigen::VectorXd Q0 = op.CMat[vid0] * f0;
  // Eigen::VectorXd Q1 = op.CMat[vid1] * f1;
  // Eigen::VectorXd Q2 = op.CMat[vid2] * f2;
  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f ku;
  vec3f kv;
  vec3f fpi;
  for (auto i = 0; i < 3; ++i) {
    auto Qt = reparametrized_quadric_full(Q[i], T[i]);
    ku[i] = evaluate_quadric_du(Qt, zero2f);
    kv[i] = evaluate_quadric_dv(Qt, zero2f);
    fpi[i] = evaluate_quadric(Qt, zero2f);
  }

  auto ep = get_basis_at_p(data.triangles, data.positions, tid);
  auto Tp = tranformation_matrix(
      ep[0], ep[1], ep[2], eval_position(data.triangles, data.positions, p));
  // auto pos_p_3d = switch_reference_frame(
  //     Tp, eval_position(data.triangles, data.positions, p));
  // auto pos_p = vec2f{pos_p_3d.x, pos_p_3d.y};
  auto fu = (1 - p.uv.x - p.uv.y) * ku[0] + p.uv.x * ku[1] + p.uv.y * ku[2];
  auto fv = (1 - p.uv.x - p.uv.y) * kv[0] + p.uv.x * kv[1] + p.uv.y * kv[2];
  auto fp = (1 - p.uv.x - p.uv.y) * fpi[0] + p.uv.x * fpi[1] + p.uv.y * fpi[2];
  Eigen::Vector2f diff_f;
  diff_f << fu, fv;
  Eigen::Vector2f components = g.inverse() * diff_f;

  auto grad = xuv[0] * components(0) + xuv[1] * components(1);

  // grad = switch_reference_frame_vector(Tp, grad);
  // grad[2] = 0;
  // grad = switch_reference_frame_vector(Tp.inverse(), grad);
  if (squared)
    grad /= (2 * std::sqrt(fp));

  return grad;
}
vec3f gradient_at_p_lerp(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const vector<float> &f,
                         const mesh_point &p, const bool &squared = true) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto g0 = gradient_at_vid(data, geometry, op, f, vid0);
  auto g1 = gradient_at_vid(data, geometry, op, f, vid1);
  auto g2 = gradient_at_vid(data, geometry, op, f, vid2);
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, g0,
                  data.normals, vid0, p.face, V2T);
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, g1,
                  data.normals, vid1, p.face, V2T);
  parallel_transp(geometry.angles, geometry.total_angles, data.triangles,
                  data.positions, geometry.adjacencies, geometry.v2t, g2,
                  data.normals, vid2, p.face, V2T);

  auto grad = (1 - p.uv.x - p.uv.y) * g0 + p.uv.x * g1 + p.uv.y * g2;
  auto fp = field_at_p_lerp(data, geometry, op, f, p);
  if (squared)
    grad /= (2 * std::sqrt(fp));

  return grad;
}
std::tuple<vec3f, float>
gradient_at_p_w_magnitude(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f, const mesh_point &p) {
  auto [g, xuv, T] = metric_tensor(data, geometry, op, p);
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];

  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);

  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f ku;
  vec3f kv;
  vec3f fpi;
  for (auto i = 0; i < 3; ++i) {
    auto Qt = reparametrized_quadric_full(Q[i], T[i]);
    ku[i] = evaluate_quadric_du(Qt, zero2f);
    kv[i] = evaluate_quadric_dv(Qt, zero2f);
    fpi[i] = evaluate_quadric(Qt, zero2f);
  }

  auto fu = (1 - p.uv.x - p.uv.y) * ku[0] + p.uv.x * ku[1] + p.uv.y * ku[2];
  auto fv = (1 - p.uv.x - p.uv.y) * kv[0] + p.uv.x * kv[1] + p.uv.y * kv[2];
  auto fp = (1 - p.uv.x - p.uv.y) * fpi[0] + p.uv.x * fpi[1] + p.uv.y * fpi[2];
  Eigen::Vector2f diff_f;
  diff_f << fu, fv;
  Eigen::Vector2f components = g.inverse() * diff_f;

  auto grad = xuv[0] * components(0) + xuv[1] * components(1);
  auto mag = g(0, 0) * components(0) * components(0) +
             2 * g(0, 1) * components(0) * components(1) +
             g(1, 1) * components(1) * components(1);
  return {grad, mag};
}

float laplacian_at_vid(const shape_data &data, const shape_geometry &geometry,
                       const shape_op &op, const vector<float> &f,
                       const int &vid) {
  auto c = op.quadrics[vid];
  auto xu = vec3d{c(0, 1), c(1, 1), c(2, 1)};
  auto xv = vec3d{c(0, 2), c(1, 2), c(2, 2)};
  auto xuu = vec3d{2 * c(0, 3), 2 * c(1, 3), 2 * c(2, 3)};
  auto xuv = vec3d{c(0, 4), c(1, 4), c(2, 4)};
  auto xvv = vec3d{2 * c(0, 5), 2 * c(1, 5), 2 * c(2, 5)};
  Eigen::Matrix2d g;
  g << dot(xu, xu), dot(xu, xv), dot(xu, xv), dot(xv, xv);
  auto det = g.determinant();

  auto field_at_neighbors = sample_field(data, geometry, f, vid);
  // Eigen::VectorXd Q = op.CMat[vid] * field_at_neighbors;
  Eigen::VectorXd Q = op.CMat[vid].solve(op.Rhs[vid] * field_at_neighbors);
  auto fu = Q(1);
  auto fv = Q(2);
  auto fuu = 2 * Q(3);
  auto fvv = 2 * Q(5);
  auto fuv = Q(4);
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

  return g_deltau * fu + g_deltav * fv + g_deltauu * fuu + g_deltavv * fvv +
         g_deltauv * fuv;
}
float laplacian_at_p(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f,
                     const mesh_point &p) {
  auto [g, Xuv, T] = metric_tensor(data, geometry, op, p);
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];

  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);
  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f fui;
  vec3f fvi;
  vec3f fuui;
  vec3f fvvi;
  vec3f fuvi;
  // auto T = reparametrization_matrices(data, geometry, p);
  for (auto i = 0; i < 3; ++i) {
    auto Qt = reparametrized_quadric_full(Q[i], T[i]);
    fui[i] = evaluate_quadric_du(Qt, zero2f);
    fvi[i] = evaluate_quadric_dv(Qt, zero2f);
    fuui[i] = evaluate_quadric_duu(Qt, zero2f);
    fuvi[i] = evaluate_quadric_duv(Qt, zero2f);
    fvvi[i] = evaluate_quadric_dvv(Qt, zero2f);
  }
  auto fu = (1 - p.uv.x - p.uv.y) * fui[0] + p.uv.x * fui[1] + p.uv.y * fui[2];
  auto fv = (1 - p.uv.x - p.uv.y) * fvi[0] + p.uv.x * fvi[1] + p.uv.y * fvi[2];
  auto fuu =
      (1 - p.uv.x - p.uv.y) * fuui[0] + p.uv.x * fuui[1] + p.uv.y * fuui[2];
  auto fvv =
      (1 - p.uv.x - p.uv.y) * fvvi[0] + p.uv.x * fvvi[1] + p.uv.y * fvvi[2];
  auto fuv =
      (1 - p.uv.x - p.uv.y) * fuvi[0] + p.uv.x * fuvi[1] + p.uv.y * fuvi[2];
  auto det = g.determinant();

  auto xu = Xuv[0];
  auto xv = Xuv[1];
  auto xuu = Xuv[2];
  auto xvv = Xuv[3];
  auto xuv = Xuv[4];
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

  return g_deltau * fu + g_deltav * fv + g_deltauu * fuu + g_deltauv * fuv +
         g_deltavv * fvv;
}
Eigen::Matrix2f Hessian_at_p(const shape_data &data,
                             const shape_geometry &geometry, const shape_op &op,
                             const vector<float> &f, const mesh_point &p) {
  auto [g, Xuv, T] = metric_tensor(data, geometry, op, p);
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];

  auto f0 = sample_field(data, geometry, f, vid0);
  auto f1 = sample_field(data, geometry, f, vid1);
  auto f2 = sample_field(data, geometry, f, vid2);
  Eigen::VectorXd Q0 = op.CMat[vid0].solve(op.Rhs[vid0] * f0);
  Eigen::VectorXd Q1 = op.CMat[vid1].solve(op.Rhs[vid1] * f1);
  Eigen::VectorXd Q2 = op.CMat[vid2].solve(op.Rhs[vid2] * f2);
  vector<Eigen::VectorXd> Q = {Q0, Q1, Q2};
  vec3f fui;
  vec3f fvi;
  vec3f fuui;
  vec3f fvvi;
  vec3f fuvi;
  // auto T = reparametrization_matrices(data, geometry, p);
  for (auto i = 0; i < 3; ++i) {
    auto Qt = reparametrized_quadric_full(Q[i], T[i]);
    fui[i] = evaluate_quadric_du(Qt, zero2f);
    fvi[i] = evaluate_quadric_dv(Qt, zero2f);
    fuui[i] = evaluate_quadric_duu(Qt, zero2f);
    fuvi[i] = evaluate_quadric_duv(Qt, zero2f);
    fvvi[i] = evaluate_quadric_dvv(Qt, zero2f);
  }
  auto fu = (1 - p.uv.x - p.uv.y) * fui[0] + p.uv.x * fui[1] + p.uv.y * fui[2];
  auto fv = (1 - p.uv.x - p.uv.y) * fvi[0] + p.uv.x * fvi[1] + p.uv.y * fvi[2];
  auto fuu =
      (1 - p.uv.x - p.uv.y) * fuui[0] + p.uv.x * fuui[1] + p.uv.y * fuui[2];
  auto fvv =
      (1 - p.uv.x - p.uv.y) * fvvi[0] + p.uv.x * fvvi[1] + p.uv.y * fvvi[2];
  auto fuv =
      (1 - p.uv.x - p.uv.y) * fuvi[0] + p.uv.x * fuvi[1] + p.uv.y * fuvi[2];

  auto inv_g = g.inverse();
  auto xu = Xuv[0];
  auto xv = Xuv[1];
  auto xuu = Xuv[2];
  auto xvv = Xuv[3];
  auto xuv = Xuv[4];

  Eigen::Matrix2f gjm0;
  Eigen::Matrix2f gjm1;
  Eigen::Vector2f f_I;
  Eigen::Vector2f f_II0;
  Eigen::Vector2f f_II1;
  f_I << fu, fv;
  f_II0 << fuu, fuv;
  f_II1 << fuv, fvv;
  gjm0 << dot(xu, xuu), dot(xv, xuu), dot(xu, xuv), dot(xv, xuv);
  gjm1 << dot(xu, xuv), dot(xv, xuv), dot(xu, xvv), dot(xv, xvv);

  // auto A = -inv_g.row(0) * gjm0 * inv_g * f_I + inv_g.row(0) * f_II;
  // auto C = -inv_g.row(0) * gjm1 * inv_g * f_I + inv_g.row(0) * f_II;
  // f_II << fuv, fvv;
  // auto B = -inv_g.row(1) * gjm0 * inv_g * f_I + inv_g.row(1) * f_II;
  // auto D = -inv_g.row(1) * gjm1 * inv_g * f_I + inv_g.row(1) * f_II;
  Eigen::Matrix2f f_II;
  Eigen::Matrix2f g_ijm;
  Eigen::Vector2f g_f_I = inv_g * f_I;
  Eigen::Matrix2f II;
  Eigen::Vector2f gj00;
  Eigen::Vector2f gj01;
  Eigen::Vector2f gj11;

  gj00 << dot(xu, xuu), dot(xv, xuu);
  gj01 << dot(xu, xuv), dot(xv, xuv);
  gj11 << dot(xu, xvv), dot(xv, xvv);
  f_II << fuu, fuv, fuv, fvv;
  II << gj00.dot(g_f_I), gj01.dot(g_f_I), gj01.dot(g_f_I), gj11.dot(g_f_I);

  // auto A = -inv_g.row(0) * gjm0 * inv_g * f_I + inv_g.row(0) * f_II0;
  // auto B = -inv_g.row(1) * gjm0 * inv_g * f_I + inv_g.row(1) * f_II0;
  // auto C = -inv_g.row(0) * gjm1 * inv_g * f_I + inv_g.row(0) * f_II1;
  // auto D = -inv_g.row(1) * gjm1 * inv_g * f_I + inv_g.row(1) * f_II1;

  Eigen::Matrix2f H = inv_g * (f_II - II);
  // H << A, B, C, D;

  return H;
}
float laplacian_at_p_lerp(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f, const mesh_point &p) {
  auto tid = p.face;
  auto vid0 = data.triangles[tid][0];
  auto vid1 = data.triangles[tid][1];
  auto vid2 = data.triangles[tid][2];
  auto verts = vector<int>{vid0, vid1, vid2};
  auto l0 = laplacian_at_vid(data, geometry, op, f, vid0);
  auto l1 = laplacian_at_vid(data, geometry, op, f, vid1);
  auto l2 = laplacian_at_vid(data, geometry, op, f, vid2);

  return (1 - p.uv.x - p.uv.y) * l0 + p.uv.x * l1 + p.uv.y * l2;
}

vector<vec3f> metric_tensor_eigen_decomposition(const shape_data &data,
                                                const shape_geometry &geometry,
                                                const shape_op &op,
                                                const mesh_point &p) {

  // auto [is_vert, kv] = point_is_vert(p);
  // Eigen::Matrix2f g;
  // vector<vec3f> xuv;
  // vector<Eigen::Matrix3d> Ti;
  // if (is_vert) {
  //   auto vid = data.triangles[p.face][kv];
  //   std::tie(g, xuv) = metric_tensor(data, geometry, op, vid);
  // } else {
  //   std::tie(g, xuv, Ti) = metric_tensor(data, geometry, op, p);
  // }

  // Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> dec(g);
  // if (dec.info() != Eigen::Success)
  //   std::cerr << "Error in computing the eigenvalues of the metric tensor";

  // auto lambda = dec.eigenvalues();
  // auto U = dec.eigenvectors();
  // if (lambda(0) < 0 || lambda(1) < 0)
  //   std::cerr << "Error in computing the eigenvalues of the metric "
  //                "tensor:negatives eigenvalues";
  // auto a = 1 / std::sqrt(lambda(0));
  // auto b = 1 / std::sqrt(lambda(1));

  // auto v0 = xuv[0] * U(0, 0) + xuv[1] * U(1, 0);
  // auto v1 = xuv[0] * U(0, 1) + xuv[1] * U(1, 1);

  // return {v0 / a, v1 / b};
}
vector<mesh_point> internal_points(const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto points = vector<mesh_point>(N);
  for (auto i = 0; i <= n; ++i) {
    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      points[entry] = mesh_point{tid, {(float)j / n, (float)i / n}};
    }
  }

  return points;
}
float interpolate_field_in_tid(const vector<vec3i> &triangles,
                               const vector<float> &f, const mesh_point &p) {
  return (1 - p.uv.x - p.uv.y) * f[triangles[p.face].x] +
         p.uv.x * f[triangles[p.face].y] + p.uv.y * f[triangles[p.face].z];
}
vector<vec3f> PCE_at_PS(const shape_data &data, const vector<float> &f,
                        const vector<mesh_point> &PS) {
  auto result = vector<vec3f>(PS.size());
  for (auto i = 0; i < PS.size(); ++i) {
    result[i] =
        compute_PCE(data.triangles, data.positions, f, PS[i]) /
        (2 * std::sqrt(interpolate_field_in_tid(data.triangles, f, PS[i])));
  }

  return result;
}

std::tuple<vector<vec3f>, vector<vec3f>>
pos_and_normals_PS(const shape_data &data, const vector<mesh_point> &PS) {
  auto pos = vector<vec3f>(PS.size());
  auto normals = vector<vec3f>(PS.size());

  for (auto i = 0; i < pos.size(); ++i) {
    pos[i] = eval_position(data.triangles, data.positions, PS[i]);
    normals[i] = tid_normal(data.triangles, data.positions, PS[i].face);
  }

  return {pos, normals};
}
void export_field(const vector<float> &field, const string &filename) {
  string name = "/Users/claudiomancinelli/Documents/GitHub/"
                "InterpolatingQuadric/build/";
  name.append(filename);
  std::ofstream outfile;
  outfile.open(name);
  assert(outfile.is_open());
  outfile << "SCALARFIELD"
          << " " << field.size() << "\n";
  for (auto i = 0; i < field.size(); ++i) {
    outfile << field[i] << "\n";
  }
  outfile.close();
}
void export_vector_field(const vector<vec3f> &grad, const string &name) {
  // string name = "/Users/claudiomancinelli/Documents/GitHub/"
  //               "InterpolatingQuadric/build/";
  // name.append(filename);
  std::ofstream outfile;
  outfile.open(name);
  assert(outfile.is_open());
  outfile << "VECTOR_FIELD"
          << " " << grad.size() << "\n";
  for (auto i = 0; i < grad.size(); ++i) {
    outfile << grad[i].x << " " << grad[i].y << " " << grad[i].z << "\n";
  }
  outfile.close();
}
std::tuple<vector<vec3f>, vector<vec3i>>
subdivide_tid(const shape_data &data, const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);

  auto pos = vector<vec3f>(N);
  auto faces = vector<vec3i>((int)pow(n, 2));

  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};

      pos[entry] = eval_position(data.triangles, data.positions, curr_point);
      if (j == 0 || i == n)
        continue;

      auto face_entry = 2 * entry - 2 * (i + 1) - i;
      faces[face_entry] = {entry, entry + n - i, entry - 1};
      if (j != n - i)
        faces[face_entry + 1] = {entry, entry + n - i + 1, entry + n - i};
    }
  }

  return {pos, faces};
}
std::tuple<vector<vec3f>, vector<vec3i>>
bumpy_triangle(const shape_data &data, const shape_geometry &geometry,
               const shape_op &op, const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);

  auto pos = vector<vec3f>(N);
  auto faces = vector<vec3i>((int)pow(n, 2));

  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};

      pos[entry] = quadric_at_p(data, geometry, op, curr_point);
      if (j == 0 || i == n)
        continue;

      auto face_entry = 2 * entry - 2 * (i + 1) - i;
      faces[face_entry] = {entry, entry + n - i, entry - 1};
      if (j != n - i)
        faces[face_entry + 1] = {entry, entry + n - i + 1, entry + n - i};
    }
  }

  return {pos, faces};
}
std::tuple<vector<vec3f>, vector<vec3i>, vector<vec3f>>
bumpy_triangle_w_normals(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);

  auto pos = vector<vec3f>(N);
  auto normals = vector<vec3f>(N);
  auto faces = vector<vec3i>((int)pow(n, 2));

  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};

      pos[entry] = quadric_at_p(data, geometry, op, curr_point);
      normals[entry] = normal_at_p(data, geometry, op, curr_point);
      if (j == 0 || i == n)
        continue;

      auto face_entry = 2 * entry - 2 * (i + 1) - i;
      faces[face_entry] = {entry, entry + n - i, entry - 1};
      if (j != n - i)
        faces[face_entry + 1] = {entry, entry + n - i + 1, entry + n - i};
    }
  }

  return {pos, faces, normals};
}

vector<float> interpolate_field(const shape_data &data,
                                const shape_geometry &geometry,
                                const shape_op &op, const vector<float> &f,
                                const int tid, const int k, const bool lerp) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      if (!lerp)
        field[entry] = field_at_p(data, geometry, op, f, curr_point);
      else
        field[entry] = field_at_p_lerp(data, geometry, op, f, curr_point);
    }
  }

  return field;
}
vector<vector<vec3f>> interpolate_tangent_basis(const shape_data &data,
                                                const shape_geometry &geometry,
                                                const shape_op &op,
                                                const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<vector<vec3f>>(N, vector<vec3f>(2));
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      auto [is_vert, kv] = point_is_vert(curr_point);
      if (is_vert) {
        auto [g, xuv] = metric_tensor(data, geometry, op,
                                      data.triangles[curr_point.face][kv]);
        field[entry] = {xuv[0], xuv[1]};
      } else {
        auto [g, xuv, T] = metric_tensor(data, geometry, op, curr_point);
        field[entry] = {xuv[0], xuv[1]};
      }
    }
  }

  return field;
}
vector<vec3f> interpolate_eigendecompositions(const shape_data &data,
                                              const shape_geometry &geometry,
                                              const shape_op &op, const int tid,
                                              const int k,
                                              const int eigen_entry = 0) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<vec3f>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      field[entry] = metric_tensor_eigen_decomposition(data, geometry, op,
                                                       curr_point)[eigen_entry];
    }
  }

  return field;
}
vector<float> interpolate_curvature(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const shape_op &op, const int tid,
                                    const int k, const bool lerp) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      if (lerp)
        field[entry] =
            gaussian_curvature_at_p_lerp(data, geometry, op, curr_point);
      else {
        auto [is_vert, kv] = point_is_vert(curr_point);
        if (is_vert)
          field[entry] = gaussian_curvature_at_vid(
              data, geometry, op, data.triangles[curr_point.face][kv]);
        else
          field[entry] =
              gaussian_curvature_at_p(data, geometry, op, curr_point);
      }

      assert(!isnan(field[entry]));
    }
  }

  return field;
}
vector<vec3f> interpolate_normals(const shape_data &data,
                                  const shape_geometry &geometry,
                                  const shape_op &op, const int tid,
                                  const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<vec3f>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      auto [is_vert, kv] = point_is_vert(curr_point);
      if (is_vert)
        field[entry] = data.normals[data.triangles[tid][kv]];
      else
        field[entry] = normal_at_p(data, geometry, op, curr_point);
    }
  }

  return field;
}
vector<vec3f> interpolate_gradient(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const vector<float> &f,
                                   const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<vec3f>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      auto [is_vert, kv] = point_is_vert(curr_point);
      if (is_vert)
        field[entry] = gradient_at_vid(data, geometry, op, f,
                                       data.triangles[curr_point.face][kv]);
      else
        field[entry] = gradient_at_p(data, geometry, op, f, curr_point);
    }
  }

  return field;
}
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
interpolate_gradient(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f,
                     const vector<mesh_point> &points) {
  auto N = points.size();
  auto field = vector<vec3f>(N);
  auto pos = vector<vec3f>(N);
  auto normals = vector<vec3f>(N);
  for (auto i = 0; i < N; ++i) {

    field[i] = gradient_at_p(data, geometry, op, f, points[i]);
    pos[i] = quadric_at_p(data, geometry, op, points[i]);
    normals[i] = normal_at_p(data, geometry, op, points[i]);
  }

  return {field, pos, normals};
}
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
interpolate_gradient_lerp(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f,
                          const vector<mesh_point> &points) {
  auto P = PCE_matrix(data.triangles, data.positions);
  auto dist_field = wrapper(f);
  Eigen::VectorXd PCE_grad = P * dist_field;
  auto N = points.size();
  auto field = vector<vec3f>(N);
  auto pos = vector<vec3f>(N);
  auto normals = vector<vec3f>(N);
  for (auto i = 0; i < N; ++i) {
    auto tid = points[i].face;
    auto field_at_tid = interpolate_triangle(
        f[data.triangles[tid][0]], f[data.triangles[tid][1]],
        f[data.triangles[tid][2]], vec2f{0.33, 0.33});
    field[i] = vec3f{(float)PCE_grad(3 * tid), (float)PCE_grad(3 * tid + 1),
                     (float)PCE_grad(3 * tid + 2)};
    field[i] /= (2 * std::sqrt(field_at_tid));
    pos[i] = eval_position(data.triangles, data.positions, points[i]);
    normals[i] = tid_normal(data.triangles, data.positions, points[i].face);
  }

  return {field, pos, normals};
}
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_within_triangles(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f,
                          const vector<mesh_point> &points, const bool lerp) {
  if (lerp)
    return interpolate_gradient_lerp(data, geometry, op, f, points);
  else
    return interpolate_gradient(data, geometry, op, f, points);
}
vector<float> interpolate_laplacian(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const shape_op &op, const vector<float> &f,
                                    const int tid, const int k,
                                    const bool lerp) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<float>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      if (lerp)
        field[entry] = laplacian_at_p_lerp(data, geometry, op, f, curr_point);
      else {
        auto [is_vert, kv] = point_is_vert(curr_point);
        if (is_vert)
          field[entry] = laplacian_at_vid(data, geometry, op, f,
                                          data.triangles[curr_point.face][kv]);
        else
          field[entry] = laplacian_at_p(data, geometry, op, f, curr_point);
      }
    }
  }

  return field;
}
std::tuple<vector<float>, vector<float>>
interpolate_hessian(const shape_data &data, const shape_geometry &geometry,
                    const shape_op &op, const vector<float> &f, const int tid,
                    const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto lambda0 = vector<float>(N);
  auto lambda1 = vector<float>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      auto H = Hessian_at_p(data, geometry, op, f, curr_point);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> dec(H);
      if (dec.info() != Eigen::Success)
        std::cerr << "Error in computing the eigenvalues of the metric tensor";

      auto lambda = dec.eigenvalues();
      lambda0[entry] = std::min(lambda(0), lambda(1));
      lambda1[entry] = std::max(lambda(0), lambda(1));
    }
  }

  return {lambda0, lambda1};
}
// vector<vec3f> interpolate_parametrization(const shape_data &data,
//                                           const shape_geometry &geometry,
//                                           const shape_op &op, const int tid,
//                                           const int k) {

//   auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
//   auto N = (n + 1) * (n + 2) / 2;
//   auto field = vector<vec3f>(N, zero3f);
//   for (auto i = 0; i <= n; ++i) {

//     for (auto j = 0; j <= n - i; ++j) {
//       auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
//       auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
//       auto [is_vert, kv] = point_is_vert(curr_point);
//       if (is_vert)
//         field[entry] = parametrization_at_vid(
//             data, geometry, op, data.triangles[curr_point.face][kv]);
//       else
//         field[entry] = parametrization_at_p(data, geometry, op, curr_point);
//     }
//   }

//   return field;
// }
std::tuple<vector<vec3f>, vector<float>> interpolate_gradient_and_magnitude(
    const shape_data &data, const shape_geometry &geometry, const shape_op &op,
    const vector<float> &f, const int tid, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto field = vector<vec3f>(N);
  auto magnitud = vector<float>(N);
  for (auto i = 0; i <= n; ++i) {

    for (auto j = 0; j <= n - i; ++j) {
      auto entry = N - (n - i + 2) * (n - i + 1) / 2 + j;
      auto curr_point = mesh_point{tid, {(float)j / n, (float)i / n}};
      auto [is_vert, kv] = point_is_vert(curr_point);
      if (is_vert)
        std::tie(field[entry], magnitud[entry]) = gradient_at_vid_w_magnitude(
            data, geometry, op, f, data.triangles[curr_point.face][kv]);
      else
        std::tie(field[entry], magnitud[entry]) =
            gradient_at_p_w_magnitude(data, geometry, op, f, curr_point);
    }
  }

  return {field, magnitud};
}
void export_subdivided_mesh(const shape_data &data,
                            const shape_geometry &geometry, const shape_op &op,
                            const int tid, const int k, const bool bumpy) {
  auto [points, faces] = (bumpy) ? bumpy_triangle(data, geometry, op, tid, k)
                                 : subdivide_tid(data, tid, k);
  string tri_name = "/Users/claudiomancinelli/Documents/GitHub/"
                    "InterpolatingQuadric/build/triangle_";
  if (bumpy)
    tri_name.append("bumpy_");
  tri_name.append(std::to_string(tid));
  tri_name.append(".obj");
  auto curr_shape = shape_data{};
  curr_shape.positions = points;
  curr_shape.triangles = faces;
  save_shape(tri_name, curr_shape);
}
void save_bumpy_mesh(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto F_per_tri = (int)pow(n, 2);
  auto positions = vector<vec3f>(N * (int)data.triangles.size());
  auto normals = vector<vec3f>(N * (int)data.triangles.size());
  auto faces = vector<vec3i>(F_per_tri * (int)data.triangles.size());
  for (auto i = 0; i < data.triangles.size(); ++i) {
    auto [points, tris, normal] =
        bumpy_triangle_w_normals(data, geometry, op, i, k);
    for (auto j = 0; j < N; ++j) {
      positions[N * i + j] = points[j];
      normals[N * i + j] = normal[j];
    }
    for (auto j = 0; j < F_per_tri; ++j) {
      faces[F_per_tri * i + j] = tris[j] + N * i;
    }
  }
  auto curr_shape = shape_data{};
  curr_shape.positions = positions;
  curr_shape.triangles = faces;
  curr_shape.normals = normals;
  string shape_name = "bumpy_bunny.obj";
  save_shape(shape_name, curr_shape);
}

void export_distance_field(const shape_data &data,
                           const shape_geometry &geometry, const shape_op &op,
                           const vector<float> &f, const int tid, const int k,
                           const bool lerp) {
  string name = "field_";
  name.append(std::to_string(tid));
  if (lerp)
    name.append("_lerp");

  auto field = interpolate_field(data, geometry, op, f, tid, k, lerp);

  export_field(field, name);
}
void export_curvature(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const int tid, const int k,
                      const bool lerp) {

  auto field = interpolate_curvature(data, geometry, op, tid, k, lerp);

  string name = "curvature_";
  if (lerp)
    name.append("lerp_");
  name.append(std::to_string(tid));
  export_field(field, name);
}
void export_gradient(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f, const int tid,
                     const int k) {
  auto grad = interpolate_gradient(data, geometry, op, f, tid, k);

  string name = "gradient_";
  name.append(std::to_string(tid));
  export_vector_field(grad, name);
}
void export_gradient_on_bumpy_mesh(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const vector<float> &f,
                                   const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto grad = vector<vec3f>(N * (int)data.triangles.size());
  for (auto i = 0; i < data.triangles.size(); ++i) {
    auto curr_grad = interpolate_gradient(data, geometry, op, f, i, k);
    for (auto j = 0; j < N; ++j) {
      grad[i * N + j] = curr_grad[j];
    }
  }

  string name = "gradient_on_bumpy_mesh";
  export_vector_field(grad, name);
}

void export_gradient(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f,
                     const vector<mesh_point> &PS) {
  auto [grad, pos, normals] = interpolate_gradient(data, geometry, op, f, PS);

  string name = "gradient_PS";
  export_vector_field(grad, name);
  string normal_name = "normals_PS";
  export_vector_field(normals, normal_name);
  string PS_name = "PS.obj";
  export_vector_field(pos, PS_name);
}
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_on_bumpy_mesh() {

  string name = "gradient_PS";
  std::ifstream f;
  f.open(name);
  assert(f.is_open());
  std::string dummy;
  int size = 0;
  f >> dummy >> size;
  auto grad = vector<vec3f>{};
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    grad.push_back(curr);
  }
  f.close();
  string PS_name = "PS.obj";
  f.open(PS_name);
  f >> dummy >> size;
  auto pos = vector<vec3f>{};
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    pos.push_back(curr);
  }
  f.close();

  f.open("normals_PS");
  f >> dummy >> size;
  auto normals = vector<vec3f>{};
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    normals.push_back(curr);
  }
  f.close();

  return {grad, pos, normals};
}
void export_laplacian(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f, const int tid,
                      const int k, const bool lerp) {
  auto lap = interpolate_laplacian(data, geometry, op, f, tid, k, lerp);

  string name = "laplacian_";
  if (lerp)
    name.append("lerp_");

  name.append(std::to_string(tid));
  export_field(lap, name);
}

void export_hessian(const shape_data &data, const shape_geometry &geometry,
                    const shape_op &op, const vector<float> &f, const int tid,
                    const int k) {
  auto [l0, l1] = interpolate_hessian(data, geometry, op, f, tid, k);

  string name = "Hessian_min_";
  name.append(std::to_string(tid));
  export_field(l0, name);
  name = "Hessian_max_";
  name.append(std::to_string(tid));
  export_field(l1, name);
}

void export_gradient_and_magnitude(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const vector<float> &f,
                                   const int tid, const int k) {
  auto [grad, mag] =
      interpolate_gradient_and_magnitude(data, geometry, op, f, tid, k);
  string name = "gradient_";
  name.append(std::to_string(tid));
  export_vector_field(grad, name);

  name = "gradient_mag_";
  name.append(std::to_string(tid));
  export_field(mag, name);
}

// std::tuple<float, float, int, int, int>
// check_gradient(const shape_data &data, const shape_geometry &geometry,
//                const shape_op &op, const vector<float> &f, const int k) {
//   // vector<vector<vec3f>> check(data.triangles.size());
//   vector<vector<float>> check(data.triangles.size());
//   vector<vector<vec3f>> pos(data.triangles.size());
//   auto err = flt_min;
//   auto avg_err = 0.f;
//   auto count = 0;
//   auto tid = 0;
//   auto vid = 0;
//   auto entry = 0;
//   auto adj_tid = 0;
//   auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
//   auto N = (n + 1) * (n + 2) / 2;
//   for (auto i = 0; i < data.triangles.size(); ++i) {

//     check[i] = interpolate_curvature(data, geometry, op, i, 4);
//     // check[i] = interpolate_eigendecompositions(data, geometry, op, i, 4);
//     // check[i] = interpolate_parametrization(data, geometry, op, i, 4);

//     // auto [grad, mag] =
//     //     interpolate_gradient_and_magnitude(data, geometry, op, f, i, 4);
//     // check[i] = mag;
//     //   check[i] = interpolate_gradient(data, geometry, op, f, i, 4);
//     // auto xuv = interpolate_tangent_basis(data, geometry, op, i, k);
//     // auto normals = vector<vec3f>(xuv.size());
//     // for (auto j = 0; j < xuv.size(); ++j) {
//     //   normals[j] = normalize(cross(xuv[j][0], xuv[j][1]));
//     // }
//     // check[i] = normals;
//     auto [pos_tid, faces] = subdivide_tid(data, i, 4);
//     pos[i] = pos_tid;
//   }

//   for (auto i = 0; i < data.triangles.size(); ++i) {
//     for (auto j = 0; j < 3; ++j) {
//       auto adj = geometry.adjacencies[i][j];
//       if (adj < i)
//         continue;

//       auto check_at_tid = vector<float>(n + 1);
//       // auto check_at_tid = vector<vec3f>(n + 1);
//       auto pos_at_tid = vector<vec3f>(n + 1);
//       if (j == 0) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_tid[s] = check[i][s];
//           pos_at_tid[s] = pos[i][s];
//         }
//       } else if (j == 1) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_tid[s] = check[i][N - (n - s + 2) * (n - s + 1) / 2 + n -
//           s]; pos_at_tid[s] = pos[i][N - (n - s + 2) * (n - s + 1) / 2 + n -
//           s];
//         }
//       } else if (j == 2) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_tid[s] = check[i][N - (s + 2) * (s + 1) / 2];
//           pos_at_tid[s] = pos[i][N - (s + 2) * (s + 1) / 2];
//         }
//       }
//       auto check_at_nei = vector<float>(n + 1);
//       // auto check_at_nei = vector<vec3f>(n + 1);
//       auto pos_at_nei = vector<vec3f>(n + 1);

//       auto offset = find(geometry.adjacencies[adj], i);

//       if (offset == 0) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_nei[s] = check[adj][s];
//           pos_at_nei[s] = pos[adj][s];
//         }
//       } else if (offset == 1) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_nei[s] =
//               check[adj][N - (n - s + 2) * (n - s + 1) / 2 + n - s];
//           pos_at_nei[s] = pos[adj][N - (n - s + 2) * (n - s + 1) / 2 + n -
//           s];
//         }
//       } else if (offset == 2) {
//         for (auto s = 0; s <= n; ++s) {
//           check_at_nei[s] = check[adj][N - (s + 2) * (s + 1) / 2];
//           pos_at_nei[s] = pos[adj][N - (s + 2) * (s + 1) / 2];
//         }
//       }

//       reverse(check_at_nei.begin(), check_at_nei.end());
//       reverse(pos_at_nei.begin(), pos_at_nei.end());

//       for (auto s = 0; s <= n; ++s) {
//         // auto value = length(check_at_nei[s] - check_at_tid[s]);
//         auto value = std::abs(check_at_nei[s] - check_at_tid[s]);
//         avg_err += value;
//         ++count;
//         if (value > err) {
//           err = value;
//           tid = i;
//           adj_tid = adj;
//           entry = s;
//         }

//         // if (std::abs(check_at_nei[s] - check_at_tid[s]) > err) {
//         //   err = std::abs(check_at_nei[s] - check_at_tid[s]);
//         //   tid = i;
//         //   entry = s;
//         // }
//         auto len = length(pos_at_nei[s] - pos_at_tid[s]);
//         if (len > 1e-10)
//           std::cerr << "These should be the same points";
//       }
//     }
//   }
//   avg_err /= count;
//   return {err, avg_err, tid, vid, entry};
// }
std::tuple<float, float, int, int, int>
check_gradient(const shape_data &data, const shape_geometry &geometry,
               const shape_op &op, const vector<float> &f, const int k) {
  vector<vector<vec3f>> check(data.triangles.size());
  // vector<vector<float>> check(data.triangles.size());
  vector<vector<vec3f>> pos(data.triangles.size());
  auto err = flt_min;
  auto avg_err = 0.f;
  auto count = 0;
  auto tid = 0;
  auto adj_tid = 0;
  auto entry = 0;
  Eigen::MatrixXd Q0;
  Eigen::MatrixXd Q1;

  for (auto i = 0; i < data.triangles.size(); ++i) {

    // auto points = internal_points(i, k);

    for (auto j = 0; j < 3; ++j) {

      auto adj = geometry.adjacencies[i][j];
      if (adj < i)
        continue;
      // auto vid = data.triangles[i][j];
      // auto Q = op.quadrics[vid];
      // for (auto s = 0; s < points.size(); ++s) {
      for (auto s = 0; s <= 10; ++s) {

        // auto T = reparametrization_matrix_from_vert_to_point(data, geometry,
        // Q,
        //                                                      vid, points[s]);
        // auto Qt = reparametrized_quadric_full(Q, T);

        auto bary = zero3f;
        bary[j] = (float)s / 10;
        bary[(j + 1) % 3] = 1 - bary[j];
        auto p0 = mesh_point{i, {bary.y, bary.z}};

        auto offset = find(geometry.adjacencies[adj], i);
        auto bary0 = zero3f;
        bary0[offset] = bary[(j + 1) % 3];
        bary0[(offset + 1) % 3] = bary[j];
        auto p1 = mesh_point{adj, {bary0.y, bary0.z}};
        // auto k0 = curvature_at_p_on_edge(
        //     data, geometry, op, p0,
        //     vec2i{data.triangles[i][j], data.triangles[i][(j + 1) % 3]});

        // auto k1 = curvature_at_p_on_edge(
        //     data, geometry, op, p1,
        //     vec2i{data.triangles[i][j], data.triangles[i][(j + 1) % 3]});

        auto n0 = normal_at_p_on_edge(
            data, geometry, op, p0,
            vec2i{data.triangles[i][j], data.triangles[i][(j + 1) % 3]});

        auto n1 = normal_at_p_on_edge(
            data, geometry, op, p1,
            vec2i{data.triangles[i][j], data.triangles[i][(j + 1) % 3]});

        //   // if (length(eval_position(data.triangles, data.positions, p0) -
        //   //            eval_position(data.triangles, data.positions, p1)) >
        //   //            1e-7)
        //   //   std::cerr << "These should be the same points";

        // auto value = std::abs(k0 - k1);
        // auto p_uv =
        //     coordinates_in_tangent_space_of_p(data, geometry, vid,
        //     points[s]);

        // auto n0 = evaluate_quadric_n(Q, p_uv);
        // auto n1 = evaluate_quadric_n(Qt, zero2f);
        auto value = length(n0 - n1);
        // auto value = std::abs(k0 - k1);
        assert(!isnan(value));
        if (value > err) {
          err = value;
          tid = i;
          adj_tid = j;
          entry = s;
        }
        // auto quadrics = quadric_at_p_on_edge(data, geometry, op, p0, p1);
        // auto q_w_rep =
        //     quadric_at_p_on_edge_w_repara(data, geometry, op, p0, p1);
        // for (auto h = 0; h < 3; ++h) {
        //   for (auto p = 0; p < 6; ++p) {
        //     auto value = std::abs(quadrics[0](h, p) - q_w_rep[0](h, p));
        //     avg_err += value;
        //     ++count;
        //     if (value > err) {
        //       err = value;
        //       tid = i;
        //       adj_tid = adj;
        //       entry = s;
        //       Q0 = quadrics[0];
        //       Q1 = quadrics[1];
        //     }
        //   }
        // }

        //   // if (std::abs(check_at_nei[s] - check_at_tid[s]) > err) {
        //   //   err = std::abs(check_at_nei[s] - check_at_tid[s]);
        //   //   tid = i;
        //   //   entry = s;
        //   // }
        //   // auto len = length(pos_at_nei[s] - pos_at_tid[s]);
        //   // if (len > 1e-10)
        //   //   std::cerr << "These should be the same points";
      }
    }
  }
  avg_err /= count;

  return {err, avg_err, tid, adj_tid, entry};
}

std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_inside_triangles(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f, const int k, const bool lerp,
                          const bool &squared) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto result = vector<vec3f>(N * data.triangles.size());
  auto pos = vector<vec3f>(N * data.triangles.size());
  auto normals = vector<vec3f>(N * data.triangles.size());
  for (auto i = 0; i < data.triangles.size(); ++i) {
    auto points = internal_points(i, k);
    auto n = tid_normal(data.triangles, data.positions, i);
    for (auto j = 0; j < N; ++j) {
      auto [is_vert, kv] = point_is_vert(points[j]);
      if (lerp)
        result[N * i + j] =
            gradient_at_p_lerp(data, geometry, op, f, points[j]);
      else {
        if (is_vert)
          result[N * i + j] = gradient_at_vid(
              data, geometry, op, f, data.triangles[points[j].face][kv]);
        else
          result[N * i + j] =
              gradient_at_p(data, geometry, op, f, points[j], squared);
        pos[N * i + j] =
            eval_position(data.triangles, data.positions, points[j]);
        normals[N * i + j] = n;
      }
    }
  }

  return {result, pos, normals};
}

std::tuple<vector<vec3f>, vector<vec3f>>
normals_inside_triangles(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const int k) {
  auto n = (k % 2) ? (int)pow(2, k) - 1 : (int)pow(2, k);
  auto N = (n + 1) * (n + 2) / 2;
  auto result = vector<vec3f>(N * data.triangles.size());
  auto pos = vector<vec3f>(N * data.triangles.size());
  for (auto i = 0; i < data.triangles.size(); ++i) {
    auto points = internal_points(i, k);
    for (auto j = 0; j < N; ++j) {
      result[N * i + j] = normal_at_p(data, geometry, op, points[j]);
      pos[N * i + j] = eval_position(data.triangles, data.positions, points[j]);
    }
  }

  return {result, pos};
}

vector<vec3f> compute_glyph(const shape_data &data,
                            const shape_geometry &geometry, const shape_op &op,
                            const mesh_point &p) {
  auto [g, xuv, Tr] = metric_tensor(data, geometry, op, p);
  std::cout << "\n";
  printf("curr point is {%d,{%f,%f}}", p.face, p.uv.x, p.uv.y);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> dec(g);
  if (dec.info() != Eigen::Success)
    std::cerr << "Error in computing the eigenvalues of the metric tensor";

  auto lambda = dec.eigenvalues();
  auto U = dec.eigenvectors();
  if (lambda(0) < 0 || lambda(1) < 0)
    std::cerr << "Error in computing the eigenvalues of the metric "
                 "tensor:negatives eigenvalues";

  auto theta = subdivide_angles(50);
  auto ellipse = vector<vec3f>(50);
  auto a = 1 / std::sqrt(lambda(0));
  auto b = 1 / std::sqrt(lambda(1));
  auto T = tranformation_matrix(vec2f{U(0, 0), U(1, 0)},
                                vec2f{U(0, 1), U(1, 1)}, zero2f);
  for (auto i = 0; i < 50; ++i) {
    auto ell = switch_reference_frame(
        T, vec2f{a / 100 * std::cos(theta[i]), b / 100 * std::sin(theta[i])});
    ellipse[i].x = ell.x;
    ellipse[i].y = ell.y;
    ellipse[i].z = 0;
    auto ep = get_basis_at_p(data.triangles, data.positions, p.face);

    auto Tp = tranformation_matrix(
        ep[0], ep[1], ep[2], eval_position(data.triangles, data.positions, p));
    ellipse[i] = switch_reference_frame(Tp.inverse(), ellipse[i]);
    // ellipse[i] /= 100;
  }

  return ellipse;
}
vector<vector<vec3f>> compute_glyph_cross(const shape_data &data,
                                          const shape_geometry &geometry,
                                          const shape_op &op,
                                          const mesh_point &p) {
  auto [g, xuv, Tr] = metric_tensor(data, geometry, op, p);

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> dec(g);
  if (dec.info() != Eigen::Success)
    std::cerr << "Error in computing the eigenvalues of the metric tensor";

  auto lambda = dec.eigenvalues();
  auto U = dec.eigenvectors();

  if (lambda(0) < 0 || lambda(1) < 0)
    std::cerr << "Error in computing the eigenvalues of the metric "
                 "tensor:negatives eigenvalues";
  auto result0 = vector<vec3f>(10);
  auto result1 = vector<vec3f>(10);
  auto a = 1 / std::sqrt(lambda(0));
  auto b = 1 / std::sqrt(lambda(1));
  auto u0 = a / 100 * vec2f{U(0, 0), U(1, 0)};
  auto u1 = b / 100 * vec2f{U(0, 1), U(1, 1)};
  for (auto i = -5; i < 5; ++i) {

    auto pos0 = (float)i / (10.f) * u0;
    auto pos1 = (float)i / (10.f) * u1;

    auto ep = get_basis_at_p(data.triangles, data.positions, p.face);
    auto Tp = tranformation_matrix(
        ep[0], ep[1], ep[2], eval_position(data.triangles, data.positions, p));
    result0[i + 5] = switch_reference_frame(Tp, vec3f{pos0.x, pos0.y, 0});
    result1[i + 5] = switch_reference_frame(Tp, vec3f{pos1.x, pos1.y, 0});
  }

  return vector<vector<vec3f>>{result0, result1};
}
std::tuple<vec3f, vec3f> compute_glyph_cross(const shape_data &data,
                                             const shape_geometry &geometry,
                                             const shape_op &op,
                                             const vector<float> &f,
                                             const mesh_point &p) {
  auto [g, xuv, Tr] = metric_tensor(data, geometry, op, p);
  auto H = Hessian_at_p(data, geometry, op, f, p);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> dec(H);
  if (dec.info() != Eigen::Success)
    std::cerr << "Error in computing the eigenvalues of the metric tensor";

  auto lambda = dec.eigenvalues();
  auto U = dec.eigenvectors();

  auto result0 = zero3f;
  auto result1 = zero3f;
  auto a = std::sqrt(std::abs(lambda(0)));
  auto b = std::sqrt(std::abs(lambda(1)));
  auto u0 = a * normalize(vec2f{U(0, 0), U(1, 0)});
  auto u1 = b * normalize(vec2f{U(0, 1), U(1, 1)});
  // auto n = normalize(cross(xuv[0], xuv[1]));
  // auto ep = get_basis_at_p(data.triangles, data.positions, p.face);
  // auto Tp = tranformation_matrix(xuv[0], xuv[1], n,
  //                                quadric_at_p(data, geometry, op, p));

  result0 = xuv[0] * u0.x + xuv[1] * u0.y;
  result1 = xuv[0] * u1.x + xuv[1] * u1.y;

  // result0 = switch_reference_frame(Tp.inverse(), vec3f{u0.x, u0.y, 0});
  // result1 = switch_reference_frame(Tp.inverse(), vec3f{u1.x, u1.y, 0});

  return {result0, result1};
}
void export_glyphs(const shape_data &data, const shape_geometry &geometry,
                   const shape_op &op, const vector<float> &f,
                   const vector<mesh_point> &p) {
  std::ofstream outfile0, outfile1, outfile2, outfile3;

  outfile0.open("alpha");
  outfile1.open("beta");
  outfile2.open("pos");
  outfile3.open("normals");

  for (auto i = 0; i < p.size(); ++i) {
    auto [alpha, beta] = compute_glyph_cross(data, geometry, op, f, p[i]);
    auto n = normal_at_p(data, geometry, op, p[i]);
    auto pos_p = quadric_at_p(data, geometry, op, p[i]);
    outfile0 << alpha.x << " " << alpha.y << " " << alpha.z << "\n";
    outfile1 << beta.x << " " << beta.y << " " << beta.z << "\n";
    outfile2 << pos_p.x << " " << pos_p.y << " " << pos_p.z << "\n";
    outfile3 << n.x << " " << n.y << " " << n.z << "\n";
  }

  outfile0.close();
  outfile1.close();
  outfile2.close();
  outfile3.close();
}
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>, vector<vec3f>>
import_glyph() {
  std::ifstream f;
  auto alpha = vector<vec3f>{};
  auto beta = vector<vec3f>{};
  auto pos = vector<vec3f>{};
  auto normals = vector<vec3f>{};
  f.open("alpha");
  assert(f.is_open());
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    alpha.push_back(curr);
  }
  f.close();

  f.open("beta");
  assert(f.is_open());
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    beta.push_back(curr);
  }
  f.close();

  f.open("pos");
  assert(f.is_open());
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    pos.push_back(curr);
  }
  f.close();

  f.open("normals");
  assert(f.is_open());
  while (!f.eof()) {
    vec3f curr;
    f >> curr[0] >> curr[1] >> curr[2];
    normals.push_back(curr);
  }
  f.close();

  return {alpha, beta, pos, normals};
}
vector<mesh_point> load_Poisson_sampling(const string &filename) {
  auto curr_shape = load_shape(filename);
  auto result = vector<mesh_point>(curr_shape.PS_faces.size());
  for (auto i = 0; i < result.size(); ++i) {
    result[i] =
        mesh_point{(int)curr_shape.PS_faces[i],
                   vec2f{curr_shape.PS_bary[i].y, curr_shape.PS_bary[i].z}};
  }
  return result;
}

void write_frames(const shape_data &data, const shape_geometry &geometry,
                  const int vid, const vector<uint> &vertices) {
  std::ofstream outfile;
  outfile.open("frames" + std::to_string(vid) + ".csv");
  auto e = polar_basis(data.triangles, data.positions, geometry.v2t,
                       data.normals, vid);
  auto n = data.normals[vid];
  outfile << e.x << " " << e.y << " " << e.z << "\n";
  outfile << n.x << " " << n.y << " " << n.z << "\n";
  for (auto i = 0; i < vertices.size(); ++i) {
    auto e = polar_basis(data.triangles, data.positions, geometry.v2t,
                         data.normals, vertices[i]);
    auto n = data.normals[vertices[i]];
    outfile << e.x << " " << e.y << " " << e.z << "\n";
    outfile << n.x << " " << n.y << " " << n.z << "\n";
  }

  outfile.close();
}
void write_neighborhood_coordinates(const shape_data &data,
                                    const shape_geometry &geometry,
                                    const vector<uint> &vertices,
                                    const int vid) {
  auto vert = data.positions[vid];
  auto nbr = one_ring(data.triangles, geometry.adjacencies, geometry.v2t, vid);
  auto coords = vector<vec2f>(vertices.size());
  for (auto i = 0; i < nbr.size(); ++i) {
    auto curr = nbr[i];
    auto it = std::find(vertices.begin(), vertices.end(), curr);
    if (it == vertices.end())
      std::cerr << "This should not happen";
    auto entry = std::distance(vertices.begin(), it);
    auto theta = geometry.angles[vid][i];
    auto r = length(data.positions[*it] - vert);
    coords[entry] = vec2f{r * std::cos(theta), r * std::sin(theta)};
  }
  std::ofstream outfile;
  outfile.open("coordinates" + std::to_string(vid) + ".csv");
  for (auto i = 0; i < coords.size(); ++i) {
    outfile << coords[i].x << " " << coords[i].y << "\n";
  }

  outfile.close();

  auto s = vertices.size();
  auto nei_coords = vector<vector<vec2f>>(s);

  for (auto i = 0; i < s; ++i) {
    auto curr_vert = vertices[i];
    auto nbr =
        one_ring(data.triangles, geometry.adjacencies, geometry.v2t, curr_vert);
    nei_coords[i].resize(3);
    for (auto j = 0; j < nbr.size(); ++j) {
      if (nbr[j] == vid) {
        auto theta = geometry.angles[curr_vert][j];
        auto r = length(data.positions[curr_vert] - vert);
        nei_coords[i][0] = vec2f{r * std::cos(theta), r * std::sin(theta)};
      }
      auto it = std::find(vertices.begin(), vertices.end(), nbr[j]);
      if (it == vertices.end())
        continue;
      auto entry = std::distance(vertices.begin(), it);

      if (entry == (s - 1 + i) % s) {
        auto theta = geometry.angles[curr_vert][j];
        auto r = length(data.positions[curr_vert] - data.positions[*it]);
        nei_coords[i][1] = vec2f{r * std::cos(theta), r * std::sin(theta)};
      } else if (entry == (i + 1) % s) {
        auto theta = geometry.angles[curr_vert][j];
        auto r = length(data.positions[curr_vert] - data.positions[*it]);
        nei_coords[i][2] = vec2f{r * std::cos(theta), r * std::sin(theta)};
      }
    }
  }
  outfile.open("coordinates_neighbors" + std::to_string(vid) + ".csv");
  for (auto i = 0; i < nei_coords.size(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      outfile << nei_coords[i][j].x << " " << nei_coords[i][j].y << " ";
    }
    outfile << "\n";
  }
  outfile.close();
}

void write_quadrics(const shape_op &op, const int vid,
                    const vector<uint> &vertices) {
  std::ofstream outfile;
  auto name = "quadrics_" + std::to_string(vid) + ".csv";
  outfile.open(name);

  auto c = op.quadrics[vid];
  outfile << c(0, 0) << " " << c(0, 1) << " " << c(0, 2) << " " << c(0, 3)
          << " " << c(0, 4) << " " << c(0, 5) << " \n";
  outfile << c(1, 0) << " " << c(1, 1) << " " << c(1, 2) << " " << c(1, 3)
          << " " << c(1, 4) << " " << c(1, 5) << " \n";
  outfile << c(2, 0) << " " << c(2, 1) << " " << c(2, 2) << " " << c(2, 3)
          << " " << c(2, 4) << " " << c(2, 5) << " \n";
  for (auto i = 0; i < vertices.size(); ++i) {
    c = op.quadrics[vertices[i]];
    outfile << c(0, 0) << " " << c(0, 1) << " " << c(0, 2) << " " << c(0, 3)
            << " " << c(0, 4) << " " << c(0, 5) << " \n";
    outfile << c(1, 0) << " " << c(1, 1) << " " << c(1, 2) << " " << c(1, 3)
            << " " << c(1, 4) << " " << c(1, 5) << " \n";
    outfile << c(2, 0) << " " << c(2, 1) << " " << c(2, 2) << " " << c(2, 3)
            << " " << c(2, 4) << " " << c(2, 5) << " \n";
  }

  outfile.close();
}

void write_radii(const shape_data &data, const shape_geometry &geometry,
                 const int vid, const vector<uint> &vertices) {
  std::ofstream outfile;
  outfile.open("radii" + std::to_string(vid) + ".csv");
  auto len =
      nbr_avg_edge_length(data.triangles, data.positions, geometry.v2t, vid);
  outfile << len << "\n";
  for (auto i = 0; i < vertices.size(); ++i) {
    len = nbr_avg_edge_length(data.triangles, data.positions, geometry.v2t,
                              vertices[i]);
    outfile << len << "\n";
  }
  outfile.close();
}

void export_quadrics(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const int vid,
                     const bool all_the_mesh) {
  auto star = geometry.v2t[vid];
  auto vertices = vector<uint>(star.size());
  auto pos = vector<vec3f>{};
  auto tri = vector<vec3i>{};
  if (all_the_mesh) {
    save_shape("mesh.obj", data);

    vertices =
        one_ring(data.triangles, geometry.adjacencies, geometry.v2t, vid);

  } else {
    pos.resize(star.size() + 1);
    tri.resize(star.size());
    pos[0] = data.positions[vid];
    for (auto i = 0; i < star.size(); ++i) {
      auto tid = star[i];
      auto k = find(data.triangles[tid], vid);
      auto curr_vid = data.triangles[tid][(k + 1) % 3];
      vertices[i] = curr_vid;
      pos[i + 1] = data.positions[curr_vid];
      auto curr_tri = zero3i;
      curr_tri[k] = 0;
      curr_tri[(k + 1) % 3] = i + 1;
      curr_tri[(k + 2) % 3] = (i == star.size() - 1) ? 1 : i + 2;
      tri[i] = curr_tri;
    }

    auto shape = shape_data{};
    shape.positions = pos;
    shape.triangles = tri;
    save_shape("neighborhood" + std::to_string(vid) + ".obj", shape);
  }

  write_quadrics(op, vid, vertices);
  write_frames(data, geometry, vid, vertices);
  write_neighborhood_coordinates(data, geometry, vertices, vid);
  write_radii(data, geometry, vid, vertices);
}
void export_local_coordinates(const shape_data &data,
                              const shape_geometry &geometry,
                              const shape_op &op, const int vid, const int vid0,
                              const int vid1) {
  auto coords = get_local_coordinates(data.triangles, data.positions,
                                      geometry.adjacencies, geometry.angles,
                                      geometry.v2t, vid, vid0, vid1);

  std::ofstream outfile;
  outfile.open("coordinates" + std::to_string(vid1) + ".csv");
  for (auto i = 0; i < coords.size(); ++i) {
    outfile << coords[i].x << " " << coords[i].y << "\n";
  }
}
