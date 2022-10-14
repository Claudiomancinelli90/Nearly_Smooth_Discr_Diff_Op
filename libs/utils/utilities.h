#ifndef UTILITIES_H
#define UTILITIES_H

#include "Surazhsky/geodesic_algorithm_base.h"
#include "Surazhsky/geodesic_algorithm_exact.h"
#include <diff_geo/diff_geo.h>
#include <fstream>
#include <utils/mesh_io.h>
#include <yocto/yocto_sceneio.h>
vector<mesh_point> internal_points(const int tid, const int k);

// std::tuple<float, float, int, int, int>
// check_gradient(const shape_data &data, const shape_geometry &geometry,
//                const shape_op &op, const vector<float> &f, const int k);
std::tuple<float, float, int, int, int>
check_gradient(const shape_data &data, const shape_geometry &geometry,
               const shape_op &op, const vector<float> &f, const int k);

// vec3f parametrization_at_p(const shape_data &data,
//                            const shape_geometry &geometry, const shape_op
//                            &op, const mesh_point &p);

void export_subdivided_mesh(const shape_data &data,
                            const shape_geometry &geometry, const shape_op &op,
                            const int tid, const int k, const bool bumpy);

float field_at_p(const shape_data &data, const shape_geometry &geometry,
                 const shape_op &op, const vector<float> &f,
                 const mesh_point &p);
float laplacian_at_p(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f,
                     const mesh_point &p);
float laplacian_at_p_lerp(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f, const mesh_point &p);
vector<mesh_point> load_Poisson_sampling(const string &filename);

float field_at_p_lerp(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f,
                      const mesh_point &p);
float gaussian_curvature_at_p_lerp(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const mesh_point &p);
void save_bumpy_mesh(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, int k);

void export_glyphs(const shape_data &data, const shape_geometry &geometry,
                   const shape_op &op, const vector<float> &f,
                   const vector<mesh_point> &p);
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>, vector<vec3f>>
import_glyph();
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_on_bumpy_mesh();

vector<vec3f> PCE_at_PS(const shape_data &data, const vector<float> &f,
                        const vector<mesh_point> &PS);
std::tuple<vector<vec3f>, vector<vec3f>>
pos_and_normals_PS(const shape_data &data, const vector<mesh_point> &PS);
void export_gradient_on_bumpy_mesh(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const vector<float> &f,
                                   const int k);
void export_field(const vector<float> &field, const string &filename);

void export_distance_field(const shape_data &data,
                           const shape_geometry &geometry, const shape_op &op,
                           const vector<float> &f, const int tid, const int k,
                           const bool lerp);
void export_gradient(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f, const int tid,
                     const int k);
void export_gradient(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const vector<float> &f,
                     const vector<mesh_point> &PS);
void export_gradient_magnitude(const shape_data &data,
                               const shape_geometry &geometry,
                               const shape_op &op, const vector<float> &f,
                               const int tid, const int k);
void export_gradient_and_magnitude(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const vector<float> &f,
                                   const int tid, const int k);
void export_laplacian(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f, const int tid,
                      const int k, const bool lerp);
void export_hessian(const shape_data &data, const shape_geometry &geometry,
                    const shape_op &op, const vector<float> &f, const int tid,
                    const int k);
void export_curvature(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const int tid, const int k,
                      const bool lerp);
std::tuple<vector<float>, geodesic::Mesh, geodesic::GeodesicAlgorithmExact>
exact_geodesic_distance_surazhsky(const vector<vec3i> &triangles,
                                  const vector<vec3f> &positions,
                                  const mesh_point &source);
vector<vec3f> compute_glyph(const shape_data &data,
                            const shape_geometry &geometry, const shape_op &op,
                            const mesh_point &p);
std::tuple<vector<vec3f>, vector<vec3f>>
normals_inside_triangles(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const int k);
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_inside_triangles(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f, const int k, const bool lerp,
                          const bool &squared);
std::tuple<vector<vec3f>, vector<vec3f>, vector<vec3f>>
gradient_within_triangles(const shape_data &data,
                          const shape_geometry &geometry, const shape_op &op,
                          const vector<float> &f,
                          const vector<mesh_point> &points, const bool lerp);
vec3f gradient_at_vid(const shape_data &data, const shape_geometry &geometry,
                      const shape_op &op, const vector<float> &f,
                      const int &vid);
void export_quadrics(const shape_data &data, const shape_geometry &geometry,
                     const shape_op &op, const int vid,
                     const bool all_the_mesh = false);
vec3f normal_at_p_in_tangent_space(const shape_data &data,
                                   const shape_geometry &geometry,
                                   const shape_op &op, const int vid,
                                   const mesh_point &p);
vec3f normal_at_p_in_tid(const shape_data &data, const shape_geometry &geometry,
                         const shape_op &op, const mesh_point &p,
                         const int vid);
void export_local_coordinates(const shape_data &data,
                              const shape_geometry &geometry,
                              const shape_op &op, const int vid, const int vid0,
                              const int vid1);
#endif