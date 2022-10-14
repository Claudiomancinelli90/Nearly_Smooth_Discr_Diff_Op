#include <stdio.h>
#include <thread>
#include <vector>
#include <yocto/yocto_cli.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
using namespace std;
#include "app.h"
#include <utils/drawing_circle.h>
#include <utils/utilities.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_opengl.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

//
#include "editing.h"
#include "playback.h"

void set_common_uniforms(const App &app, const ogl_program *program) {
  auto &view = app.matrices.view;
  auto &projection = app.matrices.projection;
  set_uniform(program, "frame", identity4x4f);
  set_uniform(program, "view", view);
  set_uniform(program, "projection", projection);
  set_uniform(program, "eye", app.camera->frame.o);
  set_uniform(program, "envlight", (int)app.envlight);
  set_uniform(program, "gamma", app.shade_params.gamma);
  set_uniform(program, "exposure", app.shade_params.exposure);
  // set_uniform(program, "size", app.line_size);
  if (app.scene->environments.size()) {
    auto &env = app.scene->environments.front();
    if (env->envlight_diffuse)
      set_uniform(program, "envlight_irradiance", env->envlight_diffuse, 6);
    if (env->envlight_specular)
      set_uniform(program, "envlight_reflection", env->envlight_specular, 7);
    if (env->envlight_brdflut)
      set_uniform(program, "envlight_brdflut", env->envlight_brdflut, 8);
  }
}

void draw_scene(const App &app, const vec4i &viewport) {
  clear_ogl_framebuffer(vec4f{0, 0, 0, 1});

  // Draw mesh and environment.
  draw_scene(app.scene, app.camera, viewport, app.shade_params);

  if (app.show_points) {
    auto program = &app.shaders.at("points");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "size", 3.0f * 0.0015f * app.line_size);

    set_uniform(program, "color", vec3f{0, 1, 0});

    auto draw_mesh_point = [&](const mesh_point &point, const vec3f &color) {
      if (point.face == -1)
        return;
      auto p = eval_position(app.mesh.triangles, app.mesh.positions, point);
      static auto shape = ogl_shape{};
      set_points_shape(&shape, {p});
      set_uniform(program, "color", color);
      draw_shape(&shape);
    };

    for (int i = 0; i < app.eval_points.size(); i++) {
      draw_mesh_point(app.eval_points[i], {1, 1, 1});
    }
  }

  if (app.temp_levels > 0)
    draw_shape(app.temp_points[app.temp_levels]);
  auto camera_aspect = (float)viewport.z / (float)viewport.w;
  auto camera_yfov =
      camera_aspect >= 0
          ? (2 * yocto::atan(app.camera->film /
                             (camera_aspect * 2 * app.camera->lens)))
          : (2 * yocto::atan(app.camera->film / (2 * app.camera->lens)));
  auto view = frame_to_mat(inverse(app.camera->frame));
  auto projection = perspective_mat(
      camera_yfov, camera_aspect, app.shade_params.near, app.shade_params.far);

  if (app.gpu_shapes.find("edges") != app.gpu_shapes.end())
    gpu::draw_shape(app.gpu_shapes.at("edges"), app.gpu_shaders.at("points"),
                    gpu::Uniform("color", vec3f{0, 0, 0}));
  gpu::set_point_size(10);
  if (app.gpu_shapes.find("selected_points") != app.gpu_shapes.end()) {

    gpu::draw_shape(
        app.gpu_shapes.at("selected_points"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));
  }

  if (app.gpu_shapes.find("vector_field") != app.gpu_shapes.end())
    gpu::draw_shape(
        app.gpu_shapes.at("vector_field"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{0, 0, 1}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));

  if (app.gpu_shapes.find("vector_field_2") != app.gpu_shapes.end())
    gpu::draw_shape(
        app.gpu_shapes.at("vector_field_2"), app.gpu_shaders.at("points"),
        gpu::Uniform("color", vec3f{1, 0, 0}),
        gpu::Uniform("frame", identity4x4f), gpu::Uniform("view", view),
        gpu::Uniform("projection", projection));

  if (app.show_edges) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{0, 0, 0});
    draw_shape(&app.edges_shape);
  }

  if (app.show_branches) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);

    set_uniform(program, "color", vec3f{0, 0, 1});
    set_uniform(program, "size", 0.005f);
    draw_shape(&app.branches_shape);
  }

  if (app.show_cobranches) {
    auto program = &app.shaders.at("lines");
    bind_program(program);
    set_common_uniforms(app, program);
    set_uniform(program, "color", vec3f{1, 0, 0});
    set_uniform(program, "size", 0.01f);
    draw_shape(&app.co_branches_shape);
  }
}

inline void sleep(int ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

inline bool is_releasing(gui_button button) {
  return button.state == gui_button::state::releasing;
}
inline bool is_down(gui_button button) {
  return button.state == gui_button::state::down ||
         button.state == gui_button::state::pressing;
}

bool process_camera_move(App &app, const gui_input &input) {
  auto &camera = *app.camera;

  auto update_camera_frame = [&](frame3f &frame, float &focus, bool rotating,
                                 bool panning, bool zooming) {
    auto last_pos = input.mouse_last;
    auto mouse_pos = input.mouse_pos;
    auto mouse_left = is_down(input.mouse_left);
    auto mouse_right = is_down(input.mouse_right);
    // handle mouse and keyboard for navigation
    if (mouse_left) {
      auto dolly = 0.0f;
      auto pan = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left)
          rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right)
          dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_right)
          pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    } else if (mouse_right) {
      auto dolly = 0.0f;
      auto pan = zero2f;
      auto rotate = zero2f;
      if (rotating) {
        if (mouse_left)
          rotate = (mouse_pos - last_pos) / 100.0f;
      }
      if (zooming) {
        if (mouse_right)
          dolly = (mouse_pos.y - last_pos.y) / 100.0f;
      }
      if (panning) {
        if (mouse_right)
          pan = (mouse_pos - last_pos) * focus / 200.0f;
      }
      pan.x = -pan.x;
      rotate.y = -rotate.y;
      update_turntable(frame, focus, rotate, dolly, pan);
    }
  };

  if (is_down(input.mouse_left)) {
    update_camera_frame(camera.frame, app.camera_focus, true, false, false);

    return true;
  }
  if (is_down(input.mouse_right)) {
    update_camera_frame(camera.frame, app.camera_focus, false, true, false);

    return true;
  }

  // Zoom-in/out by scrolling;
  float zoom = input.scroll.y * 0.1;
  if (zoom != 0) {
    update_turntable(camera.frame, app.camera_focus, zero2f, zoom, zero2f);
    return true;
  }

  return false;
}

bool process_user_input(App &app, const gui_input &input) {
  //  static bool yyy = false;

  auto mouse = input.mouse_pos;
  auto size = vec2f{(float)input.window_size.x, (float)input.window_size.y};
  mouse = vec2f{2 * (mouse.x / size.x) - 1, 1 - 2 * (mouse.y / size.y)};

  if (input.modifier_shift && is_releasing(input.mouse_left)) {
    // Here if pressing, but not clicked on an existing control point.
    auto point = intersect_mesh(app, mouse);
    if (point.face != -1) {
      app.control_points.push_back(point);
      std::cout << "\n";
      printf("selected point is: {%d,{%f,%f}}", point.face, point.uv.x,
             point.uv.y);
      // auto vid = forced_vert_from_point(app.mesh.triangles, point);
      // auto [nbr, tetas, lens] =
      //     uniform_stencil(app.mesh, app.topology, vid, 36);

      // add_points_shape(app, nbr, 0.0015, {1, 0, 0});
      app.source_shape_entry = app.added_points.size();
      add_points_shape(app, {point}, 0.0015, {0, 0, 1});

      return true;
    }
  }
  if (process_camera_move(app, input)) {
    update_camera_info(app, input);
    return false;
  }
  return false;
}

void update_app(App &app, const gui_input &input) {
  // process_gui_input(app, input); TODO(giacomo)

  if (is_active(app.widget))
    return;

  app.window_size = input.window_size;

  process_user_input(app, input);

  auto tasks = vector<vec2i>{};
}

void draw(const gui_input &input, void *data) {
  auto &app = *(App *)data;
  app.started = true;

  update_camera_info(app, input);

  // Do everything
  auto &t = app.playback_tick;
  if (app.playback && t < app.input_record.size()) {
    update_app(app, app.input_record[t]);
    t += 1;
  } else {
    update_app(app, input);
  }

  draw_scene(app, input.framebuffer_viewport);

  auto widget = app.widget;
  begin_widget(widget, "Interpolation");
  static vector<string> field_names = {"Distance", "Gradient", "Curvature",
                                       "Laplacian", "Hessian"};

  draw_bullet_text(widget, "I/O");
  draw_checkbox(widget, "LERP", app.lerp);
  draw_checkbox(widget, "No Geometry", app.no_geometry);
  draw_checkbox(widget, "Bumpy Geometry", app.bumpy);
  draw_combobox(widget, "Exported Field", app.exported_field, field_names);
  continue_line(widget);
  if (draw_button(widget, "Export")) {
    switch (app.exported_field) {
    case fields::dist: {
      if (app.field.size() > 0) {
        string name = "distance_field";
        export_field(app.field, name);
      }
    }; break;

    case fields::grad: {
      if (app.field.size() > 0) {
        string name = "gradient_field";

        export_field(app.field, name);
      }

    }; break;
    }
  }
  if (draw_button(widget, "Export Interpolated")) {
    switch (app.exported_field) {
    case fields::dist: {
      if (app.field.size() > 0) {
        for (auto i = 0; i < app.mesh.triangles.size(); ++i) {

          export_distance_field(app.mesh, app.topology, app.operators,
                                app.field, i, 4, app.lerp);
          if (!app.no_geometry)
            export_subdivided_mesh(app.mesh, app.topology, app.operators, i, 4,
                                   app.bumpy);
        }
      }
    }; break;

    case fields::grad: {
      if (app.field.size() > 0) {
        if (app.PS.size() > 0)
          export_gradient(app.mesh, app.topology, app.operators, app.field,
                          app.PS);
        else {
          for (auto i = 0; i < app.mesh.triangles.size(); ++i) {
            export_gradient(app.mesh, app.topology, app.operators, app.field, i,
                            4);
            if (!app.no_geometry)
              export_subdivided_mesh(app.mesh, app.topology, app.operators, i,
                                     4, app.bumpy);
          }
        }
      }

    }; break;

    case fields::lap: {

      for (auto i = 0; i < app.mesh.triangles.size(); ++i) {
        export_laplacian(app.mesh, app.topology, app.operators, app.field, i, 4,
                         app.lerp);
      }

    }; break;

    case fields::hess: {

      for (auto i = 0; i < app.mesh.triangles.size(); ++i) {
        export_hessian(app.mesh, app.topology, app.operators, app.field, i, 4);
      }

    }; break;
    }
  }

  if (draw_button(widget, "Export Bumpy Mesh")) {
    save_bumpy_mesh(app.mesh, app.topology, app.operators, 4);
  }
  if (draw_button(widget, "Export Mesh")) {
    for (auto i = 0; i < app.mesh.triangles.size(); ++i) {

      export_subdivided_mesh(app.mesh, app.topology, app.operators, i, 4,
                             false);
    }
  }
  static vector<string> method_names = {"VTP", "Heat", "Graph"};
  draw_bullet_text(widget, "Metric Tensor");

  if (draw_button(widget, "Export Glyph")) {

    // app.control_points.push_back(mesh_point{0, vec2f{0.33, 0.33}});
    if (app.PS.size() > 0 && app.field.size() > 0) {
      export_glyphs(app.mesh, app.topology, app.operators, app.field, app.PS);

      std::tie(app.alpha, app.beta, app.glyph_pos, app.glyph_normals) =
          import_glyph();
      app.alpha_shape = add_glyph_shape(
          app, app.alpha, app.glyph_pos, app.glyph_normals, app.vector_size,
          app.vector_thickness, {1, 0, 0}, app.lift_factor);
      app.beta_shape = add_glyph_shape(
          app, app.beta, app.glyph_pos, app.glyph_normals, app.vector_size,
          app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
  }
  if (draw_button(widget, "Load Glyph")) {

    std::tie(app.alpha, app.beta, app.glyph_pos, app.glyph_normals) =
        import_glyph();
    app.alpha_shape = add_glyph_shape(
        app, app.alpha, app.glyph_pos, app.glyph_normals, app.vector_size,
        app.vector_thickness, {1, 0, 0}, app.lift_factor);
    app.beta_shape = add_glyph_shape(
        app, app.beta, app.glyph_pos, app.glyph_normals, app.vector_size,
        app.vector_thickness, {0, 0, 1}, app.lift_factor);
  }

  draw_textinput(widget, "PS filename", app.PS_name);
  if (draw_button(widget, "Import PS")) {
    app.PS = load_Poisson_sampling(app.PS_name);
  }
  if (draw_button(widget, "Interpolate Normals")) {
    std::tie(app.vector_field, app.vector_field_pos) =
        normals_inside_triangles(app.mesh, app.topology, app.operators, 3);

    app.vector_field_shape = add_generic_vector_field_shape(
        app, app.vector_field, app.vector_field_pos, app.vector_field_normals,
        app.vector_size, app.vector_thickness, {1, 0, 0}, app.lift_factor);
  }

  if (draw_button(widget, "Export Curvature")) {

    for (auto i = 0; i < app.mesh.triangles.size(); ++i) {
      export_curvature(app.mesh, app.topology, app.operators, i, 4, app.lerp);
      if (!app.no_geometry)
        export_subdivided_mesh(app.mesh, app.topology, app.operators, i, 4,
                               app.bumpy);
    }
  }
  draw_separator(widget);
  if (draw_button(widget, "Compute Field")) {
    if (app.control_points.size() > 0) {
      if (app.control_points.size() == 1)
        app.field = exact_geodesic_distance(
            app.mesh.triangles, app.mesh.positions, app.control_points[0]);
      else
        app.field = compute_geodesic_distances(
            app.solver, app.mesh.triangles, app.mesh.positions,
            app.topology.adjacencies, app.control_points);
      // for (auto i = 0; i < app.field.size(); ++i) {
      //   app.field[i] = std::sin(10 * app.field[i]);
      // }
      // std::transform(app.field.begin(), app.field.end(), app.field.begin(),
      //                [](float lambda) { return lambda * lambda; });
    }
  }

  draw_separator(widget);
  draw_bullet_text(widget, "Gradient Field");
  if (draw_button(widget, "Compute Gradient Inside Triangles")) {
    if (app.control_points.size() > 0 && app.field.size() > 0) {
      if (app.PS.size() > 0)
        std::tie(app.vector_field, app.vector_field_pos,
                 app.vector_field_normals) =
            gradient_within_triangles(app.mesh, app.topology, app.operators,
                                      app.field, app.PS, app.lerp);
      else
        std::tie(app.vector_field, app.vector_field_pos,
                 app.vector_field_normals) =
            gradient_inside_triangles(app.mesh, app.topology, app.operators,
                                      app.field, 4, app.lerp, true);
      // app.vector_field2.resize(app.mesh.positions.size());
      // for (auto i = 0; i < app.mesh.positions.size(); ++i) {
      //   app.vector_field2[i] = gradient_at_vid(app.mesh, app.topology,
      //                                          app.operators, app.field, i);
      // }
      app.vector_field_shape = add_generic_vector_field_shape(
          app, app.vector_field, app.vector_field_pos, app.vector_field_normals,
          app.vector_size, app.vector_thickness, {1, 0, 0}, app.lift_factor);
      // app.vector_field_shape2 =
      //     add_vector_field_shape(app, app.vector_field2, app.vector_size,
      //                            app.vector_thickness, {0, 0, 1});
    }
  }
  if (draw_button(widget, "Load Gradient")) {

    std::tie(app.vector_field, app.vector_field_pos, app.vector_field_normals) =
        gradient_on_bumpy_mesh();
    app.vector_field_shape = add_generic_vector_field_shape(
        app, app.vector_field, app.vector_field_pos, app.vector_field_normals,
        app.vector_size, app.vector_thickness, {1, 0, 0}, app.lift_factor);
  }

  item_size(15);
  if (draw_button(widget, "Show Gradient")) {
    if (app.field.size() > 0 && app.PS.size() > 0) {

      app.vector_field = PCE_at_PS(app.mesh, app.field, app.PS);
      std::tie(app.vector_field_pos, app.vector_field_normals) =
          pos_and_normals_PS(app.mesh, app.PS);
      app.vector_field_shape = add_generic_vector_field_shape(
          app, app.vector_field, app.vector_field_pos, app.vector_field_normals,
          app.vector_size, app.vector_thickness, {1, 0, 0}, app.lift_factor);
      // std::cout << "here" << std::endl;

      // auto [err_max, avg_err, tid, adj_tid, entry] =
      //     check_gradient(app.mesh, app.topology, app.operators, app.field,
      //     4);
      // printf("Max error is %f, which has been computed at the edge shared by
      // "
      //        "triangle %d and triangle %d the is entry %d",
      //        err_max, tid, adj_tid, entry);
      // std::cout << "here" << std::endl;
      // std::cout << "\n";
      // std::cout << "first quadric is" << std::endl;
      // std::cout << Q0 << std::endl;
      // std::cout << "second quadric is" << std::endl;
      // std::cout << Q1 << std::endl;
      // auto vid = app.mesh.triangles[948][2];
      // auto [nbr, tetas, lens] =
      //     uniform_stencil(app.mesh, app.topology, vid, 36);
      // auto vid0 = 461;
      // auto vid1 = 485;
      // auto tid = 948;
      // auto len = 0.0568701439;
      // auto alpha = 1.62358284;
      // auto v = rot_vect(app.mesh.positions[vid0] - app.mesh.positions[vid],
      //                   tid_normal(app.mesh.triangles, app.mesh.positions,
      //                   tid), alpha);
      // auto path = straightest_path_from_vert_in_triangle(
      //     app.mesh, app.topology, tid, vid, len, normalize(v));
      // add_path_shape(app,
      //                polyline_pos(app.mesh.triangles, app.mesh.positions,
      //                path), 0.001, zero3f);
      // for (auto i = 24; i < 36; ++i) {
      //   add_points_shape(app, {nbr[i]}, 0.0015, {1, 0, 0});
      // }
      // add_points_shape(app, nbr, 0.0015, {1, 0, 0});

      // add_points_shape(app, {app.mesh.positions[vid]}, 0.0015, {0, 0, 0});
      // vid = app.mesh.triangles[948][1];
      // std::tie(nbr, tetas, lens) =
      //     uniform_stencil(app.mesh, app.topology, vid, 36);

      // add_points_shape(app, nbr, 0.0015, {0, 1, 0});
      // vid = app.mesh.triangles[948][2];
      // std::tie(nbr, tetas, lens) =
      //     uniform_stencil(app.mesh, app.topology, vid, 36);

      // add_points_shape(app, nbr, 0.0015, {0, 0, 1});
    }
  }
  if (draw_slider(widget, "vectors size", app.scale_factor, 0, 20)) {
    app.vector_size = 0.001 * app.scale_factor;
    if (app.vector_field_shape != nullptr) {
      update_generic_vector_field_shape(
          app.vector_field_shape->instance->shape, app.vector_field,
          app.vector_field_pos, app.vector_field_normals, app.vector_size,
          app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
    if (app.vector_field_shape2 != nullptr) {
      update_vector_field_shape(app.vector_field_shape2->instance->shape,
                                app.mesh, app.vector_field2, app.vector_size,
                                app.vector_thickness, {1, 0, 0});
    }
    if (app.alpha_shape != nullptr) {
      update_glyph_shape(app.alpha_shape->instance->shape, app.alpha,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {1, 0, 0}, app.lift_factor);
    }

    if (app.beta_shape != nullptr) {
      update_glyph_shape(app.beta_shape->instance->shape, app.beta,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
  }

  if (draw_slider(widget, "vectors thickness",
                  app.vector_thickness_scale_factor, 1, 10)) {
    app.vector_thickness = 0.0001 * app.vector_thickness_scale_factor;
    if (app.vector_field_shape != nullptr) {
      update_generic_vector_field_shape(
          app.vector_field_shape->instance->shape, app.vector_field,
          app.vector_field_pos, app.vector_field_normals, app.vector_size,
          app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
    if (app.vector_field_shape2 != nullptr) {
      update_vector_field_shape(app.vector_field_shape2->instance->shape,
                                app.mesh, app.vector_field2, app.vector_size,
                                app.vector_thickness, {1, 0, 0});
    }
    if (app.alpha_shape != nullptr) {
      update_glyph_shape(app.alpha_shape->instance->shape, app.alpha,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {1, 0, 0}, app.lift_factor);
    }

    if (app.beta_shape != nullptr) {
      update_glyph_shape(app.beta_shape->instance->shape, app.beta,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
  }

  if (draw_slider(widget, "lift vector", app.lift_factor, 0, 0.05)) {

    if (app.vector_field_shape != nullptr && app.lift_factor != 0) {

      update_generic_vector_field_shape(
          app.vector_field_shape->instance->shape, app.vector_field,
          app.vector_field_pos, app.vector_field_normals, app.vector_size,
          app.vector_thickness, {1, 0, 0}, app.lift_factor);
    }

    if (app.vector_field_shape2 != nullptr && app.lift_factor != 0) {
      for (auto i = 0; i < app.vector_field.size(); ++i) {
        app.vector_field2[i] += app.lift_factor * app.mesh.normals[i];
      }
      update_vector_field_shape(app.vector_field_shape2->instance->shape,
                                app.mesh, app.vector_field2, app.vector_size,
                                app.vector_thickness, {1, 0, 0});
    }

    if (app.alpha_shape != nullptr) {
      update_glyph_shape(app.alpha_shape->instance->shape, app.alpha,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {1, 0, 0}, app.lift_factor);
    }

    if (app.beta_shape != nullptr) {
      update_glyph_shape(app.beta_shape->instance->shape, app.beta,
                         app.glyph_pos, app.glyph_normals, app.vector_size,
                         app.vector_thickness, {0, 0, 1}, app.lift_factor);
    }
  }
  draw_textinput(widget, "Scene Name", app.exported_scene_name);
  if (draw_button(widget, "Export Scene")) {
    export_scene(app, app.mesh_material->color, app.exported_scene_name);
  }
  draw_checkbox(widget, "show edges", app.show_edges);
  draw_coloredit(widget, "Mesh Color", app.mesh_material->color);

  app.shade_params.faceted = app.show_edges;
  draw_separator(widget);
  if (draw_button(widget, " Reset")) {
    app.control_points.clear();
    for (auto path : app.added_paths)
      clear_shape(path->instance->shape);
    for (auto point : app.added_points)
      clear_shape(point->instance->shape);
    app.vector_field.clear();
    app.vector_field_pos.clear();
    app.added_paths.clear();
    app.added_points.clear();
  }

  end_widget(widget);
}

int main(int num_args, const char *args[]) {
  auto app = App();

  bool grad = false;
  bool lap = false;
  bool hess = false;
  int msaa = 1;

  auto cli = make_cli("bezier", "interactive viewer for mesh processing");
  add_argument(cli, "mesh", app.filename, "Model filenames");
  add_option(cli, "Gradient", grad, "Compute the Gradient Matrix");
  add_option(cli, "Laplacian", lap, "Compute the Laplacian Matric");
  add_option(cli, "Hessian", hess, "Compute the Hessian matrix");
  parse_cli_and_handle_errors(cli, vector<string>{args, args + num_args});

  // Load model and init bvh for fast click-intersection.
  if (!load_mesh(app.filename, app.mesh, app.topology, app.operators,
                 app.solver, app.dual_solver, app.error, grad, lap, hess))
    print_fatal(app.error);
  init_bvh(app);

  // Init window.
  auto win = new gui_window();
  win->msaa = msaa;
  init_window(win, {1080, 720}, "mesh viewer", true);
  win->user_data = &app;

  init_gpu(app, app.envlight);

  init_widget(app.widget, win);

  if (msaa > 1)
    set_ogl_msaa();

  run_ui(win, draw);

  // TODO(giacomo): delete app
  clear_window(win);
}
