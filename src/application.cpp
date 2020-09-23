#include <iostream>

#include "application.h"
#include "cloth.h"

namespace CGL {

Application::Application(AppConfig config) { std::cout<<"constructing\n";i=0;this->config = config; }

Application::~Application() {}

void Application::init() {
  // Enable anti-aliasing and circular points.
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glPointSize(8);
  glLineWidth(4);

  glColor3f(1.0, 1.0, 1.0);

  ropeEuler = new Cloth(Vector2D(-48, -48), Vector2D(48, 48), 16, config.mass,
                       config.ks, {0, 15});
}

void Application::render() {
  for (int i = 0; i < config.steps_per_frame; i++) {
    ropeEuler->simulateEuler(1 / config.steps_per_frame, config.gravity);
  }

  Cloth *rope;

  for (int i = 0; i < 2; i++) {
    if (i == 0) {
      glColor3f(0.0, 0.0, 1.0);
      rope = ropeEuler;
    } else {
      glColor3f(0.0, 1.0, 0.0);
    }

    glBegin(GL_POINTS);

    for (auto &m : rope->masses) {
      Vector3D p = m->position;
      glVertex2d(p.x, p.z);
    }

    glEnd();

    glBegin(GL_LINES);

    for (auto &s : rope->springs) {
      Vector3D p1 = s->m1->position;
      Vector3D p2 = s->m2->position;
      glVertex2d(p1.x, p1.z);
      glVertex2d(p2.x, p2.z);
    }

    // -------------------------------------------------------------
    // This line generates a .obj file for every frame.
    // The file contains the postiion (x, y, z) of each of the 256 masses
    // at the given timestep.
    // It also contains the indices of the three masses comprising
    // each individual triangle face.
    // -------------------------------------------------------------
    // At the moment it is commented out to prevent the generation
    // of hundreds of .obj files.
    // You may test it out by uncommenting it and seeing the .obj
    // files that are generated.
    // -------------------------------------------------------------

    // rope->outputOBJ("obj_" + std::to_string(this->i) + ".obj");

    this->i++;
    glEnd();

    glFlush();
  }
}

void Application::resize(size_t w, size_t h) {
  screen_width = w;
  screen_height = h;

  float half_width = (float)screen_width / 2;
  float half_height = (float)screen_height / 2;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-half_width, half_width, -half_height, half_height, 1, 0);
}

void Application::keyboard_event(int key, int event, unsigned char mods) {
  switch (key) {
  case '-':
    if (config.steps_per_frame > 1) {
      config.steps_per_frame /= 2;
    }
    break;
  case '=':
    config.steps_per_frame *= 2;
    break;
  }
}

string Application::name() { return "Rope Simulator"; }

string Application::info() {
  ostringstream steps;
  steps << "Steps per frame: " << config.steps_per_frame;

  return steps.str();
}
}
