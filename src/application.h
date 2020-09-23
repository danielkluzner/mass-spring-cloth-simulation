#ifndef CGL_APPLICATION_H
#define CGL_APPLICATION_H

// STL
#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// libCGL
#include "CGL/CGL.h"
#include "CGL/osdtext.h"
#include "CGL/renderer.h"

#include "cloth.h"

using namespace std;

namespace CGL {

struct AppConfig {
  AppConfig() {
    // Rope config variables
    mass = 6;
    ks = 100;

    // Environment variables
    gravity = Vector3D(0, -9.8, 0);
    steps_per_frame = 23;
  }

  float mass;
  float ks;

  float steps_per_frame;
  Vector3D gravity;
};

class Application : public Renderer {
public:
  Application(AppConfig config);
  ~Application();

  void init();
  void render();
  void resize(size_t w, size_t h);

  std::string name();
  std::string info();

  void keyboard_event(int key, int event, unsigned char mods);
  // void cursor_event(float x, float y);
  // void scroll_event(float offset_x, float offset_y);
  // void mouse_event(int key, int event, unsigned char mods);

private:
  AppConfig config;

  Cloth *ropeEuler;
  Cloth *ropeVerlet;

  size_t screen_width;
  size_t screen_height;
  int i;
}; // class Application

} // namespace CGL

#endif // CGL_APPLICATION_H
