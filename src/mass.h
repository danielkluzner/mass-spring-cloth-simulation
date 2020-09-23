#ifndef MASS_H
#define MASS_H

#include "CGL/CGL.h"
#include "CGL/vector2D.h"

using namespace CGL;

struct Mass {
  Mass(Vector3D position, double mass, bool pinned)
      : start_position(position), position(position), last_position(position),
        mass(mass), pinned(pinned) {}

  double mass;
  bool pinned;

  Vector3D start_position;
  Vector3D position;

  // explicit Verlet integration

  Vector3D last_position;

  // explicit Euler integration

  Vector3D velocity;
  Vector3D forces;
};

#endif /* MASS_H */
