#ifndef SPRING_H
#define SPRING_H

#include "CGL/CGL.h"
#include "mass.h"

using namespace std;

namespace CGL {

struct Spring {
  Spring(Mass *a, Mass *b, double k, int i, int j)
      : m1(a), m2(b), k(k), rest_length((a->position - b->position).norm()) {this->n[0] = i; this->n[1] = j;}
  
  ~Spring(){
    delete m1;
    delete m2;
  }

  double k;
  double rest_length;
  static constexpr double kd = 0.8;

  int n[2];
  Mass *m1;
  Mass *m2;
  Matrix3x3 Jx;
  Matrix3x3 Jv;
}; // struct Spring
}
#endif /* SPRING_H */
