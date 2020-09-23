#ifndef CLOTH_H
#define CLOTH_H

#include "CGL/CGL.h"
#include "mass.h"
#include "spring.h"
#include <array>

using namespace std;

namespace CGL {

class Cloth {
public:
  Cloth(Vector2D start, Vector2D end, int num_nodes, float node_mass, double k,
       vector<int> pinned_nodes);
  
  ~Cloth(){
    for(Spring *s : springs) delete s;
    for(Mass *m : masses) delete m;
  }

  void simulateEuler(double delta_t, Vector3D gravity);
  
  void calculateForces(Spring* s);
  void computeJacobians();
  void multiplyDfDx(std::array<Vector3D, 256> src, std::array<Vector3D, 256>& dst);
  void multiplyA(std::array<Vector3D, 256> src, std::array<Vector3D, 256>& dst, double delta_t); 
  std::array<Vector3D, 256> calculateB(double delta_t);
  std::array<Vector3D, 256>  CG(double delta_t, int iMax, double eps);
  void outputOBJ(string file_name);
   
  vector<Mass *> masses;
  vector<Spring *> springs;
  vector<vector<int>> faces;
}; // struct Cloth
}
#endif /* CLOTH_H */
