#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <algorithm>
#include <numeric>
#include <vector>

#include "CGL/vector2D.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"

#include "mass.h"
#include "cloth.h"
#include "spring.h"

typedef std::array<Vector3D, 256> vArr;


namespace CGL {
    int matrix_to_list(int i, int j, int num_nodes){
        return (i*num_nodes) + j;
    }

    Cloth::Cloth(Vector2D start, Vector2D end, int num_nodes, float node_mass, double k, vector<int> pinned_nodes)
    {
//generates in x-z plane
//note: face indices seem right
//GOES FROM LEFT TO RIGHT, THEN DOWN AND REPEATS
          vector<Mass *> m;
          vector<Spring *> s;
          vector<vector<int>> f;
          double step_size_x = (end.x - start.x) / ((num_nodes - 1) * 1.0);
          double step_size_y = (end.y - start.y) / ((num_nodes - 1) * 1.0);
          for (int i = 0; i < num_nodes; i++){
              for(int j = 0; j < num_nodes; j++){
                  Vector3D position = Vector3D(start.x + (i * 1.0 * step_size_x), 0, start.y + (j * 1.0 * step_size_y));
                  Mass* mass = new Mass(position, node_mass, false);
                  int s0 = matrix_to_list(i, j, num_nodes);
                  m.push_back(mass);  
                  if(i + 1 < num_nodes && j + 1 < num_nodes){
                      //add in faces
                      int v0 = matrix_to_list(i, j, num_nodes);
                      int v1 = matrix_to_list(i, j + 1, num_nodes);
                      int v2 = matrix_to_list(i + 1, j, num_nodes);
                      int v3 = matrix_to_list(i + 1, j + 1, num_nodes);
                      f.push_back({v0, v1, v2});

                      f.push_back({v3, v1, v2});
                  } 
                  if( i == 0){
                      //add line/top row of springs in
                      if(j > 0){
                          Spring* spring = new Spring(m.at(j - 1), mass, k, j - 1, s0);
                          s.push_back(spring);
                      }
                  } 
                  else{
                      //add connections to previous row
                      if(j == 0){
                          //add prev row same col, prev row next col
                          int s1 = matrix_to_list(i - 1, j, num_nodes);
                          int s2 = matrix_to_list(i - 1, j + 1, num_nodes);
                          Spring* spring1 = new Spring(m.at(s1), mass, k, s1, s0);
                          Spring* spring2 = new Spring(m.at(s2), mass, k, s2, s0);
                          s.push_back(spring1);
                          s.push_back(spring2);
                      }
                      else if ( j == num_nodes - 1){
                          //add prev row prev col, prev row same col, same row prev col
                          int s1 = matrix_to_list(i, j - 1, num_nodes);
                          int s2 = matrix_to_list(i - 1, j - 1, num_nodes);
                          int s3 = matrix_to_list(i - 1, j, num_nodes);
                          Spring* spring1 = new Spring(m.at(s1), mass, k, s1, s0);
                          Spring* spring2 = new Spring(m.at(s2), mass, k, s2, s0);
                          Spring* spring3 = new Spring(m.at(s3), mass, k, s3, s0);
                          s.push_back(spring1);
                          s.push_back(spring2);  
                          s.push_back(spring3);                    
                      }
                      else{
                          //add same row prev col, prev row prev col, prev row same col, prev row next col 
                          int s1 = matrix_to_list(i, j - 1, num_nodes);
                          int s2 = matrix_to_list(i - 1, j - 1, num_nodes);
                          int s3 = matrix_to_list(i - 1, j, num_nodes); 
                          int s4 = matrix_to_list(i - 1, j + 1, num_nodes);
                          Spring* spring1 = new Spring(m.at(s1), mass, k, s1, s0);
                          Spring* spring2 = new Spring(m.at(s2), mass, k, s2, s0);
                          Spring* spring3 = new Spring(m.at(s3), mass, k, s2, s0);
                          Spring* spring4 = new Spring(m.at(s4), mass, k, s4, s0);
                          s.push_back(spring1); 
                          s.push_back(spring2); 
                          s.push_back(spring3);
                          s.push_back(spring4);
                      }
                  }
              }
          }
          this->masses = m;
          this->faces = f;
          this->springs = s;
//        Comment-in this part when you implement the constructor
         for (auto &i : pinned_nodes) {
              masses[i]->pinned = true;
         }
    }

    void Cloth::simulateEuler(double delta_t, Vector3D gravity)
    {
        this->computeJacobians();
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            this->calculateForces(s);    
        }
        for (auto &m : masses)
        {   
            if (!m->pinned)
            {
                m->forces += gravity;
            }
        }

        vArr dv = CG(delta_t, 1000, 0.0001);
//        std::cout<<dv<<std::endl;
        for (int i = 0; i < masses.size(); i ++)
        {
            Mass* m = masses.at(i);
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->velocity = m->velocity + dv[i];
                m->position = m->position + m->velocity * delta_t;
            }

            // Reset all forces on each mass
            m->forces = Vector3D(0, 0, 0);
        }
    }
    
    void Cloth::outputOBJ(string file_name){
        ofstream MyFile(file_name); 
        for( auto & m : this->masses){
             MyFile <<"v "<<m->position.x << " " << m->position.y <<" "<< m->position.z <<std::endl;
        }
        for(auto & f: this->faces){
            MyFile<< "f " <<f.at(0) + 1<<" "<< f.at(1) + 1<< " "<< f.at(2) + 1<<"\n";
        }
        MyFile.close();
    }
    
    void Cloth::calculateForces(Spring* s){
        //Hooke's Law
        Vector3D diff = (s->m1)->position - (s->m2)->position;
        double length = diff.norm();
        Vector3D hooke = s->k * diff / length * (length - s->rest_length);
        (s->m1)->forces -= hooke;
        (s->m2)->forces += hooke;

        //Damping force
        Vector3D vDiff = (s->m1)->velocity - (s->m2)->velocity;
        (s->m1)->forces -= vDiff * Spring::kd;
        (s->m2)->forces += vDiff * Spring::kd;
    }

  void Cloth::computeJacobians(){
    for (auto &s : springs){
      Vector3D dx = (s->m1)->position - (s->m2)->position;
      Matrix3x3 dxtdx, I;
      dxtdx = outer(dx, dx);
      I = Matrix3x3::identity();

      double l = dx.norm();
      if(l != 0.0) l = 1.0 / l;

      dxtdx = dxtdx * l * l;

      Matrix3x3 temp = dxtdx;
      dxtdx += (I - temp) * (1.0 - s->rest_length * l);
      s->Jx = dxtdx * s->k;
      s->Jv = Matrix3x3::identity();
      s->Jv = s->Jv * Spring::kd;
    }
  }

  void Cloth::multiplyDfDx(vArr src, vArr& dst){
//    dst = new Vector3D[masses.size()];
    for(int i = 0; i < masses.size(); i ++){
      dst[i] = Vector3D(0, 0, 0);
    }
    for (auto &s : springs){
      Vector3D temp = s->Jx * (src[(s->n)[0]] - src[(s->n)[1]]);
      dst[(s->n)[0]] -= temp;
      dst[(s->n)[1]] += temp;
    }
  }

  void Cloth::multiplyA(vArr src, vArr &dst, double delta_t){

//HAVE NOT TESTED
      this->multiplyDfDx(src, dst);
      for(int i = 0; i < masses.size(); i ++){
      dst[i] = dst[i] * -1.0 * delta_t * delta_t;
      }
      for (auto &s : springs){
          Vector3D temp = delta_t * s->Jv * (src[(s->n)[0]] - src[(s->n)[1]]);
          //changed signs below cuz minus
          dst[(s->n)[0]] += temp;
          dst[(s->n)[1]] -= temp;
      }
      for(int i = 0; i < this->masses.size(); i++){
          dst[i] = dst[i] + this->masses.at(i)->mass * src[i];
      }
  }
  vArr Cloth::calculateB(double delta_t){
      
      vArr f0;
      vArr v0;
      vArr r;
      
      for (int i = 0; i<256; i ++){
          v0[i]=masses[i]->velocity;
          f0[i]=masses[i]->forces;
      }
      multiplyDfDx(v0,r);
      for (int i = 0; i<256; i ++){
          r[i] = r[i] * delta_t;
      }
      
      transform(f0.begin(), f0.end(),r.begin(),r.begin(),plus<Vector3D>());
      
      for (int i = 0; i<256; i ++){
          r[i] = r[i] * delta_t;
          
      }
      return r;
  }

  double dot_helper( vArr a,  vArr b){
     double r = inner_product(a.begin(), a.end(),b.begin(), 0.0, std::plus<double>(), [](Vector3D l, Vector3D r){ return dot(l,r);});
     return r;
  }

  vArr multByScalar(const vArr& input, double scalar){
    vArr temp;
    for(int i = 0; i < 256; i++){
      temp[i] = input[i] * scalar;
    }
    return temp;
  }

   vArr Cloth::CG(double delta_t, int iMax, double eps){
      vArr dv;
      for(int i = 0; i< this->masses.size(); i++){
          dv[i] = Vector3D(0,0,0);
      }
      int m = 0;
      vArr r = this->calculateB(delta_t);
      vArr d;
      d = r;
      vArr q;
      double epsNew = dot_helper(r, r);
      double eps0 = epsNew;
      while ( m < iMax && epsNew > eps * eps0){
          this->multiplyA(d, q, delta_t);
          double alpha = epsNew / dot_helper(d, q);
          vArr temp = multByScalar(d, alpha);
          transform(dv.begin(), dv.end(), temp.begin(), dv.begin(), std::plus<Vector3D>());
          temp = multByScalar(q, alpha);
          transform(r.begin(), r.end(), temp.begin(), r.begin(), std::minus<Vector3D>());
          double epsOld = epsNew;
          epsNew = dot_helper(r, r);
          double beta = epsNew / epsOld;
          temp = multByScalar(d, beta);
          transform(r.begin(), r.end(), temp.begin(), d.begin(), std::plus<Vector3D>());
          m++;
      }
      return dv;
  }

}
