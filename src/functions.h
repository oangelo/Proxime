#include <vector>
#include <math.h>
#include <pthread.h>
#include <iostream>
#include <map>
#include <string>


typedef double type_data;
typedef std::vector<type_data> type_container;

class functions_capsule {
public:
  functions_capsule();
  virtual ~functions_capsule(){};
  typedef std::pair<std::string, unsigned> name_item;
  typedef std::map<std::string,unsigned> items; 

  virtual void set(type_data &t, type_container & variables, type_container & parameters)=0;
  type_data get_result(unsigned i) const;
  const type_container & get_result() const;
  unsigned size() const;
 
  //map of the variables and index associated to them
  items variable_name_index,parameter_name_index;     
protected:
  std::string func_name;
  type_container __result;
};


class RosslerFunction : public functions_capsule {
public:
  RosslerFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);

 protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data a,b,c;

};

class Jacobian_RosslerFunction : public functions_capsule {
public:
  Jacobian_RosslerFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};

class LorenzFunction : public functions_capsule {
public:
  LorenzFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data sigma,gamma,beta;
};

class Jacobian_LorenzFunction : public functions_capsule {
public:
  Jacobian_LorenzFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dx();
  type_data dy();
  type_data dz();
  
  type_data X,Y,Z;
  type_data sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

class DoublePendulumFunction : public functions_capsule {
public:
  DoublePendulumFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  type_data dTheta1();
  type_data dTheta2();
  type_data dOmega1();
  type_data dOmega2();
  
  type_data theta1,theta2,omega1,omega2;
  type_data l1,l2,m1,m2,g;

};

class Jacobian_DoublePendulumFunction : public functions_capsule {
public:
  Jacobian_DoublePendulumFunction();
  void set(type_data &t, type_container & variables, type_container & parameters);
protected:
  void Matrix_Jacob(type_data theta1, type_data theta2,type_data omega1, type_data omega2);
  type_data JdTheta1();
  type_data JdTheta2();
  type_data JdOmega1();
  type_data JdOmega2();

  type_data theta1,theta2,omega1,omega2;
  type_data l1,l2,m1,m2,g;
  type_data Jacobian[4][4];
};

enum variables {
  V_THETA1, V_THETA2, V_OMEGA1, V_OMEGA2
};

enum parameters {
  P_L1, P_L2, P_M1, P_M2, P_G, P_OMEGA0, P_A, P_SIGMA, P_THETA1, P_THETA2, P_OMEGA1, P_OMEGA2
};

enum lorenz_enum {
  P_SIGMA_LORENZ,P_GAMMA_LORENZ,P_BETA_LORENZ,P_X,P_Y,P_Z
};


