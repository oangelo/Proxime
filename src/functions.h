#include <vector>
#include <math.h>
#include <pthread.h>
#include <iostream>
#include <map>
#include <string>



class functions_capsule {
public:
  typedef std::pair<std::string, unsigned> name_item;
  typedef std::map<std::string,unsigned> items; 

  virtual void set(double &t, std::vector<double> & variables, std::vector<double> & parameters)=0;
  const double get_result(unsigned i) const;
  const std::vector<double> & get_result() const;
  const unsigned size() const;
 
  //map of the variables and index associated to them
  items variable_name_index,parameter_name_index;     
protected:
  std::string func_name;
  std::vector<double> __result;
};


class rossler_func : public functions_capsule {
public:
  rossler_func();
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);

 protected:
  double dx();
  double dy();
  double dz();
  
  double X,Y,Z;
  double a,b,c;

};

class Jacobian_rossler_func : public functions_capsule {
public:
  Jacobian_rossler_func(){__result.clear();__result.resize(3);};
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);
protected:
  double dx();
  double dy();
  double dz();
  
  double X,Y,Z;
  double a,b,c,X_fiducial,Y_fiducial,Z_fiducial;
};

class lorenz_func : public functions_capsule {
public:
  lorenz_func();
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);
protected:
  double dx();
  double dy();
  double dz();
  
  double X,Y,Z;
  double sigma,gamma,beta;
};

class Jacobian_lorenz_func : public functions_capsule {
public:
  Jacobian_lorenz_func(){__result.clear();__result.resize(3);};
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);
protected:
  double dx();
  double dy();
  double dz();
  
  double X,Y,Z;
  double sigma,gamma,beta,X_fiducial,Y_fiducial,Z_fiducial;
};

class double_pendulum_func : public functions_capsule {
public:
  double_pendulum_func(){__result.clear();__result.resize(4);};
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);
protected:
  double dTheta1();
  double dTheta2();
  double dOmega1();
  double dOmega2();
  
  double theta1,theta2,omega1,omega2;
  double l1,l2,m1,m2,g;

};

class jacobian_double_pendulum_func : public functions_capsule {
public:
  jacobian_double_pendulum_func(){__result.clear();__result.resize(4);};
  void set(double &t, std::vector<double> & variables, std::vector<double> & parameters);
protected:
  void Matrix_Jacob(double theta1, double theta2,double omega1, double omega2);
  double JdTheta1();
  double JdTheta2();
  double JdOmega1();
  double JdOmega2();

  double theta1,theta2,omega1,omega2;
  double l1,l2,m1,m2,g;
  double Jacobian[4][4];
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


