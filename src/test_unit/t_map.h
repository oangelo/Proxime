TEST(logistic_map) {
  ofstream data_file;
  data_file.open("l_map_zoom.out");
  for(double a=3.8;a<3.9;a+=0.0001){
    logistic_map m(0.3,a);
    m.next(10000);
    for(unsigned i=0;i<50;i++){
      m.next();
      data_file << a <<" "<< m.get_variable(0) << std::endl;
    }
  }
  data_file.close();
}

TEST(logistic_map_lyapunov) {
  ofstream data_file;
  data_file.open("l_map_lyapunov.out");
  for(double a=2.4;a<4.0;a+=0.001){
    logistic_map m(0.5,a);
    data_file << a;
    m.next(10000);
    data_file <<" "<< m.lyapunov(200) ;
    data_file << std::endl;
  }
  data_file.close();
}

TEST(logistic_map_supertrack) {
  ofstream data_file;
  data_file.open("l_map_supertrac.out");
  for(double a=2.4;a<4.0;a+=0.001){
    logistic_map m(0.5,a);
    data_file << a;
    for(unsigned i=0;i<8;i++){
      m.next();
      data_file <<" "<< m.get_variable(0) ;
    }
    data_file << std::endl;
  }
  data_file.close();
}

TEST(logistic_map_series) {
  ofstream data_file;
  data_file.open("l_map_s_a3.555.out");
  double a=3.555;
    logistic_map m(0.3,a);
    m.next(10000);
    for(unsigned i=0;i<400;i++){
      m.next();
      data_file << m.get_variable(0) << std::endl;
    }
  data_file.close();
}

TEST(logistic_map_geo_series) {
  ofstream data_file;
  data_file.open("l_map_geo_a3.51.out");
  double a=3.51, aux;
    logistic_map m(0.82,a);
    data_file << m.get_variable(0) <<" "<< 0 << std::endl;
    for(unsigned i=0;i<400;i++){
      aux=m.get_variable(0);
      m.next();
      data_file << aux << " "<< m.get_variable(0) << std::endl;
      data_file << m.get_variable(0) << " "<< m.get_variable(0) << std::endl;
    }
  data_file.close();
}