#include <gsl/gsl_integration.h>
#include <tuple>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <omp.h>
using namespace std;

//Constant definition
const double Dt = 0.001;
const double Max_time = 3.0;
const double KB_T = 0.1;
const double Epsilon_0 = 0.0;
const double Deltaa = 5.0;
const double Delta_L = 10.0;
const double Delta_R = 0.0;
const double Gamma_L = 0.5;
const double Gamma_R = 0.5;
const double U = 10.0;
const double Mu_left = 0.0;
const double Mu_right = 0.0;
const complex<double> im(0.0, 1.0);
const double lenght = Max_time/Dt + 2.; 
//Pulse duration
double S = 3.0;
double Gamma = Gamma_L + Gamma_R;

//Delta(t) function
double Delta(double t, double s, double delta){
return (0.0 < t && t <= s) ? delta : 0.0;
}

//Delta_{L/R}(t) function
pair<double, double> Delta_L_R(double t, double s, double delta_L, double delta_R) {
  double Del_L = (0.0 < t && t <= s) ? delta_L : 0.0;
  double Del_R = (0.0 < t && t <= s) ? delta_R : 0.0;
return make_pair(Del_L, Del_R);
}

//A_(L/R)^(0/U) function
tuple<vector<complex<double>>,vector<complex<double>>,vector<complex<double>>, vector<complex<double>>> A_L_R_0_U(double max_time, double epsilon, double s, double delta, double delta_L, double delta_R, double dt) {

    complex<double> A_L_0_dt = 1. / (epsilon - Epsilon_0 + 0.5 * im * Gamma);
    complex<double> A_L_U_dt = 1. / (epsilon - Epsilon_0 - U + 0.5 * im * Gamma);
    complex<double> A_R_0_dt;
    complex<double> A_R_U_dt;

    vector<complex<double>> A_L_0_array, A_L_U_array, A_R_0_array, A_R_U_array;

    A_L_0_array.push_back(A_L_0_dt);
    A_L_U_array.push_back(A_L_U_dt);
    A_R_0_array.push_back(A_L_0_dt);
    A_R_U_array.push_back(A_L_U_dt);

  for (double t = 0.0; t <= max_time; t += dt) {
    double t_mid = t + 0.5 * dt;
    double t_mid2 = t + 0.75 * dt;

    A_L_0_dt = A_L_0_dt * exp(im * (epsilon - Epsilon_0) * dt) * exp(- im * dt *
      (Delta(t_mid, s, delta) - Delta_L_R(t_mid, s, delta_L, delta_R).first)) * exp(-Gamma * 0.5 * dt) - im* dt * exp(im * epsilon * 0.5 * dt) *
      exp(-im * dt * 0.5 * Epsilon_0) * exp(-im * dt * 0.5 * (Delta(t_mid2, s, delta) -
      Delta_L_R(t_mid2, s, delta_L, delta_R).first)) * exp(-Gamma * 0.25 * dt);

    A_L_0_array.push_back(A_L_0_dt);

    A_L_U_dt = A_L_U_dt * exp(im * (epsilon - Epsilon_0 - U) * dt) * exp(-im * dt *
      (Delta(t_mid, s, delta) - Delta_L_R(t_mid, s, delta_L, delta_R).first)) * exp(-Gamma * 0.5 * dt) - im* dt * exp(im * epsilon * 0.5 * dt) *
      exp(-im * dt * 0.5 * (Epsilon_0 + U)) * exp(-im * dt * 0.5 * (Delta(t_mid2,s, delta) -
      Delta_L_R(t_mid2, s, delta_L, delta_R).first)) * exp(-Gamma * 0.25 * dt);

    A_L_U_array.push_back(A_L_U_dt);

    A_R_0_dt = A_L_0_dt * exp(im * (epsilon - Epsilon_0) * dt) * exp(-im * dt *
      (Delta(t_mid, s, delta) - Delta_L_R(t_mid, s, delta_L, delta_R).second)) * exp(-Gamma * 0.5 * dt) - im * dt * exp(im * epsilon * 0.5 * dt) *
      exp(-im * dt * 0.5 * Epsilon_0) * exp(-im * dt * 0.5 * (Delta(t_mid2, s,delta) -
      Delta_L_R(t_mid2, s, delta_L, delta_R).second)) * exp(-Gamma * 0.25 * dt);

    A_R_0_array.push_back(A_R_0_dt);

    A_R_U_dt = A_L_U_dt * exp(im * (epsilon - Epsilon_0 - U) * dt) * exp(-im * dt * (Delta(t_mid, s, delta) - 
      Delta_L_R(t_mid, s, delta_L, delta_R).second)) * exp(-Gamma * 0.5 * dt) - im * dt * exp(im * epsilon * 0.5 * dt) *
      exp(-im * dt * 0.5 * (Epsilon_0 + U)) * exp(-im* dt * 0.5 * (Delta(t_mid2,s, delta) -
      Delta_L_R(t_mid2, s, delta_L, delta_R).second)) * exp(-Gamma * 0.25 * dt);

    A_R_U_array.push_back(A_R_U_dt);
  }
return make_tuple(A_L_0_array, A_L_U_array, A_R_0_array, A_R_U_array);
}

//Fermi functions
pair<double, double> fermi_functions(double energy, double mu_left, double mu_right, double kB_T) {
  double fermi_left = 1.0 / (1.0 + exp((energy - mu_left) / (kB_T)));
  double fermi_right = 1.0 / (1.0 + exp((energy - mu_right) / (kB_T)));
return make_pair(fermi_left, fermi_right);
}

//First integrand
double Integrand1(double epsilonn, void * params){
  double j= *(int *) params;
  double Integrand1 =(1 / (2 * M_PI)) * ( (fermi_functions(epsilonn, Mu_left, Mu_right, KB_T).first *
    pow(abs((get<0>(A_L_R_0_U(Max_time, epsilonn, S, Deltaa, Delta_L, Delta_R, Dt)))[j]), 2) + fermi_functions(epsilonn, Mu_left, Mu_right, KB_T).second *
    pow(abs((get<2>(A_L_R_0_U(Max_time, epsilonn, S, Deltaa, Delta_L, Delta_R, Dt)))[j]), 2)));
return Integrand1;
}

//Second integrand
double Integrand2(double epsilonn, void * params){
  double k= *(int *) params;
  double Integrand2 =(1 / (2 * M_PI)) * ((fermi_functions(epsilonn, Mu_left, Mu_right,KB_T).first *
    pow(abs((get<1>(A_L_R_0_U(Max_time, epsilonn, S, Deltaa, Delta_L, Delta_R, Dt)))[k]), 2) + fermi_functions(epsilonn, Mu_left, Mu_right, KB_T).second *
    pow(abs((get<3>(A_L_R_0_U(Max_time, epsilonn, S, Deltaa, Delta_L, Delta_R, Dt)))[k]), 2)));
return Integrand2;
}

//Integration function (Gauss-Kronrod method)
vector<double> integration_gauss(){
  //Assign memory workspace for gsl 
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

  //Integration results and errors declaration
  double n_0, err0;
  double n_U, errU;
  vector<double> t_values;
  //Time array generation
  for (double l=0 ; l<=lenght*Dt; l+=Dt){
    t_values.push_back(l);
  };
  //Declaration of array with sigma values
  vector<double> n_sigma;

  //Memory limits
  size_t numU = 100, num0 = 100;

    ofstream coquito;
    coquito.open("n_sigma.dat", ios::out | ios::trunc);

  //Paralelized for cycle
  //#pragma omp parallel for
  for (int i=0; i <=lenght; i++){
    int j = i;
    int k = i;

    //First integrand integration
    gsl_function I1;
    I1.function = &Integrand1;
    I1.params = &j;
    gsl_integration_qng(&I1,-50.0,50.0,1e-2,1e-4, &n_0, &err0, &num0);
  
    //Second integrand integration 
    gsl_function I2;
    I2.function = &Integrand2;
    I2.params =&k;
    gsl_integration_qng(&I2,-50.0,50.0,1e-2,1e-4, &n_U, &errU, &numU);
  
    //Sigma values assignation
    double n_sigmaa = (0.5*n_0) /(1 + 0.5*n_0 - 0.5*n_U);

    //Adding sigma values to array
    n_sigma.push_back(n_sigmaa);

    //quitar esta parte (es debugging)


  coquito<< t_values[i] <<" "<< n_sigmaa << endl;
  }
  coquito.close();
  //Freeing workspace
  gsl_integration_workspace_free(w);
  return n_sigma;
};


tuple<vector<complex<double>>,vector<complex<double>>,vector<complex<double>>, vector<complex<double>>> B_L_R_0_U(vector<double> times, vector<double> n_values, double epsilon, double s, double delta, double delta_L, double delta_R, double dt) {
  complex<double> B_L_0_dt = (1.0 - n_values[0]) / (epsilon - Epsilon_0 + 0.5 * im * Gamma );
  complex<double> B_L_U_dt = n_values[0] / (epsilon - Epsilon_0 - U + 0.5 * im * Gamma );
  complex<double> B_R_0_dt;
  complex<double> B_R_U_dt;

  vector<complex<double>> B_L_0_array, B_L_U_array, B_R_0_array, B_R_U_array;

  B_L_0_array.push_back(B_L_0_dt);
  B_L_U_array.push_back(B_L_U_dt);
  B_R_0_array.push_back(B_L_0_dt);
  B_R_U_array.push_back(B_L_U_dt);

  double current_t = 0.0;
  double t_mid, t_mid2;

  for(size_t a=0; a < times.size(); ++a) {
      double current_t = times[a];

      t_mid = current_t + 0.5 * dt;
      t_mid2 = current_t + 0.75 * dt;

    B_L_0_dt = B_L_0_dt * exp(im * (epsilon - Epsilon_0) * dt) * exp(-im * dt * (Delta(t_mid, s, delta)-
              Delta_L_R(t_mid, s, delta_L, delta_R).first)) * exp(-im * Gamma * 0.5 * dt) - im * dt * exp(im * epsilon * 0.5 * dt) *
              (0.5 * (2.0 - n_values[a] - n_values[a+1])) * exp(-im * dt * 0.5 * Epsilon_0) * exp(-im * dt * 0.5 * (Delta(t_mid2, s, delta) -
              Delta_L_R(t_mid2, s, delta_L, delta_R).first)) * exp(-im * Gamma * 0.25 * dt ) ;

    B_L_U_dt = B_L_U_dt * exp(im * (epsilon - Epsilon_0 - U) * dt) * exp(-im * dt * (Delta(t_mid, s, delta) -
              Delta_L_R(t_mid, s, delta_L, delta_R).first)) * exp(-im * Gamma * 0.5 * dt) - im * dt * exp(im * epsilon * 0.5 * dt) *
              (0.5 * (n_values[a] + n_values[a+1])) * exp(-im * dt * 0.5 * (Epsilon_0 + U)) * exp(-im * dt * 0.5 * (Delta(t_mid2, s, delta) -
              Delta_L_R(t_mid2, s, delta_L, delta_R).first)) * exp(-im * Gamma * 0.25 * dt );

    B_R_0_dt = B_L_0_dt * exp(im * (epsilon - Epsilon_0) * dt) * exp(-im * dt * (Delta(t_mid, s, delta) -
              Delta_L_R(t_mid, s, delta_L, delta_R).second)) * exp(-im * Gamma * 0.5 * dt) - im * dt * exp(im * epsilon * 0.5 * dt) *
              (0.5 * (2.0 - n_values[a] - n_values[a+1])) * exp(-im * dt * 0.5 * Epsilon_0) * exp(-im * dt * 0.5 * (Delta(t_mid2, s, delta) -
              Delta_L_R(t_mid2, s, delta_L, delta_R).second)) * exp(-im * Gamma * 0.25 * dt );

    B_R_U_dt = B_L_U_dt * exp(im * (epsilon - Epsilon_0 - U) * dt) * exp(-im * dt * (Delta(t_mid, s, delta) -
              Delta_L_R(t_mid, s, delta_L, delta_R).second)) * exp(-im * Gamma * 0.5 * dt) - im * dt * std::exp(im * epsilon * 0.5 * dt) *
              (0.5 * (n_values[a] + n_values[a+1])) * exp(-im * dt * 0.5 * (Epsilon_0 + U)) * exp(-im * dt * 0.5 * (Delta(t_mid2, s, delta) -
              Delta_L_R(t_mid2, s, delta_L, delta_R).second)) * exp(-im * Gamma * 0.25 * dt );

    B_L_0_array.push_back(B_L_0_dt);
    B_L_U_array.push_back(B_L_U_dt);
    B_R_0_array.push_back(B_R_0_dt);
    B_R_U_array.push_back(B_R_U_dt);

    current_t += dt;
    };

  return make_tuple(B_L_0_array, B_L_U_array, B_R_0_array, B_R_U_array);
};

int main() {

  //File assignation for writing n_sigma results
  //ofstream coquito;
  //coquito.open("n_sigma.dat", ios::out | ios::trunc);
  
  vector<double> t_values;
  //Time array generation
  for (double l=0 ; l<=30.002; l+=Dt){
    t_values.push_back(l);
  };
  
  vector<double> n_sigma_val = integration_gauss();

  //for (int i=0; i <=30002; i++){
    //Writing sigma value with time in file
    //coquito <<t_values[i]<<" "<< n_sigma_val[i]<< endl;
  //};
/*To test the values of A_(L/R)^(0/U)
 
double Epsilon_1 = 30.0;

auto gen_a = A_L_R_0_U(Max_time, Epsilon_1, S, Deltaa, Delta_L, Delta_R, Dt);

//Selection of the array elements of the tuple
vector<complex<double>> A_L_0_array = get<0>(gen_a);
vector<complex<double>> A_L_U_array = get<1>(gen_a);
vector<complex<double>> A_R_0_array = get<2>(gen_a);
vector<complex<double>> A_R_U_array = get<3>(gen_a);

//Printing the first 20 results
cout << " t" << "              " << "A_L_O" << "                   " << "A_L_U" << "                     " << "A_R_O" << "                    " << "A_R_U" << endl;
for (size_t i = 0; i < 20; ++i) {
cout << t_values[i] << "    " << A_L_0_array[i] << "    " << A_L_U_array[i] << "    " << A_R_0_array[i] << "    " << A_R_U_array[i] << endl;
};

*/

//


return 0;
};

