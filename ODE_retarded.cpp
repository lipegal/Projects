#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <lapack.h>
#include <matrixtypes.h>
#include <iomanip>
#include <complex.h>
#include <fstream>
#include <iostream>
#include <cmath>
using namespace ula;
using namespace std;

//Runge-Kutta-4 function.
ComplexMatrix rk4(ComplexMatrix GR,double dt, ComplexMatrix (*func)(ComplexMatrix, double, double), double t ){
ComplexMatrix k1(4,4);
ComplexMatrix k2(4,4);
ComplexMatrix k3(4,4);
ComplexMatrix k4(4,4);
k1 = dt*func(GR, dt, t);
k2 = dt*func(GR + 0.5*dt*k1, dt, t + dt/2 );
k3 = dt*func(GR + 0.5*dt*k2, dt, t + dt/2 );
k4 = dt*func(GR + dt*k3, dt, t + dt );
return GR+((k1+2*k2+2*k3+k4)/6);
}

//eta value definition
double eta = 0.0001;

//Function of the matrix containing e_m elements (diagonalized hamiltonian).
ComplexMatrix E(){
ComplexMatrix E(3,3);
E(0,0) = 12410.000; E(0,1) = 0.;        E(0,2) = 0.;
E(1,0) = 0.;        E(1,1) = 12529.380; E(1,2) = 0.;
E(2,0) = 0.;        E(2,1) = 0.;        E(2,2) = 12209.924;
return E;
}

//Declaration of the inicial Greens retarded function matrix.
ComplexMatrix initial_GR(){
ComplexMatrix GR(3,3);
GR(0,0) = 1.; GR(0,1) = 0.; GR(0,2) = 0.;
GR(1,0) = 0.; GR(1,1) = 1.; GR(1,2) = 0.;
GR(2,0) = 0.; GR(2,1) = 0.; GR(2,2) = 1.;
return GR;
}

//Dirac delta declaration.
double diracDelta(double t,double tau, double N){
double peak = (1 / (sqrt(2 * M_PI) * N)) * exp(-(t - tau) * (t - tau) / (2 * N * N));
return peak;
}

//Definition of f(GR,t).

ComplexMatrix f(ComplexMatrix GR, double dt, double t){
std::complex<double> j(0.,1.);
double hbar = 1.;

//diracDelta parameters
//N is a big number to represent the peak behaviour of the dirac delta;
double N = 3.0;
//tau value.
double tau = 1.0;

ComplexMatrix identity = boost::numeric::ublas::identity_matrix<std::complex<double>>(3);

ComplexMatrix result = -j*(1/hbar)*diracDelta(t,tau,N)*identity-eta*(1/hbar)*GR-j*(1/hbar)*prod(E(),GR);

return result;
}

int main(){
//Definition of the .dat file that will have all the data.
ofstream coquito;
coquito.open("Green_retarded.dat", ios::out | ios::trunc);

//Definition of the time step and the max time.
double dt = 0.00001;
int t_max = 10;
double t = 0.00001;
//Definition of the initial greens function matrix.
ComplexMatrix GR_ev(3,3);
GR_ev = initial_GR();

//Rk4 iteration.
for (int i = 1; i<t_max; i++){
coquito << i <<"  "<< GR_ev(0,0).real() <<"  "<< GR_ev(1,1).real() <<"  "<< GR_ev(2,2).real() << endl;
cout<<rk4(GR_ev,dt,f,t)<<endl;
GR_ev = rk4(GR_ev,dt,f,i*dt);
t+=i*dt;
}

//Closing the file to write the data.
coquito.close();

return 0;
}
