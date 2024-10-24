#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <lapack.h>
#include <matrixtypes.h>
#include <iomanip>
#include <complex.h>
#include <fstream>
#include <iostream>
using namespace ula;
using namespace std;

//Runge-Kutta-4 function.
ComplexMatrix rk4(ComplexMatrix p,double dt, ComplexMatrix (*func)(ComplexMatrix) ){
ComplexMatrix k1(4,4);
ComplexMatrix k2(4,4);
ComplexMatrix k3(4,4);
ComplexMatrix k4(4,4);
k1 = dt*func(p);
k2 = dt*func(p+0.5*k1);
k3 = dt*func(p+0.5*k2);
k4 = dt*func(p+k3);
return p+((k1+2*k2+2*k3+k4)/6);
}


//Function of the M1 matrix.
ComplexMatrix M1(ComplexMatrix p){
//a is used to test the whole found equations (1) or just the literature part (0).
double a= 0.;
//Definition of constant values implemented in the matrices.
double e1 = 0.188*12410., e2 = 0.188*12530., e3 = 0.188*12210.;
complex <double> t12 (0.188*-87.8,0), t23 (0.188*30.8,0);
//Definition of the M1 matrix with its components.
complex <double> j (0,1);
ComplexMatrix m1(4,4);
m1(0,0) = 0;
m1(0,1) = -e1*p(0,1)-a*(t12*p(0,2));
m1(0,2) = -e2*p(0,2)-a*(t12*p(0,1)+t23*p(0,3));
m1(0,3) = -e3*p(0,3)-a*(t23*p(0,2));
m1(1,0) = e1*p(1,0)+a*(t12*p(2,0));
m1(1,1) = t12*(p(2,1)-p(1,2));
m1(1,2) = (e1-e2)*p(1,2)+t12*(p(2,2)-p(1,1))-t23*p(1,3);
m1(1,3) = (e1-e3)*p(1,3)+t12*p(2,3)-t23*p(1,2);
m1(2,0) = e2*p(2,0)+a*(t12*p(1,0)+t23*p(3,0));
m1(2,1) = (e2-e1)*p(2,1)+t12*(p(1,1)-p(2,2))+t23*p(3,1);
m1(2,2) = t12*(p(1,2)-p(2,1))+t23*(p(3,2)-p(2,3));
m1(2,3) = (e2-e3)*p(2,3)+t12*p(1,3)+t23*(p(3,3)-p(2,2));
m1(3,0) = e3*p(3,0)+a*(t23*p(2,0));
m1(3,1) = (e3-e1)*p(3,1)-t12*p(3,2)+t23*p(2,1);
m1(3,2) = (e3-e2)*p(3,2)-t12*p(3,1)+t23*(p(2,2)-p(3,3));
m1(3,3) = t23*(p(2,3)-p(3,2));
return -j*m1;
}


//Function of the M2 matrix (Injection).
ComplexMatrix M2(ComplexMatrix p){
//Definition of the gamma (Injection) parameter.
double y_inj = 400;
//Definition of the M2 matrix with its components.
ComplexMatrix m2(4,4);
m2(0,0) = 2.*p(0,0);
m2(0,1) = p(0,1);
m2(0,2) = p(0,2);
m2(0,3) = p(0,3);
m2(1,0) = p(1,0);
m2(1,1) = -2.*p(0,0);
m2(1,2) = 0;
m2(1,3) = 0;
m2(2,0) = p(2,0);
m2(2,1) = 0;
m2(2,2) = 0;
m2(2,3) = 0;
m2(3,0) = p(3,0);
m2(3,1) = 0;
m2(3,2) = 0;
m2(3,3) = 0;
return -0.5*y_inj*m2;
}


//Function of the M3 matrix (Extraction).
ComplexMatrix M3(ComplexMatrix p){
//Definition of the gamma (Injection) parameter.
double y_ext = 0.1;
//Definition of the M3 matrix with its components.
ComplexMatrix m3(4,4);
m3(0,0) = -2.*p(3,3); 
m3(0,1) = 0;
m3(0,2) = 0;
m3(0,3) = p(0,3);
m3(1,0) = 0;
m3(1,1) = 0;
m3(1,2) = 0;
m3(1,3) = p(1,3);
m3(2,0) = 0;
m3(2,1) = 0;
m3(2,2) = 0;
m3(2,3) = p(2,3);
m3(3,0) = p(3,0);
m3(3,1) = p(3,1);
m3(3,2) = p(3,2);
m3(3,3) = 2.*p(3,3);
return -0.5*y_ext*m3;
}


//Function of the M4 matrix (Dephasing).
ComplexMatrix M4(ComplexMatrix p){
//Definition of the gamma (Dephasing) parameter.
double y_deph = 30;
//Definition of the M2 matrix with its components.
ComplexMatrix m4(4,4);
m4(0,0) = 0; 
m4(0,1) = p(0,1);
m4(0,2) = p(0,2);
m4(0,3) = p(0,3);
m4(1,0) = p(1,0);
m4(1,1) = 0;
m4(1,2) = p(1,2);
m4(1,3) = p(1,3);
m4(2,0) = p(2,0);
m4(2,1) = p(2,1);
m4(2,2) = 0;
m4(2,3) = p(2,3);
m4(3,0) = p(3,0);
m4(3,1) = p(3,1);
m4(3,2) = p(3,2);
m4(3,3) = 0;
return -y_deph*m4;
}

//Sum of all the matrices.
ComplexMatrix f(ComplexMatrix p){
return M1(p)+M2(p)+M3(p)+M4(p);
}

//Declaration of the inicial p state.
ComplexMatrix initial_p(){
ComplexMatrix p(4,4);
p(0,0) = 0.; p(0,1) = 0.; p(0,2) = 0.; p(0,3) = 0.;
p(1,0) = 0.; p(1,1) = 0.; p(1,2) = 0.; p(1,3) = 0.;
p(2,0) = 0.; p(2,1) = 0.; p(2,2) = 1.; p(2,3) = 0.;
p(3,0) = 0.; p(3,1) = 0.; p(3,2) = 0.; p(3,3) = 0.;
return p;
}


int main(){
//Definition of the .dat file that will have all the data.
ofstream coquito;
coquito.open("FMO.dat", ios::out | ios::trunc);

//Definition of the time step and the max time.
double dt = 0.001;
int t_max = 5000;

//Definition of the p-state matrix.
ComplexMatrix state(4,4);
state = initial_p();
//Rk4 iteration.
for (int i = 1; i<t_max; i++){
coquito << i <<"  "<< state(0,0).real() <<"  "<< state(1,1).real() <<"  "<< state(2,2).real() <<"  "<< state(3,3).real() << endl;
cout<<rk4(state,dt,&f)<<endl;
state = rk4(state,dt,&f);
}

//Closing the file to write the data.
coquito.close();

return 0;
}
