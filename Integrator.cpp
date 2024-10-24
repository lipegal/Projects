//Used packages (No need for extra installations or makefiles).
#include <iomanip>
#include <fstream>
#include <iostream>
using namespace std;
//Definition of the time step (MUST BE THE SAME AS IN THE state_rk4.cpp that generated the data)
double dt=0.001;
//Function for counting file lenght.
int count_lines(string filename){
  int count=-1;
  string line;
  ifstream file(filename);
    while (getline(file,line))
    {
      count++;
    }
return count;
}
//Integration function (column from 0 to 3).
double simpson_integrator(string filename){
  //Opening of input file.
  ifstream in(filename);
  //File lenght (N).
  int N= count_lines(filename);
  //Definition of time and columns arrays.
  double index[N], c0[N], c1[N], c2[N], c3[N];
  //Reading and assignation of file columns to arrays.
  for(int i=0;i<N+1;i++){
    in >> index[i] >> c0[i] >> c1[i] >> c2[i] >> c3[i];
  }
  double sum,x;
   //YOU CAN CHANGE THE COLUMN TO INTEGRATE BY MODIFYNG THE NUMBER (i.e.c2 for the second column).
      sum=c3[0]+c3[N];
      for(int i=1;i<N-1;i+=2){
        sum += 4.*c3[i]+2.*c3[i+1];
      } 
    return sum*(dt/3.);
}

int main(){
  //Peak isolation of the choosed column with the specified file.
  cout<<simpson_integrator("FMO_proof.dat")<<endl;

  return 0;
}
