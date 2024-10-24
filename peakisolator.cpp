//Used packages (No need for extra installations or makefiles).
#include <iomanip>
#include <fstream>
#include <iostream>
using namespace std;
//Function for counting file lenght.
int count_lines(string filename){
  int count=0;
  string line;
  ifstream file(filename);
    while (getline(file,line)){
     count++;
    }
return count;
}
//Peak isolating function (column from 0 to 3).
void peak_isolator(string filename, int column){
  //Opening of input file.
  ifstream in(filename);
  //File lenght (N).
  int N= count_lines(filename);
  //Definition of time and columns arrays.
  double t[N], c0[N], c1[N],c2[N],c3[N] ;
  //Reading and assignation of file columns to arrays.
  for(int i=1;i<N;i++){
    in >> t[i] >> c0[i] >> c1[i] >> c2[i] >> c3[i];
  }
  //Opening of output file.
  ofstream coquito;
  coquito.open("maxmin.dat", ios::out | ios::trunc);
  //Cases of peak_isolator depending of the column parameter.
  switch(column){
    //column 0.
    case 0:

      for(int i=1;i<N;i++){

        if(c0[i-1]<c0[i] && c0[i+1]<c0[i]){
          coquito<<i<<" "<< c0[i] <<endl;}

        if (c0[i-1]>c0[i] && c0[i+1]>c0[i]){
          coquito<<i<<" "<< c0[i]<< endl;}
      }
      break; 
    //column 1.
    case 1:

      for(int i=1;i<N;i++){

        if(c1[i-1]<c1[i] && c1[i+1]<c1[i]){
          coquito<<i<<" "<< c1[i] << endl;}

        if (c1[i-1]>c1[i] && c1[i+1]>c1[i]){
          coquito<<i<<" "<< c1[i] << endl;}
      }
      break;
    //column 2.
    case 2:

      for(int i=1;i<N;i++){

        if(c2[i-1]<c2[i] && c2[i+1]<c2[i]){
          coquito<<i<<" "<< c2[i] << endl;}

        if (c2[i-1]>c2[i] && c2[i+1]>c2[i]){
          coquito<<i<<" "<< c2[i] << endl;}
      }
      break;
    //column 3.
    case 3:

      for(int i=1;i<N;i++){

        if(c3[i-1]<c3[i] && c3[i+1]<c3[i]){
          coquito<<i<<" "<< c3[i] << endl;}

        if (c3[i-1]>c3[i] && c3[i+1]>c3[i]){
          coquito<<i<<" "<< c3[i] << endl;}
      }
      break;
  }
}
int main(){

  //Definition of the column with the data that is wanted to peak isolate.
  int column=1;
  //Peak isolation of the choosed column with the specified file.
  peak_isolator("FMO_proof.dat", column);

  return 0;
}
