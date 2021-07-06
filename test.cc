#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;


int main(){
    Matrix<double, 1, 7> arr[10];
    Matrix<double,1,7> return_one;
    for(int i=0;i<10;i++){
        for(int j=0;j<7;j++){
            arr[i](0,j) = j;
        }
    }
    add(7,arr,&return_one);
    cout<<return_one<<"\n";

}