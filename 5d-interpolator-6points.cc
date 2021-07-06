
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Eigen/Eigen"

// we use the slope avg with an approximation of second derivative using only four points

using namespace std;
int x_size = 91;
int q_size = 20;
double q_arr[20], x_arr[91]; 
double xf[20*91][3];


int x_binary(int min, int max, double val){
    while(max-min>1){
        int mid  = (min+max)/2; 
        if(x_arr[mid]>val){
            max = mid;
        }
        else{
            min = mid;
        }
    }
    return min;
}

int x_index(double x){
    return x_binary(0,x_size-1,x);
}

int q_binary(int min, int max, double val){
    while(max-min>1){
        int mid  = (min+max)/2; 
        if(q_arr[mid]>val){
            max = mid;
        }
        else{
            min = mid;
        }
    }
    return min;
}

int q_index(double x){
    return q_binary(0,q_size-1,x);
}

int return_index(int index, int max, int min, int add){
    if(index+add < min){
        return min;
    }
    else if(index+add > max){
        return max;
    }
    else{
        return index+add;
    }

}   

double get_val(int x_index, int q_index, int flavour){
    int index = x_index*q_size + q_index;
    return xf[index][flavour];
}

double solve_system(double x,double x1, double x2, double y1, double y2, double dy1, double dy2, double ddy1, double ddy2){
    Eigen::Matrix <double, 1,7> equations[6];
    
    double pow_x[6];
    pow_x[0] = 1;
    for(int i=1;i<6;i++){
        pow_x[i] = pow_x[i-1]*x;
    }
    
    double pow_x1[6];
    pow_x1[0] = 1;
    for(int i=1;i<6;i++){
        pow_x1[i] = pow_x1[i-1]*x1;
    }
    double pow_x2[6];
    pow_x2[0] = 1;
    for(int i=1;i<6;i++){
        pow_x2[i] = pow_x2[i-1]*x2;
    }

    // normal function equations
    for(int i=0;i<6;i++){
        equations[0](0,i) = pow_x1[5-i];
    }
    equations[0](0,6) = y1;
    for(int i=0;i<6;i++){
        equations[1](0,i) = pow_x2[5-i];
    }
    equations[1](0,6) = y2;


    // first derivative equations
    for(int i=0;i<6;i++){
        if(5-i-1 > 0){
            equations[2](0,i) = (5-i)*pow_x1[5-i-1];
        }
        else{
            equations[2](0,i) = 0;
        }
    }
    equations[2](0,6) = dy1;

    for(int i=0;i<6;i++){
        if(5-i-1 > 0){
            equations[3](0,i) = (5-i)*pow_x2[5-i-1];
        }
        else{
            equations[3](0,i) = 0;
        }
    }
    equations[3](0,6) = dy2;

    // second derivative equations
    for(int i=0;i<6;i++){
        if(5-i-2 > 0){
            equations[4](0,i) = (5-i)*(5-i-1)*pow_x1[5-i-2];
        }
        else{
            equations[4](0,i) = 0;
        }
    }
    equations[4](0,6) = ddy1;

    for(int i=0;i<6;i++){
        if(5-i-2 > 0){
            equations[5](0,i) = (5-i)*(5-i-1)*pow_x2[5-i-2];
        }
        else{
            equations[5](0,i) = 0;
        }
    }
    equations[5](0,6) = ddy2;

    for(int i=0;i<6;i++){
        for(int j=i+1;j<6;j++){
            equations[j] = equations[j] - (equations[i]*( equations[j](0,i)/equations[i](0,i) ));
        }
    }

    for(int i=5;i>=0;i--){
        for(int j=i-1;j>=0;j--){
            equations[j] = equations[j] - (equations[i]*( equations[j](0,i)/equations[i](0,i) ));
        }
    }

    double  coeff[6];
    for(int i=0;i<6;i++){
        coeff[i] = equations[i](0,6)/equations[i](0,i);
    }

    double ans = 0;

    for(int i=0;i<6;i++){
        ans+= coeff[i]*pow_x[5-i];
    }
    return ans;

}

double interpolate_x(double x, int q_i, int n, int flavour){
    int x_i[n];
    int index1 = x_index(x);
    for(int i=0;i<n;i++){
        x_i[i] = return_index(index1, x_size-1, 0, i-1);
    }
    // x_i contains the index at which you'll find the x's you want
    double a[n];
    for(int i=0;i<n;i++){
        a[i] = get_val(x_i[i],q_i,flavour);
    }
    // a contains the four values
    double y1 = a[1];
    double y2 = a[2];
    double dy1,dy2;
    if(x_i[0]==x_i[1]){
        double dy1_2 = (a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]]);
        dy1 = dy1_2;
    }
    else{
        double dy1_1 = (a[1] - a[0]) / (x_arr[x_i[1]] - x_arr[x_i[0]]);
        double dy1_2 = (a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]]);
        dy1 = ((dy1_2*(x_arr[x_i[1]] - x_arr[x_i[0]]))/(x_arr[x_i[2]] - x_arr[x_i[0]])) + ((dy1_1*(x_arr[x_i[2]] - x_arr[x_i[1]]))/(x_arr[x_i[2]] - x_arr[x_i[0]]));
    }

    if(x_i[2]==x_i[3]){
        double dy2_1 = (a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]]);
        dy2 = dy2_1;
    }
    else{
        double dy2_2 = (a[3] - a[2]) / (x_arr[x_i[3]] - x_arr[x_i[2]]);
        double dy2_1 = (a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]]);
        dy1 = ((dy2_2*(x_arr[x_i[2]] - x_arr[x_i[1]]))/(x_arr[x_i[3]] - x_arr[x_i[1]])) + ((dy2_1*(x_arr[x_i[3]] - x_arr[x_i[2]]))/(x_arr[x_i[3]] - x_arr[x_i[1]]));
    }

    double ddy1, ddy2;

    if(x_i[0]==x_i[1]){
        ddy1 = 0;
    }
    else{
        ddy1 = 2*(((a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]])) - ((a[1] - a[0]) / (x_arr[x_i[1]] - x_arr[x_i[0]])))/(x_arr[x_i[2]] - x_arr[x_i[0]]);
    }

    if(x_i[2]==x_i[3]){
        ddy2 = 0;
    }
    else{
        ddy2 = 2*(((a[3] - a[2]) / (x_arr[x_i[3]] - x_arr[x_i[2]])) - ((a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]])))/(x_arr[x_i[3]] - x_arr[x_i[1]]);
    }

    
    return solve_system(x,x_arr[x_i[1]],x_arr[x_i[2]],y1,y2,dy1,dy2,ddy1,ddy2);
}

double interpolate_q(double q, double* values , int n){
    int q_i[n];
    int index1 = q_index(q);
    for(int i=0;i<n;i++){
        q_i[i] = return_index(index1, q_size-1, 0, i-1);
    }
    double a[n];
    // cout<<"values are (2):\n";
    for(int i=0;i<n;i++){
        a[i] = values[i];
        // cout<<a[i]<<" ";
    }
    // a contains the four values
    double y1 = a[1];
    double y2 = a[2];
    double dy1,dy2;
    if(q_i[0]==q_i[1]){
        double dy1_2 = (a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]]);
        dy1 = dy1_2;
    }
    else{
        double dy1_1 = (a[1] - a[0]) / (q_arr[q_i[1]] - q_arr[q_i[0]]);
        double dy1_2 = (a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]]);
        dy1 = ((dy1_2*(q_arr[q_i[1]] - q_arr[q_i[0]]))/(q_arr[q_i[2]] - q_arr[q_i[0]])) + ((dy1_1*(q_arr[q_i[2]] - q_arr[q_i[1]]))/(q_arr[q_i[2]] - q_arr[q_i[0]]));
    }

    if(q_i[2]==q_i[3]){
        double dy2_1 = (a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]]);
        dy2 = dy2_1;
    }
    else{
        double dy2_2 = (a[3] - a[2]) / (q_arr[q_i[3]] - q_arr[q_i[2]]);
        double dy2_1 = (a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]]);
        dy1 = ((dy2_2*(q_arr[q_i[2]] - q_arr[q_i[1]]))/(q_arr[q_i[3]] - q_arr[q_i[1]])) + ((dy2_1*(q_arr[q_i[3]] - q_arr[q_i[2]]))/(q_arr[q_i[3]] - q_arr[q_i[1]]));
    }

    double ddy1, ddy2;

    if(q_i[0]==q_i[1]){
        ddy1 = 0;
    }
    else{
        ddy1 = 2*(((a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]])) - ((a[1] - a[0]) / (q_arr[q_i[1]] - q_arr[q_i[0]])))/(q_arr[q_i[2]] - q_arr[q_i[0]]);
    }

    if(q_i[2]==q_i[3]){
        ddy2 = 0;
    }
    else{
        ddy2 = 2*(((a[3] - a[2]) / (q_arr[q_i[3]] - q_arr[q_i[2]])) - ((a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]])))/(q_arr[q_i[3]] - q_arr[q_i[1]]);
    }

    
    return solve_system(q,q_arr[q_i[1]],q_arr[q_i[2]],y1,y2,dy1,dy2,ddy1,ddy2);
    
}

double interpolate(double x, double q, int flavour){
    int x_i = x_index(x);
    int q_i = q_index(q);
    int n = 4;
    int x_a[n],q_a[n];
    for(int i=0;i<n;i++){
        x_a[i] = return_index(x_i, x_size-1, 0, i-1);
        q_a[i] = return_index(q_i, q_size-1, 0, i-1);
    }

    // now we have all four coordinates for x, q
    double values[4];
    // we interpolate along x four times for the relevant values of q and store the obtained values
    // cout<<"values are\n";
    for(int i=0;i<4;i++){
        values[i] = interpolate_x(x,q_a[i],n,flavour);
        // cout<<values[i]<<" ";
    }
    // cout<<"\n";
    // Using the four values obtained we interpolate along q
    return interpolate_q(q,values, n);

}

int main(){

    fstream fin;

    // change these according to dataset
    // Open an existing file
    fin.open("./custom_func_data/my_info.csv", ios::in);
    vector<string> row;
    string line, word, temp;
    int line1 = 0;
    // getting xf
    while (fin >> temp) {

        row.clear();
        // cout<<temp<<"\n";
        // read an entire row and
        // store it in a string variable 'line'
        // getline(fin, line);
        // cout<<line<<"\n";
        // used for breaking words
        stringstream s(temp);

        // read every column data of a row and
        // store it in a string variable, 'word'
        int flavour  =0 ;
        while (getline(s, word, ',')) {

            // add all the column data
            // of a row to a vector
            // cout<<word<<"\n";
            // if(word=="---"){
            //     break;
            // }
            xf[line1][flavour] = stod(word);
            flavour++;
        }
        // if(word=='---'){
        //     break;
        // }
        line1++;
    }
    fstream fin1;
    fin1.open("./custom_func_data/x.csv", ios::in);
    line1 = 0;
    while (fin1 >> temp) {

        row.clear();
        // cout<<temp<<"\n";
        // read an entire row and
        // store it in a string variable 'line'
        // getline(fin, line);
        // cout<<line<<"\n";
        // used for breaking words
        stringstream s(temp);

        // read every column data of a row and
        // store it in a string variable, 'word'
        int flavour  =0 ;
        while (getline(s, word, ',')) {

            // add all the column data
            // of a row to a vector
            x_arr[line1] = stod(word);
            line1++;

        }
    }
    // for(int i=0;i<190;i++){
    //     cout<<x[i]<<" ";
    // }
    // cout<<"\n";
    fstream fin2;
    fin2.open("./custom_func_data/q.csv", ios::in);
    line1 = 0;
    while (fin2 >> temp) {

        row.clear();
        // cout<<temp<<"\n";
        // read an entire row and
        // store it in a string variable 'line'
        // getline(fin, line);
        // cout<<line<<"\n";
        // used for breaking words
        stringstream s(temp);

        // read every column data of a row and
        // store it in a string variable, 'word'
        int flavour  =0 ;
        while (getline(s, word, ',')) {

            // add all the column data
            // of a row to a vector
            q_arr[line1] = stod(word);
            line1++;

        }
    }
    //<<----------------------------------_-----------------_-----------------_-----------------_-----------------_-----------------_-----------------_-----------------_>
    // for(int i=0;i<4560;i++){
    //     cout<<xf[i][0]<<" ";
    // }
    // cout<<"\n";
    // cout<<get_val(1,1,0)<<" is the val\n";
    ofstream MyFile("./dense_grids/4-points-approx.csv");

    double x_val = -20;
    double q_val = 0.5;
    // cout<<interpolate(exp(-18),exp(1),0)<<"\n";
    MyFile<<"lnx,lnq,xf\n";
    while(q_val<10){
        x_val = -18;
        while(x_val<0){
            MyFile<<x_val<<","<<q_val<<","<<interpolate(exp(x_val),exp(q_val),0)<<"\n";
            // cout<<x_val<<"\r";
            x_val+=0.001;
        }
        q_val+=0.5;
        cout<<q_val<<"\n";
    }
    MyFile.close();


}