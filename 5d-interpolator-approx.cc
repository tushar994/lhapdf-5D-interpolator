
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Eigen/Eigen"

// type = 0 we use the slope avg with an approximation of second derivative using only four points
// type = 1, we use 6 points to find the second derivative
using namespace std;
int x_size = 91;
int q_size = 20;
double q_arr[20], x_arr[91]; 
double xf[20*91][3];

// This helps find the greatest x-value in the dataset below a given x value
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
// This finds the greatest x-value in the dataset below a given x value
int x_index(double x){
    return x_binary(0,x_size-1,x);
}
// This finds the greatest q-value in the dataset below a given q value
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
// This finds the greatest q-value in the dataset below a given q value
int q_index(double x){
    return q_binary(0,q_size-1,x);
}
// This returns index+add, unless index+add is an invalid value,
//  in which case it gives max or min depending on how indes+add overflowed
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
// This gives the xf value given the flavour, the x value's index in the 
// dataset and q value's index in the given dataset
double get_val(int x_index, int q_index, int flavour){
    int index = x_index*q_size + q_index;
    return xf[index][flavour];
}
// given a system of equations, this function solves it and gives the 
// outupt of the polynomial obtained
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
// this function, given any x value, and a q values index in the dataset, interpolates along x
double interpolate_x(double x, int q_i, int n, int flavour, int type){
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
        dy2 = ((dy2_2*(x_arr[x_i[2]] - x_arr[x_i[1]]))/(x_arr[x_i[3]] - x_arr[x_i[1]])) + ((dy2_1*(x_arr[x_i[3]] - x_arr[x_i[2]]))/(x_arr[x_i[3]] - x_arr[x_i[1]]));
    }

    // calculate second derivative
    double ddy1, ddy2;
    if(type==0){
        // if type = 0, then we get second derivative with only four points
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
        // cout<<"second\n";
        // cout<<"second derivatives are - "<<ddy1<<" "<<ddy2<<"\n";
    }
    else if(type==1){
        // we use six points to get the second derivative
        int xo_i = return_index(index1, x_size-1, 0, -2);
        int x6_i = return_index(index1, x_size-1, 0, 3);
        // cout<<"index - ";
        // cout<<xo_i<<" ";
        // for(int i=0;i<4;i++){
            // cout<<x_i[i]<<" ";
        // }
        // cout<<x6_i<<"\n";
        double dd_x1,dd_x3;
        if(xo_i==x_i[0] || x_i[0]==x_i[1]){
            dd_x1 = 0;
        }
        else{
            dd_x1 = 2*(((a[1] - a[0]) / (x_arr[x_i[1]] - x_arr[x_i[0]])) - ((a[0]-get_val(xo_i,q_i,0))/(x_arr[x_i[0]] - x_arr[xo_i])))/(x_arr[x_i[1]] - x_arr[xo_i]);
        }
        if(x_i[2]==x_i[3]){
            dd_x3 = 0;
        }
        else{
            dd_x3 = 2*(((a[3] - a[2]) / (x_arr[x_i[3]] - x_arr[x_i[2]])) - ((a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]])))/(x_arr[x_i[3]] - x_arr[x_i[1]]);
        }
        ddy1 = ((dd_x1*(x_arr[x_i[2]]-x_arr[x_i[1]]))+(dd_x3*(x_arr[x_i[1]]-x_arr[x_i[0]])))/(x_arr[x_i[2]]-x_arr[x_i[0]]);


        // for the second derivative as x3
        double dd_x2,dd_x4;
        if(x6_i==x_i[3] || x_i[3]==x_i[2]){
            dd_x4 = 0;
        }
        else{
            dd_x4 = (((get_val(x6_i,q_i,0)-a[3])/(x_arr[x6_i] - x_arr[x_i[3]])) - ((a[3] - a[2]) / (x_arr[x_i[3]] - x_arr[x_i[2]])))/((x_arr[x6_i] - x_arr[x_i[2]])/2);
        }
        if(x_i[0]==x_i[1]){
            dd_x2 = 0;
        }
        else{
            dd_x2 = 2*(((a[2] - a[1]) / (x_arr[x_i[2]] - x_arr[x_i[1]])) - ((a[1] - a[0]) / (x_arr[x_i[1]] - x_arr[x_i[0]])))/(x_arr[x_i[2]] - x_arr[x_i[0]]);
        }
        // cout<<dd_x2<<" "<<dd_x4<<" are the second derivatives things\n";
        ddy2 = ((dd_x2*(x_arr[x_i[3]]-x_arr[x_i[2]]))+(dd_x4*(x_arr[x_i[2]]-x_arr[x_i[1]])))/(x_arr[x_i[3]]-x_arr[x_i[1]]);
    }
    // cout<<"second erivative for x\n";
    // cout<<ddy1<<" "<<ddy2<<"\n";
    
    return solve_system(x,x_arr[x_i[1]],x_arr[x_i[2]],y1,y2,dy1,dy2,ddy1,ddy2);
}
// this function, given any q value, and the values that were gained from interpolation in the x direction, 
// gives the final interpolation value
double interpolate_q(double q, double* values , int n, int type){
    int q_i[4];
    int index1 = q_index(q);
    for(int i=0;i<4;i++){
        q_i[i] = return_index(index1, q_size-1, 0, i-1);
    }
    double a[4];
    // cout<<"values are (2):\n";
    for(int i=0;i<4;i++){
        a[i] = values[i];
        // cout<<a[i]<<" ";
    }
    if(type==1){
        for(int i=0;i<4;i++){
            a[i] = values[i+1];
        // cout<<a[i]<<" ";
        }   
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
        dy2 = ((dy2_2*(q_arr[q_i[2]] - q_arr[q_i[1]]))/(q_arr[q_i[3]] - q_arr[q_i[1]])) + ((dy2_1*(q_arr[q_i[3]] - q_arr[q_i[2]]))/(q_arr[q_i[3]] - q_arr[q_i[1]]));
    }

    // calculate second derivative
    double ddy1, ddy2;
    if(type==0){
        // if type = 0, then we get second derivative with only four points
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
    }
    else if(type==1){
        // we use siq points to get the second derivative
        int qo_i = return_index(index1, q_size-1, 0, -2);
        int q6_i = return_index(index1, q_size-1, 0, 3);
        double dd_q1,dd_q3;
        if(qo_i==q_i[0] || q_i[0]==q_i[1]){
            dd_q1 = 0;
        }
        else{
            dd_q1 = 2*(((a[1] - a[0]) / (q_arr[q_i[1]] - q_arr[q_i[0]])) - ((a[0]-values[0])/(q_arr[q_i[0]] - q_arr[qo_i])))/(q_arr[q_i[1]] - q_arr[qo_i]);
        }
        if(q_i[2]==q_i[3]){
            dd_q3 = 0;
        }
        else{
            dd_q3 = 2*(((a[3] - a[2]) / (q_arr[q_i[3]] - q_arr[q_i[2]])) - ((a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]])))/(q_arr[q_i[3]] - q_arr[q_i[1]]);
        }
        ddy1 = ((dd_q1*(q_arr[q_i[2]]-q_arr[q_i[1]]))+(dd_q3*(q_arr[q_i[1]]-q_arr[q_i[0]])))/(q_arr[q_i[2]]-q_arr[q_i[0]]);


        // for the second derivative as q3
        double dd_q2,dd_q4;
        if(q6_i==q_i[3] || q_i[3]==q_i[2]){
            dd_q4 = 0;
        }
        else{
            dd_q4 = -2*(((a[3] - a[2]) / (q_arr[q_i[3]] - q_arr[q_i[2]])) - ((values[5]-a[3])/(q_arr[q_i[0]] - q_arr[q6_i])))/(q_arr[q6_i] - q_arr[q_i[3]]);
        }
        if(q_i[0]==q_i[1]){
            dd_q2 = 0;
        }
        else{
            dd_q2 = 2*(((a[2] - a[1]) / (q_arr[q_i[2]] - q_arr[q_i[1]])) - ((a[1] - a[0]) / (q_arr[q_i[1]] - q_arr[q_i[0]])))/(q_arr[q_i[2]] - q_arr[q_i[0]]);
        }
        ddy2 = ((dd_q2*(q_arr[q_i[3]]-q_arr[q_i[2]]))+(dd_q4*(q_arr[q_i[2]]-q_arr[q_i[1]])))/(q_arr[q_i[3]]-q_arr[q_i[1]]);
    }

    
    return solve_system(q,q_arr[q_i[1]],q_arr[q_i[2]],y1,y2,dy1,dy2,ddy1,ddy2);
    
}

// --------------------------------------------- NEW SHIT
double six_point_inter(double x, double* values, double* x_val){
    double a[6];
    for(int i=0;i<6;i++){
        a[i] = values[i];
    }
    double y1 = a[2];
    double y2 = a[3];
    double dy1,dy2;
    
    double dy1_1 = (a[2] - a[1]) / (x_val[2] - x_val[1]);
    double dy1_2 = (a[3] - a[2]) / (x_val[3] - x_val[2]);
    dy1 = ((dy1_2*(x_val[2] - x_val[1]))/(x_val[3] - x_val[1])) + ((dy1_1*(x_val[3] - x_val[2]))/(x_val[3] - x_val[1]));

    double dy2_2 = (a[4] - a[3]) / (x_val[4] - x_val[3]);
    double dy2_1 = (a[3] - a[2]) / (x_val[3] - x_val[2]);
    dy2 = ((dy2_2*(x_val[3] - x_val[2]))/(x_val[4] - x_val[2])) + ((dy2_1*(x_val[4] - x_val[3]))/(x_val[4] - x_val[2]));

    // cout<<"index - ";
    // cout<<x_i[0]<<" ";
    // for(int i=0;i<4;i++){
        // cout<<x_i[i]<<" ";
    // }
    // cout<<x_i[5]<<"\n";
    double ddy1, ddy2;
    double dd_x1,dd_x3;
    dd_x1 = 2*(((a[2] - a[1]) / (x_val[2] - x_val[1])) - ((a[1]-a[0])/(x_val[1] - x_val[0])))/(x_val[2] - x_val[0]);
    dd_x3 = 2*(((a[4] - a[3]) / (x_val[4] - x_val[3])) - ((a[3] - a[2]) / (x_val[3] - x_val[2])))/(x_val[4] - x_val[2]);

    ddy1 = ((dd_x1*(x_val[3]-x_val[2]))+(dd_x3*(x_val[2]-x_val[1])))/(x_val[3]-x_val[1]);


    // for the second derivative as x3
    double dd_x2,dd_x4;
    dd_x4 = (((a[5]-a[4])/(x_val[5] - x_val[4])) - ((a[4] - a[3]) / (x_val[4] - x_val[3])))/((x_val[5] - x_val[3])/2);
    dd_x2 = 2*(((a[3] - a[2]) / (x_val[3] - x_val[2])) - ((a[3] - a[2]) / (x_val[2] - x_val[1])))/(x_val[3] - x_val[1]);
    
    // cout<<dd_x2<<" "<<dd_x4<<" are the second derivatives things\n";
    ddy2 = ((dd_x2*(x_val[4]-x_val[3]))+(dd_x4*(x_val[3]-x_val[2])))/(x_val[4]-x_val[2]);

    return solve_system(x,x_val[2],x_val[3],y1,y2,dy1,dy2,ddy1,ddy2);
}

double linear_inter(double x, double* values, double* x_val){
    double y1 = values[0];
    double y2 = values[1];
    double x1 = x_val[0];
    double x2 = x_val[1];
    return y1 + (y2-y1)*(x-x1)/(x2-x1);
}

double four_point_inter(double x, double* values, double* x_val){
    double a[4];
    for(int i=0;i<4;i++){
        a[i] = values[i];
    }

    double y1 = a[1];
    double y2 = a[2];
    double dy1,dy2;
    double dy1_1 = (a[1] - a[0]) / (x_val[1] - x_val[0]);
    double dy1_2 = (a[2] - a[1]) / (x_val[2] - x_val[1]);
    dy1 = ((dy1_2*(x_val[1] - x_val[0]))/(x_val[2] - x_val[0])) + ((dy1_1*(x_val[2] - x_val[1]))/(x_val[2] - x_val[0]));

    double dy2_2 = (a[3] - a[2]) / (x_val[3] - x_val[2]);
    double dy2_1 = (a[2] - a[1]) / (x_val[2] - x_val[1]);
    dy2 = ((dy2_2*(x_val[2] - x_val[1]))/(x_val[3] - x_val[1])) + ((dy2_1*(x_val[3] - x_val[2]))/(x_val[3] - x_val[1]));


    double ddy1, ddy2;

    // if type = 0, then we get second derivative with only four points
    ddy1 = 2*(((a[2] - a[1]) / (x_val[2] - x_val[1])) - ((a[1] - a[0]) / (x_val[1] - x_val[0])))/(x_val[2] - x_val[0]);
    ddy2 = 2*(((a[3] - a[2]) / (x_val[3] - x_val[2])) - ((a[2] - a[1]) / (x_val[2] - x_val[1])))/(x_val[3] - x_val[1]);

    // cout<<"second\n";
    // cout<<"second derivatives are - "<<ddy1<<" "<<ddy2<<"\n";

    return solve_system(x,x_val[1],x_val[2],y1,y2,dy1,dy2,ddy1,ddy2);

}



// ---------------------------------------------
// this is the main function from which all interpolation is called
double interpolate(double x, double q, int flavour, int type){
    int x_i = x_index(x);
    int q_i = q_index(q);
    int n = 4;
    int x_a[n],q_a[n];
    int x_a2[6],q_a2[6];
    if(type==0){
        for(int i=0;i<n;i++){
            x_a[i] = return_index(x_i, x_size-1, 0, i-1);
            q_a[i] = return_index(q_i, q_size-1, 0, i-1);
        }
    }
    else{
        for(int i=0;i<6;i++){
            x_a2[i] = return_index(x_i, x_size-1, 0, i-2);
            q_a2[i] = return_index(q_i, q_size-1, 0, i-2);
        }
    }
    // for(int i=0;i<6;i++){
    //     cout<<x_a2[i]<<" ";
    // }
    // cout<<"\n";

    // now we have all four coordinates for x, q
    double values[4];
    double values2[6];
    // we interpolate along x four times for the relevant values of q and store the obtained values
    // cout<<"values are\n";
    if(type==0){
        if(x_a[0]==x_a[1]||x_a[2]==x_a[3]){
            // degrenate to linear
            // cout<<"here\n";
            double x_vals[2];
            x_vals[0] = x_arr[x_a[1]];
            x_vals[1] = x_arr[x_a[2]];
            for(int i=0;i<4;i++){
                double vals[2];
                vals[0] = get_val(x_a[1],q_a[i],0);
                vals[1] = get_val(x_a[2],q_a[i],0);
                values[i] = linear_inter(x,vals,x_vals);
            }
        }
        else{
            double x_vals[4];
            for(int i=0;i<4;i++){
                x_vals[i]= x_arr[x_a[i]];
            }
            for(int i=0;i<4;i++){
                double vals[4];
                for(int j=0;j<4;j++){
                    vals[j] = get_val(x_a[j],q_a[i],0);
                }
                values[i] = four_point_inter(x,vals,x_vals);
                // cout<<values[i]<<" ";
            }
        }
    }
    else if(type==1){
        if(x_a2[0]==x_a2[1]||x_a2[4]==x_a2[5]){
            // degrenate to linear
            // cout<<"here\n";
            double x_vals[2];
            x_vals[0] = x_arr[x_a2[2]];
            x_vals[1] = x_arr[x_a2[3]];
            for(int i=0;i<6;i++){
                double vals[2];
                vals[0] = get_val(x_a2[2],q_a2[i],0);
                vals[1] = get_val(x_a2[3],q_a2[i],0);
                values2[i] = linear_inter(x,vals,x_vals);
            }
        }
        // normal 6 point interpolation
        else{
            double x_vals[6];
            for(int i=0;i<6;i++){
                x_vals[i]= x_arr[x_a2[i]];
            }
            for(int i=0;i<6;i++){
                double vals[6];
                for(int j=0;j<6;j++){
                    vals[j] = get_val(x_a2[j],q_a2[i],0);
                }
                values2[i] = six_point_inter(x,vals,x_vals);
                // cout<<values[i]<<" ";
            }
        }
    }
    // cout<<"\n";
    // Using the four values obtained we interpolate along q

    if(type==0){
        if(q_a[0]==q_a[1]||q_a[2]==q_a[3]){
            // degrenate to linear
            // cout<<"here\n";
            double q_vals[2];
            q_vals[0] = q_arr[q_a[1]];
            q_vals[1] = q_arr[q_a[2]];
            double vals[2];
            vals[0] = values[1];
            vals[1] = values[2];
            return linear_inter(q,vals,q_vals);
        }
        else{
            double q_vals[4];
            for(int i=0;i<4;i++){
                q_vals[i]= q_arr[q_a[i]];
            }
            return four_point_inter(q,values,q_vals);
        }
    }
    else if(type==1){
        // return interpolate_q(q,values2, 6,type);
        if(q_a2[0]==q_a2[1]||q_a2[4]==q_a2[5]){
            // degrenate to linear
            double q_vals[2];
            q_vals[0] = q_arr[q_a2[2]];
            q_vals[1] = q_arr[q_a2[3]];
            double vals[2];
            vals[0] = values2[2];
            vals[1] = values2[3];
            return linear_inter(q,vals,q_vals);
        }
        else{
            double q_vals[6];
            for(int i=0;i<6;i++){
                q_vals[i]= q_arr[q_a2[i]];
            }
            return six_point_inter(q,values2,q_vals);
        }
    }

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
    ofstream MyFile("./custom_func_data/4-points.csv");

    double x_val = -20;
    double q_val = 0.5;
    // cout<<interpolate((-0.4),(1),0,1)<<"\n";
    MyFile<<"lnx,lnq,xf\n";
    while(q_val<10){
        x_val = -18;
        while(x_val<0){
            MyFile<<x_val<<","<<q_val<<","<<interpolate((x_val),(q_val),0,0)<<"\n";
            // cout<<x_val<<"\r";
            x_val+=0.01;
        }
        q_val+=0.05;
        cout<<q_val<<"\n";
    }
    MyFile.close();


}