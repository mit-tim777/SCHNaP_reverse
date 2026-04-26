#include <iostream>
// #include <bits/stdc++.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "C:/devlibs/eigen-5.0.0/Eigen/Dense"

using namespace Eigen;
using namespace std;

typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<double, 3, 3> Matrix3d;
// typedef Vector3d Vector3d;

// Rodrigues' formula
Matrix3d rodrigues_matrix(Vector3d axis, double angle) {
    double norm = axis.norm();
    if (norm < 1e-9) return Matrix3d::Identity();
    axis /= norm;

    Matrix3d K;
    K << 0, -axis(2), axis(1),
         axis(2), 0, -axis(0),
         -axis(1), axis(0), 0;

    return Matrix3d::Identity() + sin(angle) * K + (1.0 - cos(angle)) * (K * K);
}

Matrix3d get_rotation_matrix(const vector<double>& w) {
    double L = sqrt(pow(w[3], 2) + pow(w[4], 2));
    double o = atan2(w[3], w[4]);
    double theta_axis = o - w[5] / 2.0;
    
    Vector3d u(sin(theta_axis), cos(theta_axis), 0.0);
    
    Matrix3d R_tilt = rodrigues_matrix(u, L);
    Matrix3d R_twist = rodrigues_matrix(Vector3d(0, 0, 1), w[5]);
    return R_tilt * R_twist;
}

Matrix3d get_half_rotation_matrix(const vector<double>& w) {
    double L = sqrt(pow(w[3], 2) + pow(w[4], 2));
    double o = atan2(w[3], w[4]);
    double theta_axis = o - w[5] / 2.0;
    
    Vector3d u(sin(theta_axis), cos(theta_axis), 0.0);
    
    Matrix3d R_half_tilt = rodrigues_matrix(u, L / 2.0);
    Matrix3d R_half_twist = rodrigues_matrix(Vector3d(0, 0, 1), w[5] / 2.0);
    return R_half_tilt * R_half_twist;
}

Matrix4d get_step_matrix(const vector<double>& w) {
    Matrix4d M = Matrix4d::Identity();
    M.block<3, 3>(0, 0) = get_rotation_matrix(w);
    M.block<3, 1>(0, 3) = get_half_rotation_matrix(w) * Vector3d(w[0], w[1], w[2]);
    return M;
}

double saveDouble(string s)
{   
    double d;
    try{
        return stod(s);
    }
    catch (const invalid_argument& ia )
    {
        return -676767;
    }
}

int main() {
    // Example usage for processing data
    vector<vector<double>> data; 
    // data = 

    ifstream f("pars.txt");

    string t, tt;
    while(getline(f, t))
    {
        stringstream ss(t);
        vector<double> dat;

        while(getline(ss,tt, ' '))
        {
            double val = saveDouble(tt);
            if( val != -676767 )
            {
                dat.push_back(val);
            }
        }
        if(dat.size() == 12)
        {
            data.push_back(dat);
        }
    };
    for( int i = 0; i < data.size(); i++)
    {
        data[i].erase(data[i].begin(), data[i].begin() + 6);
    }
    data.erase(data.end());


    Matrix4d T_i = Matrix4d::Identity();

    for (const auto& w : data) {
        Matrix4d M = get_step_matrix(w);
        T_i = T_i * M;
        
        // Access T_i components
        Vector3d pos = T_i.block<3, 1>(0, 3);
        cout << "Position: " << pos.transpose() << endl;
    }

    return 0;
}