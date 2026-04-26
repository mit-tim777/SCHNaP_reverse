#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <filesystem>
#include <math.h>
#include "C:/devlibs/eigen-5.0.0/Eigen/Dense"

using namespace Eigen;
using namespace std;

typedef Matrix<double, 4, 4> Matrix4d;
typedef Matrix<double, 3, 3> Matrix3d;

double safe_double(const string& s) {
    try {
        return stod(s);
    } catch (...) {
        return -676767;
    }
}

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

map<string, vector<double>> load_csv_params(const filesystem::path& filepath, int start_col, int end_col) {
    map<string, vector<double>> params;
    ifstream file(filepath);
    string line;

    while (getline(file, line)) {
        stringstream ss(line);
        string cell;
        vector<string> row;
        
        while (getline(ss, cell, ',')) {
            row.push_back(cell);
        }

        if (row.size() > (size_t)end_col) {
            string key = row[0];
            vector<double> values;
            for (int i = start_col; i < end_col; ++i) {
                values.push_back(safe_double(row[i]));
            }
            params[key] = values;
        }
    }
    return params;
}

vector<vector<double>> read_pars_file(filesystem::path filename)
{
    vector<vector<double>> data; 
    ifstream f(filename);

    string t, tt;
    while(getline(f, t))
    {
        stringstream ss(t);
        vector<double> dat;

        while(getline(ss,tt, ' '))
        {
            double val = safe_double(tt);
            if( val != -676767 )
            {
                if(dat.size() >= 9) val = val * M_PI / 180.0;
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

    return data;
}

map<string, map<string, vector<double>>> load_equilibrium_params(const filesystem::path& csv_directory) {
    map<string, map<string, vector<double>>> equilibrium_params;

    // Load Step Params (Indices 1 to 6)
    equilibrium_params["step"] = load_csv_params(csv_directory / "coords_grooves_DNA_hexamers_table.csv", 1, 7);
    
    // Load BP Params (Indices 1 to 6)
    equilibrium_params["bp"] = load_csv_params(csv_directory / "coords_grooves_DNA_heptamers_table.csv", 1, 7);

    return equilibrium_params;
}

vector<string> get_all_possible_sequences(string sequence) {
    vector<int> placeholder_inds;
    for (size_t i = 0; i < sequence.length(); ++i) {
        if (sequence[i] == '-') placeholder_inds.push_back(i);
    }

    vector<string> possible_sequences = {sequence};
    vector<char> bases = {'A', 'T', 'C', 'G'};

    for (int idx : placeholder_inds) {
        vector<string> next_batch;
        for (const string& seq : possible_sequences) {
            for (char b : bases) {
                string new_seq = seq;
                new_seq[idx] = b;
                next_batch.push_back(new_seq);
            }
        }
        possible_sequences = next_batch;
    }
    return possible_sequences;
}

vector<double> get_equalibrium_params( const string& param_type,  const string& sequence, const map<string, map<string, vector<double>>>& equilibrium_params) 
{
    // 1. Get all permutations from the helper function we wrote previously
    vector<string> possible_sequences = get_all_possible_sequences(sequence);
    
    // 2. Initialize average parameters vector with 6 zeros
    vector<double> avg_params(6, 0.0);
    size_t num_seqs = possible_sequences.size();

    if (num_seqs == 0) return avg_params;

    // 3. Sum up parameters for all possible sequences
    for (const string& seq : possible_sequences) {
        const vector<double>& params = equilibrium_params.at(param_type).at(seq);
        for (int i = 0; i < 6; ++i) {
            avg_params[i] += params[i];
        }
    }

    // 4. Calculate the average
    for (int i = 0; i < 6; ++i) {
        avg_params[i] /= static_cast<double>(num_seqs);
    }

    return avg_params;
}

string get_hex_seq(const string& seq, int step_n) {
    string hex_seq = "";
    int seq_len = static_cast<int>(seq.length());

    // Iterate through the range [step_n - 2, step_n + 4)
    for (int i = step_n - 2; i < step_n + 4; ++i) {
        // Check if index is within the valid range of the sequence
        if (i >= 0 && i < seq_len) {
            hex_seq += seq[i];
        } else {
            hex_seq += '-';
        }
    }
    return hex_seq;
}

int main() {
    filesystem::path dir = "hexamers_csv/DNA/";
    auto params = load_equilibrium_params(dir);
    // auto data = read_pars_file("pars.txt");

    vector<Vector3d> end_pos_res;

    vector<string> seqs = get_all_possible_sequences("---------"); // 9*'-'   ->    1min 30s
    for(int i = 0; i < seqs.size(); i++)
    {
        Matrix4d T_i = Matrix4d::Identity();
        for(int j = 0; j < seqs[i].size(); j++)
        {
            vector<double> w = get_equalibrium_params("step", get_hex_seq(seqs[i],j), params);
            for( int k = 3; k < 6; k++) w[k] = w[k] * M_PI / 180.0;

            Matrix4d M = get_step_matrix(w);    
            T_i = T_i * M;        
        }
        Vector3d pos = T_i.block<3, 1>(0, 3);
        end_pos_res.push_back(pos.transpose());
    }

    cout << "done" << endl;
    return 0;
}