#include "interior_point_lp.h"
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>

// New function to read LP data from a text file in the given format
void readLPData(const std::string& filename, Eigen::MatrixXd &A, Eigen::VectorXd &b, Eigen::VectorXd &c, int &numVars, int &numConstraints) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    
    std::string line;
    // Read first line: numVars and numConstraints
    std::getline(file, line);
    std::istringstream iss_header(line);
    if (!(iss_header >> numVars >> numConstraints)) {
        throw std::runtime_error("Error reading dimensions.");
    }
    
    // Read objective coefficients (line 2)
    std::getline(file, line);
    std::istringstream iss_obj(line);
    std::vector<double> c_values;
    double val;
    while (iss_obj >> val) {
        c_values.push_back(val);
    }
    if (c_values.size() != static_cast<size_t>(numVars)) {
        throw std::runtime_error("Objective coefficients count mismatch.");
    }
    c = Eigen::VectorXd::Map(c_values.data(), c_values.size());
    
    // Read constraint matrix rows (next numConstraints lines)
    A = Eigen::MatrixXd(numConstraints, numVars);
    for (int i = 0; i < numConstraints; i++) {
        if (!std::getline(file, line)) {
            throw std::runtime_error("Not enough rows for constraint matrix.");
        }
        std::istringstream iss_row(line);
        for (int j = 0; j < numVars; j++) {
            if (!(iss_row >> val)) {
                throw std::runtime_error("Error reading matrix entry.");
            }
            A(i, j) = val;
        }
    }
    
    // Read right-hand side vector (last line)
    std::getline(file, line);
    std::istringstream iss_rhs(line);
    std::vector<double> b_values;
    while (iss_rhs >> val) {
        b_values.push_back(val);
    }
    if (b_values.size() != static_cast<size_t>(numConstraints)) {
        throw std::runtime_error("RHS vector size mismatch.");
    }
    b = Eigen::VectorXd::Map(b_values.data(), b_values.size());
    
    file.close();
}

int main(int argc, char** argv) {
    // Use command-line argument if provided, else default LP file path
    std::string lp_filename = "/home/syauqirp/gsoc25/data/feasible_lp105.txt";
    if (argc > 1) {
        lp_filename = argv[1];
    }
    
    int numVars = 0, numConstr = 0;
    Eigen::MatrixXd A;
    Eigen::VectorXd b, c;
    try {
        readLPData(lp_filename, A, b, c, numVars, numConstr);
    } catch (const std::exception& ex) {
        std::cerr << "Error reading LP file: " << ex.what() << std::endl;
        return 1;
    }
    std::cout << "Read LP problem from " << lp_filename << ": " << numVars << " variables, " << numConstr << " constraints." << std::endl;
    
    // Set solver parameters
    InteriorPointLP::Parameters params;
    params.tol = 1e-5;
    params.eta = 0.8;
    params.max_iter = 20000;
    params.regularization = 1e-6;
    params.use_scaling = false;  // Adjust as needed
    params.verbose = true;
    params.debug_level = 1;
    InteriorPointLP::setParameters(params);
    
    std::cout << "Solving LP problem with " << A.rows() << " constraints and " 
              << A.cols() << " variables..." << std::endl;
              
    InteriorPointLP::Result result = InteriorPointLP::solve(A, b, c);
    if (result.success) {
        std::cout << "Optimal solution found!" << std::endl;
    } else {
        std::cerr << "Solver terminated without finding an optimal solution." << std::endl;
        std::cerr << "Final primal infeasibility: " << result.primal_infeas << std::endl;
        std::cerr << "Final dual infeasibility: " << result.dual_infeas << std::endl;
        std::cerr << "Final gap: " << result.gap << std::endl;
    }
    
    return 0;
}