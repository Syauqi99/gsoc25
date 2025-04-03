#include "lp_utils.h"
#include <limits>
#include <algorithm>
#include <iostream>

namespace LPUtils {

LPUtils::ScalingInfo scaleLP(Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& c) {
    ScalingInfo scaling;
    const int m = A.rows();
    const int n = A.cols();
    
    // Skip scaling for small problems
    if (n < 50 && m < 50) {
        scaling.row_scaling = Eigen::VectorXd::Ones(m);
        scaling.col_scaling = Eigen::VectorXd::Ones(n);
        return scaling;
    }
    
    scaling.is_scaled = true;
    scaling.row_scaling = Eigen::VectorXd::Ones(m);
    scaling.col_scaling = Eigen::VectorXd::Ones(n);
    
    // Set limits for scaling factors to avoid extreme scaling
    const double MAX_SCALING = 1e6;
    const double MIN_SCALING = 1e-6;
    
    // Iterative scaling (similar to matrix equilibration)
    for (int iter = 0; iter < 5; iter++) {
        // Scale rows
        for (int i = 0; i < m; i++) {
            double row_max = A.row(i).cwiseAbs().maxCoeff();
            if (row_max > 0) {
                double scale = std::min(std::max(1.0 / row_max, MIN_SCALING), MAX_SCALING);
                A.row(i) *= scale;
                b(i) *= scale;
                scaling.row_scaling(i) *= scale;
            }
        }
        
        // Scale columns
        for (int j = 0; j < n; j++) {
            double col_max = A.col(j).cwiseAbs().maxCoeff();
            if (col_max > 0) {
                double scale = std::min(std::max(1.0 / col_max, MIN_SCALING), MAX_SCALING);
                A.col(j) *= scale;
                c(j) *= scale;
                scaling.col_scaling(j) *= scale;
            }
        }
    }
    
    return scaling;
}

void rescaleSolution(Eigen::VectorXd& x, Eigen::VectorXd& lambda, Eigen::VectorXd& s, 
                      const ScalingInfo& scaling) {
    if (!scaling.is_scaled) {
        return;
    }
    
    // Rescale primal variables
    for (int i = 0; i < x.size(); i++) {
        x(i) *= scaling.col_scaling(i);
    }
    
    // Rescale dual variables
    for (int i = 0; i < lambda.size(); i++) {
        lambda(i) /= scaling.row_scaling(i);
    }
    
    // Rescale slack variables
    for (int i = 0; i < s.size(); i++) {
        s(i) /= scaling.col_scaling(i);
    }
}

bool containsNanOrInf(const Eigen::VectorXd& vec) {
    for (int i = 0; i < vec.size(); i++) {
        if (std::isnan(vec(i)) || std::isinf(vec(i))) {
            return true;
        }
    }
    return false;
}

void printDiagnostics(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, 
                      const Eigen::VectorXd& s, const char* label) {
    std::cout << "=== Diagnostics: " << label << " ===" << std::endl;
    std::cout << "x range: [" << x.minCoeff() << ", " << x.maxCoeff() << "]" << std::endl;
    std::cout << "lambda range: [" << lambda.minCoeff() << ", " << lambda.maxCoeff() << "]" << std::endl;
    std::cout << "s range: [" << s.minCoeff() << ", " << s.maxCoeff() << "]" << std::endl;
    std::cout << "Any NaN/Inf in x: " << (containsNanOrInf(x) ? "YES" : "no") << std::endl;
    std::cout << "Any NaN/Inf in lambda: " << (containsNanOrInf(lambda) ? "YES" : "no") << std::endl;
    std::cout << "Any NaN/Inf in s: " << (containsNanOrInf(s) ? "YES" : "no") << std::endl;
}

} // end namespace LPUtils
