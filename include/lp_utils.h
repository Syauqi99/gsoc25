#ifndef LP_UTILS_H
#define LP_UTILS_H

#include <Eigen/Dense>
#include <iostream>

namespace LPUtils {

/**
 * Structure to store scaling information for linear programming problems
 */
struct ScalingInfo {
    bool is_scaled = false;
    Eigen::VectorXd row_scaling;
    Eigen::VectorXd col_scaling;
};

/**
 * Scale a linear programming problem for numerical stability
 * @param A The constraint matrix
 * @param b The right-hand side vector
 * @param c The objective coefficient vector
 * @return Scaling information for rescaling the solution
 */
ScalingInfo scaleLP(Eigen::MatrixXd& A, Eigen::VectorXd& b, Eigen::VectorXd& c);

/**
 * Rescale the solution back to the original problem
 * @param x Primal variables
 * @param lambda Dual variables
 * @param s Slack variables
 * @param scaling Scaling information from scaleLP
 */
void rescaleSolution(Eigen::VectorXd& x, Eigen::VectorXd& lambda, Eigen::VectorXd& s, 
                      const ScalingInfo& scaling);

/**
 * Check if a vector contains NaN or Infinity values
 * @param vec The vector to check
 * @return True if the vector contains NaN or Infinity values
 */
bool containsNanOrInf(const Eigen::VectorXd& vec);

/**
 * Print diagnostic information about solution vectors
 * @param x Primal variables
 * @param lambda Dual variables
 * @param s Slack variables
 * @param label Label to identify the diagnostic output
 */
void printDiagnostics(const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, 
                      const Eigen::VectorXd& s, const char* label);

} // namespace LPUtils

#endif // LP_UTILS_H
