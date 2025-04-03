#ifndef INTERIOR_POINT_LP_H
#define INTERIOR_POINT_LP_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "lp_utils.h"

class InteriorPointLP {
public:
    // Structure to hold the result of the solve function
    struct Result {
        bool success;             // Whether the solver succeeded
        Eigen::VectorXd x;        // Optimal solution
        double optimal_value;     // Optimal objective value
        double primal_infeas;     // Primal infeasibility
        double dual_infeas;       // Dual infeasibility
        double gap;               // Complementarity gap
    };

    // Main solver function
    static Result solve(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& c);
    
    // Algorithm parameters
    struct Parameters {
        double tol = 1e-6;          // Tolerance for convergence
        int max_iter = 2000;        // Maximum iterations
        double eta = 0.9;           // Step length scaling factor 
        double regularization = 1e-8; // Regularization parameter for matrix factorization
        bool use_scaling = true;    // Whether to scale the problem
        bool verbose = false;       // Print detailed progress information
        int debug_level = 0;        // Debug level: 0=none
    };
    
    // Set algorithm parameters
    static void setParameters(const Parameters& params);

private:
    static Parameters params;
    
    // Calculate initial point
    static void computeInitialPoint(
        const Eigen::MatrixXd& A, 
        const Eigen::VectorXd& b, 
        const Eigen::VectorXd& c,
        Eigen::VectorXd& x, 
        Eigen::VectorXd& lambda, 
        Eigen::VectorXd& s);
    
    // Compute affine scaling direction (predictor step)
    static void computeAffineDirection(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& lambda, 
        const Eigen::VectorXd& s,
        Eigen::VectorXd& dx_aff, 
        Eigen::VectorXd& dlambda_aff, 
        Eigen::VectorXd& ds_aff,
        const Eigen::VectorXd& rc, 
        const Eigen::VectorXd& rb);
    
    // Compute step lengths
    static void computeStepLengths(
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& s,
        const Eigen::VectorXd& dx, 
        const Eigen::VectorXd& ds,
        double& alpha_pri, 
        double& alpha_dual);
    
    // Compute centering parameter
    static double computeCenteringParameter(
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& s,
        const Eigen::VectorXd& dx_aff, 
        const Eigen::VectorXd& ds_aff,
        double alpha_pri_aff, 
        double alpha_dual_aff,
        double mu);
    
    // Compute combined direction (corrector step)
    static void computeCombinedDirection(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& lambda, 
        const Eigen::VectorXd& s,
        const Eigen::VectorXd& dx_aff, 
        const Eigen::VectorXd& ds_aff,
        Eigen::VectorXd& dx, 
        Eigen::VectorXd& dlambda, 
        Eigen::VectorXd& ds,
        const Eigen::VectorXd& rc, 
        const Eigen::VectorXd& rb,
        double sigma, 
        double mu);
    
    // Solve the linear system with improved numerical stability
    static bool solveLinearSystem(
        const Eigen::MatrixXd& A,
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& s,
        const Eigen::VectorXd& rhs1, 
        const Eigen::VectorXd& rhs2, 
        const Eigen::VectorXd& rhs3,
        Eigen::VectorXd& dx, 
        Eigen::VectorXd& dlambda, 
        Eigen::VectorXd& ds);
    
    // Check convergence criteria
    static bool checkConvergence(
        const Eigen::VectorXd& x, 
        const Eigen::VectorXd& lambda, 
        const Eigen::VectorXd& s,
        const Eigen::VectorXd& rc, 
        const Eigen::VectorXd& rb,
        double mu,
        const Eigen::VectorXd& b_orig, 
        const Eigen::VectorXd& c_orig);
};

#endif // INTERIOR_POINT_LP_H