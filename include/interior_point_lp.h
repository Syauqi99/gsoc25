// filepath: /interior_point_lp/interior_point_lp/src/interior_point_lp.h
#ifndef INTERIOR_POINT_LP_H
#define INTERIOR_POINT_LP_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

class InteriorPointLP {
private:
    Eigen::MatrixXd A;  // Constraint matrix
    Eigen::VectorXd b;  // Right-hand side
    Eigen::VectorXd c;  // Objective coefficients
    
    // Current solution
    Eigen::VectorXd x;  // Primal variables
    Eigen::VectorXd y;  // Dual variables (for equality constraints)
    Eigen::VectorXd s;  // Dual slack variables
    
    // Parameters
    double mu;          // Barrier parameter
    double gamma;       // Centering parameter
    double eps;         // Convergence tolerance
    int max_iterations; // Maximum number of iterations
    
public:
    InteriorPointLP(const Eigen::MatrixXd& A_, const Eigen::VectorXd& b_, 
                   const Eigen::VectorXd& c_,
                   double eps_ = 1e-8, int max_iter = 100);
    
    bool initialize();
    bool solve();
    void solveNewtonSystem(const Eigen::VectorXd& r_p, const Eigen::VectorXd& r_d,
                           const Eigen::VectorXd& r_c, Eigen::VectorXd& dx,
                           Eigen::VectorXd& dy, Eigen::VectorXd& ds);
    
    // Getters for solution and objective value
    Eigen::VectorXd getPrimalSolution() const;
    Eigen::VectorXd getDualSolution() const;
    Eigen::VectorXd getDualSlacks() const;
    double getObjectiveValue() const;

    void setTolerance(double tolerance);
    void setMaxIterations(int max_iter);

    double getPrimalInfeasibility() const;
};

#endif // INTERIOR_POINT_LP_H