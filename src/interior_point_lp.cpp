// filepath: /interior_point_lp/interior_point_lp/src/interior_point_lp.cpp
#include "interior_point_lp.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

InteriorPointLP::InteriorPointLP(const Eigen::MatrixXd& A_, const Eigen::VectorXd& b_, 
                                   const Eigen::VectorXd& c_, double eps_, int max_iter) 
    : A(A_), b(b_), c(c_), eps(eps_), max_iterations(max_iter) {
    gamma = 0.1;  // Default centering parameter
}

bool InteriorPointLP::initialize() {
    int n = c.size();
    int m = b.size();

    x = Eigen::VectorXd::Ones(n).array() + 1.0;  // Start with strictly positive x
    s = Eigen::VectorXd::Ones(n).array() + 1.0;  // Start with strictly positive s

    Eigen::MatrixXd AAT = A * A.transpose();
    AAT += Eigen::MatrixXd::Identity(m, m) * 1e-8;  // Regularization for numerical stability
    Eigen::LDLT<Eigen::MatrixXd> ldlt(AAT);

    if (ldlt.info() != Eigen::Success) {
        std::cerr << "Failed to factorize A*A^T during initialization" << std::endl;
        return false;
    }

    y = ldlt.solve(A * c);
    if (ldlt.info() != Eigen::Success) {
        std::cerr << "Failed to solve for initial y" << std::endl;
        return false;
    }

    mu = x.dot(s) / n;
    return true;
}

double InteriorPointLP::getPrimalInfeasibility() const {
    return (A * x - b).norm();
}

bool InteriorPointLP::solve() {
    if (!initialize()) {
        std::cerr << "Initialization failed" << std::endl;
        return false;
    }

    int n = x.size();
    for (int iter = 0; iter < max_iterations; iter++) {
        double primal_feasibility = (A * x - b).norm();
        double dual_feasibility = (c - A.transpose() * y - s).norm();
        double complementarity = x.dot(s) / n;

        std::cout << "Iteration " << iter
                  << ": Primal infeas = " << primal_feasibility
                  << ", Dual infeas = " << dual_feasibility
                  << ", Complementarity = " << complementarity << std::endl;

        if (primal_feasibility < eps && dual_feasibility < eps && complementarity < eps) {
            std::cout << "Converged in " << iter << " iterations!" << std::endl;
            return true;
        }

        mu = gamma * x.dot(s) / n;

        Eigen::VectorXd r_p = b - A * x;
        Eigen::VectorXd r_d = c - A.transpose() * y - s;
        Eigen::VectorXd r_c = -x.array() * s.array() + mu;

        Eigen::VectorXd dx, dy, ds;
        solveNewtonSystem(r_p, r_d, r_c, dx, dy, ds);

        double alpha_p = 1.0, alpha_d = 1.0;
        for (int i = 0; i < n; i++) {
            if (dx(i) < 0) alpha_p = std::min(alpha_p, -0.99 * x(i) / dx(i));
            if (ds(i) < 0) alpha_d = std::min(alpha_d, -0.99 * s(i) / ds(i));
        }

        x += alpha_p * dx;
        y += alpha_d * dy;
        s += alpha_d * ds;
    }

    std::cout << "Maximum iterations reached without convergence" << std::endl;
    return false;
}

void InteriorPointLP::solveNewtonSystem(const Eigen::VectorXd& r_p, const Eigen::VectorXd& r_d,
                                         const Eigen::VectorXd& r_c, Eigen::VectorXd& dx,
                                         Eigen::VectorXd& dy, Eigen::VectorXd& ds) {
    int n = x.size();

    Eigen::VectorXd d2 = x.array() / s.array();
    Eigen::MatrixXd AD = A * d2.asDiagonal();
    Eigen::MatrixXd schur_complement = AD * A.transpose();
    schur_complement += Eigen::MatrixXd::Identity(schur_complement.rows(), schur_complement.cols()) * 1e-8;  // Regularization

    Eigen::VectorXd rhs = r_p + AD * (r_d - (r_c.array() / s.array()).matrix());  // Fix type mismatch
    dy = schur_complement.ldlt().solve(rhs);

    ds = r_d - A.transpose() * dy;
    dx = (r_c.array() - x.array() * ds.array()) / s.array();
}

void InteriorPointLP::setTolerance(double tolerance) {
    eps = tolerance;
}

void InteriorPointLP::setMaxIterations(int max_iter) {
    max_iterations = max_iter;
}

Eigen::VectorXd InteriorPointLP::getPrimalSolution() const { return x; }
Eigen::VectorXd InteriorPointLP::getDualSolution() const { return y; }
Eigen::VectorXd InteriorPointLP::getDualSlacks() const { return s; }
double InteriorPointLP::getObjectiveValue() const { return c.dot(x); }