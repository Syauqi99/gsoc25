#include "interior_point_lp.h"
#include "lp_utils.h"
#include <cmath>

// Initialize static parameters
InteriorPointLP::Parameters InteriorPointLP::params;

InteriorPointLP::Result InteriorPointLP::solve(const Eigen::MatrixXd& A_orig, const Eigen::VectorXd& b_orig, const Eigen::VectorXd& c_orig) {
    Result result;
    result.success = false;

    if (A_orig.rows() != b_orig.size()) {
        throw std::invalid_argument("Matrix A rows must match vector b size");
    }
    if (A_orig.cols() != c_orig.size()) {
        throw std::invalid_argument("Matrix A columns must match vector c size");
    }

    // Create working copies of the inputs
    Eigen::MatrixXd A = A_orig;
    Eigen::VectorXd b = b_orig;
    Eigen::VectorXd c = c_orig;
    
    // Problem dimensions
    const int n = c.size();  // Number of variables
    const int m = b.size();  // Number of constraints
    
    // Simple status output
    std::cout << "Solving LP problem with " << n << " variables and " << m << " constraints" << std::endl;
    
    // Scale the problem if needed
    LPUtils::ScalingInfo scaling;
    if (params.use_scaling) {
        scaling = LPUtils::scaleLP(A, b, c);
    }
    
    // Initialize variables
    Eigen::VectorXd x(n), lambda(m), s(n);
    
    // Compute initial point
    computeInitialPoint(A, b, c, x, lambda, s);
    
    // Main iteration loop
    int iter = 0;
    while (iter < params.max_iter) {
        // Check for NaN values
        if (LPUtils::containsNanOrInf(x) || LPUtils::containsNanOrInf(lambda) || LPUtils::containsNanOrInf(s)) {
            break;
        }
        
        // Current duality measure
        double mu = x.dot(s) / n;
        
        // Compute residuals
        Eigen::VectorXd rc = A.transpose() * lambda + s - c; // Dual residual
        Eigen::VectorXd rb = A * x - b;                     // Primal residual
        
        // Check convergence
        if (checkConvergence(x, lambda, s, rc, rb, mu, b_orig, c_orig)) {
            result.success = true;
            result.x = x;
            result.optimal_value = c_orig.dot(x);
            result.primal_infeas = rb.norm() / (1.0 + b_orig.norm());
            result.dual_infeas = rc.norm() / (1.0 + c_orig.norm());
            result.gap = mu;

            std::cout << "Converged after " << iter << " iterations." << std::endl;
            
            // Rescale solution if needed
            if (params.use_scaling) {
                LPUtils::rescaleSolution(result.x, lambda, s, scaling);
            }
            
            std::cout << "Optimal value: " << result.optimal_value << std::endl;
            std::cout << "Optimal solution (x): " << result.x.transpose() << std::endl;
            
            return result;
        }
        
        // Step 1: Compute affine scaling direction (predictor)
        Eigen::VectorXd dx_aff(n), dlambda_aff(m), ds_aff(n);
        computeAffineDirection(A, x, lambda, s, dx_aff, dlambda_aff, ds_aff, rc, rb);
        
        // Step 2: Compute step lengths for affine direction
        double alpha_pri_aff, alpha_dual_aff;
        computeStepLengths(x, s, dx_aff, ds_aff, alpha_pri_aff, alpha_dual_aff);
        
        // Step 3: Compute centering parameter
        double sigma = computeCenteringParameter(x, s, dx_aff, ds_aff, alpha_pri_aff, alpha_dual_aff, mu);
        
        // Step 4: Compute combined direction (corrector)
        Eigen::VectorXd dx(n), dlambda(m), ds(n);
        computeCombinedDirection(A, x, lambda, s, dx_aff, ds_aff, dx, dlambda, ds, rc, rb, sigma, mu);
        
        // Step 5: Compute step lengths for combined direction
        double alpha_pri, alpha_dual;
        computeStepLengths(x, s, dx, ds, alpha_pri, alpha_dual);
        
        // Apply step length scaling factor
        double eta_factor = params.eta;
        if (n > 1000) {
            eta_factor = std::min(0.7, eta_factor);
        }
        alpha_pri = std::min(1.0, eta_factor * alpha_pri);
        alpha_dual = std::min(1.0, eta_factor * alpha_dual);
        
        // Step 6: Update variables
        x = x + alpha_pri * dx;
        lambda = lambda + alpha_dual * dlambda;
        s = s + alpha_dual * ds;
        
        iter++;
    }
    
    std::cout << "Maximum iterations reached. Solution may not be optimal." << std::endl;
    
    // Rescale solution if needed
    if (params.use_scaling) {
        LPUtils::rescaleSolution(x, lambda, s, scaling);
    }
    
    // Recompute final residuals and duality measure
    Eigen::VectorXd rb = A * x - b;
    Eigen::VectorXd rc = A.transpose() * lambda + s - c;
    double mu = x.dot(s) / n;
    
    result.x = x;
    result.optimal_value = c_orig.dot(x);
    result.primal_infeas = rb.norm() / (1.0 + b_orig.norm());
    result.dual_infeas = rc.norm() / (1.0 + c_orig.norm());
    result.gap = mu;

    std::cout << "Current value: " << result.optimal_value << std::endl;
    std::cout << "Current solution (x): " << result.x.transpose() << std::endl;

    return result;
}

void InteriorPointLP::computeInitialPoint(
    const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& c,
    Eigen::VectorXd& x, Eigen::VectorXd& lambda, Eigen::VectorXd& s) 
{
    const int n = c.size();
    const int m = b.size();
    
    x = Eigen::VectorXd::Ones(n);
    lambda = Eigen::VectorXd::Zero(m);
    s = c - A.transpose() * lambda;
    
    double min_s = s.minCoeff();
    if (min_s <= 1e-2) {
        s = s - (min_s - 1.0) * Eigen::VectorXd::Ones(n);
    }
    
    double scaling_target = 1.0;
    if (n > 1000) {
        scaling_target = 10.0;
    }
    
    Eigen::VectorXd xs = x.cwiseProduct(s);
    double geo_mean = std::pow(xs.prod(), 1.0/n);
    
    if (!std::isfinite(geo_mean) || geo_mean < 1e-10) {
        geo_mean = 1.0;
    }
    
    double scale_factor = std::sqrt(scaling_target/geo_mean);
    x = x * scale_factor;
    s = s / scale_factor;
}

void InteriorPointLP::computeAffineDirection(
    const Eigen::MatrixXd& A, 
    const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, const Eigen::VectorXd& s,
    Eigen::VectorXd& dx_aff, Eigen::VectorXd& dlambda_aff, Eigen::VectorXd& ds_aff,
    const Eigen::VectorXd& rc, const Eigen::VectorXd& rb) 
{
    Eigen::VectorXd rhs1 = -rc;
    Eigen::VectorXd rhs2 = -rb;
    Eigen::VectorXd rhs3 = -(x.array() * s.array()).matrix();
    
    solveLinearSystem(A, x, s, rhs1, rhs2, rhs3, dx_aff, dlambda_aff, ds_aff);
}

void InteriorPointLP::computeStepLengths(
    const Eigen::VectorXd& x, const Eigen::VectorXd& s,
    const Eigen::VectorXd& dx, const Eigen::VectorXd& ds,
    double& alpha_pri, double& alpha_dual) 
{
    alpha_pri = 1.0;
    alpha_dual = 1.0;
    
    const double STEP_THRESHOLD = -1e-12;
    
    for (int i = 0; i < x.size(); i++) {
        if (dx(i) < STEP_THRESHOLD) {
            double ratio = -x(i) / dx(i);
            if (std::isfinite(ratio) && ratio < alpha_pri) {
                alpha_pri = ratio;
            }
        }
    }
    
    for (int i = 0; i < s.size(); i++) {
        if (ds(i) < STEP_THRESHOLD) {
            double ratio = -s(i) / ds(i);
            if (std::isfinite(ratio) && ratio < alpha_dual) {
                alpha_dual = ratio;
            }
        }
    }
    
    if (!std::isfinite(alpha_pri) || alpha_pri < 0) alpha_pri = 0;
    if (!std::isfinite(alpha_dual) || alpha_dual < 0) alpha_dual = 0;
}

double InteriorPointLP::computeCenteringParameter(
    const Eigen::VectorXd& x, const Eigen::VectorXd& s,
    const Eigen::VectorXd& dx_aff, const Eigen::VectorXd& ds_aff,
    double alpha_pri_aff, double alpha_dual_aff, double mu) 
{
    const int n = x.size();
    
    double mu_aff = 0.0;
    int count = 0;
    
    for (int i = 0; i < n; i++) {
        double x_new = x(i) + alpha_pri_aff * dx_aff(i);
        double s_new = s(i) + alpha_dual_aff * ds_aff(i);
        
        if (x_new > 0 && s_new > 0) {
            mu_aff += x_new * s_new;
            count++;
        }
    }
    
    if (count == 0) {
        mu_aff = mu;
    } else {
        mu_aff /= count;
    }
    
    if (mu_aff < 1e-14) mu_aff = 1e-14;
    if (mu < 1e-14) mu = 1e-14;
    
    double sigma = std::pow(mu_aff / mu, 3);
    
    if (sigma < 0.01) sigma = 0.01;
    if (sigma > 0.5) sigma = 0.5;
    
    return sigma;
}

void InteriorPointLP::computeCombinedDirection(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, const Eigen::VectorXd& s,
    const Eigen::VectorXd& dx_aff, const Eigen::VectorXd& ds_aff,
    Eigen::VectorXd& dx, Eigen::VectorXd& dlambda, Eigen::VectorXd& ds,
    const Eigen::VectorXd& rc, const Eigen::VectorXd& rb,
    double sigma, double mu) 
{
    const int n = x.size();
    
    Eigen::VectorXd corrector(n);
    
    for (int i = 0; i < n; i++) {
        if (std::fabs(dx_aff(i)) > 1e6 || std::fabs(ds_aff(i)) > 1e6) {
            corrector(i) = 0;
        } else {
            corrector(i) = dx_aff(i) * ds_aff(i);
        }
    }
    
    Eigen::VectorXd rhs1 = -rc;
    Eigen::VectorXd rhs2 = -rb;
    Eigen::VectorXd rhs3 = -(x.array() * s.array()).matrix() - corrector + sigma * mu * Eigen::VectorXd::Ones(n);
    
    solveLinearSystem(A, x, s, rhs1, rhs2, rhs3, dx, dlambda, ds);
}

bool InteriorPointLP::solveLinearSystem(
    const Eigen::MatrixXd& A,
    const Eigen::VectorXd& x, const Eigen::VectorXd& s,
    const Eigen::VectorXd& rhs1, const Eigen::VectorXd& rhs2, const Eigen::VectorXd& rhs3,
    Eigen::VectorXd& dx, Eigen::VectorXd& dlambda, Eigen::VectorXd& ds) 
{
    const int n = x.size();
    const int m = A.rows();
    
    Eigen::VectorXd d(n);
    for (int i = 0; i < n; i++) {
        if (s(i) < 1e-14) {
            d(i) = x(i) / 1e-14;
        } else {
            d(i) = x(i) / s(i);
        }
        
        if (d(i) < 1e-12) d(i) = 1e-12;
        if (d(i) > 1e12) d(i) = 1e12;
    }
    
    Eigen::MatrixXd D = d.asDiagonal();
    
    Eigen::MatrixXd AD = A * D;
    Eigen::MatrixXd M = AD * A.transpose();
    
    for (int i = 0; i < M.rows(); i++) {
        M(i, i) += params.regularization * (1.0 + M(i, i));
    }
    
    Eigen::VectorXd rhs_temp(n);
    for (int i = 0; i < n; i++) {
        if (std::fabs(s(i)) < 1e-14) {
            rhs_temp(i) = rhs3(i) / 1e-14;
        } else {
            rhs_temp(i) = rhs3(i) / s(i);
        }
        
        if (rhs_temp(i) < -1e12) rhs_temp(i) = -1e12;
        if (rhs_temp(i) > 1e12) rhs_temp(i) = 1e12;
    }
    
    Eigen::VectorXd rhs_lambda = rhs2 - A * (D * rhs1 + rhs_temp);
    
    if (LPUtils::containsNanOrInf(rhs_lambda)) {
        return false;
    }
    
    try {
        Eigen::LDLT<Eigen::MatrixXd> ldlt(M);
        if (ldlt.info() != Eigen::Success) {
            return false;
        }
        
        dlambda = ldlt.solve(rhs_lambda);
        
        if (LPUtils::containsNanOrInf(dlambda)) {
            return false;
        }
    }
    catch (const std::exception&) {
        return false;
    }
    
    dx = D * (A.transpose() * dlambda - rhs1) + rhs_temp;
    
    ds.resize(n);
    for (int i = 0; i < n; i++) {
        if (x(i) < 1e-12) {
            ds(i) = (rhs3(i) - s(i) * dx(i)) / 1e-12;
        } else {
            ds(i) = (rhs3(i) - s(i) * dx(i)) / x(i);
        }
        
        if (ds(i) < -1e12) ds(i) = -1e12;
        if (ds(i) > 1e12) ds(i) = 1e12;
    }
    
    return true;
}

bool InteriorPointLP::checkConvergence(
    const Eigen::VectorXd& x, const Eigen::VectorXd& lambda, const Eigen::VectorXd& s,
    const Eigen::VectorXd& rc, const Eigen::VectorXd& rb, double mu,
    const Eigen::VectorXd& b_orig, const Eigen::VectorXd& c_orig) 
{
    double primal_infeas = rb.norm() / (1.0 + b_orig.norm());
    double dual_infeas = rc.norm() / (1.0 + c_orig.norm());
    double gap = mu;
    
    return (primal_infeas < params.tol && dual_infeas < params.tol && gap < params.tol);
}

void InteriorPointLP::setParameters(const Parameters& p) {
    params = p;
}
