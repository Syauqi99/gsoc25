#include <iostream>
#include <Eigen/Dense>
#include "../include/interior_point_lp.h"

void solveMeatloafProblem() {
    // Meatloaf problem in standard form:
    // minimize    80x + 60y
    // subject to   x + y - s1 = 1
    //              -0.05x + 0.07y + s2 = 0
    //              x, y, s1, s2 >= 0

    // Initialize problem data
    Eigen::MatrixXd A(2, 4);
    A << 1, 1, -1, 0,
         -0.05, 0.07, 0, 1;
    
    Eigen::VectorXd b(2);
    b << 1, 0;
    
    Eigen::VectorXd c(4);
    c << 80, 60, 0, 0;
    
    // Create solver instance
    InteriorPointLP solver(A, b, c);
    solver.setTolerance(1e-6);
    solver.setMaxIterations(100);
    
    // Solve the problem
    if (solver.solve()) {
        std::cout << "Optimal solution found!" << std::endl;
        std::cout << "x = [" << solver.getPrimalSolution().transpose() << "]" << std::endl;
        std::cout << "Optimal value = " << solver.getObjectiveValue() << std::endl;
    } else {
        std::cout << "Failed to find optimal solution." << std::endl;
        std::cout << "Final primal infeasibility: " << solver.getPrimalInfeasibility() << std::endl;
    }
}

int main() {
    std::cout << "Solving the meatloaf problem:" << std::endl;
    solveMeatloafProblem();
    
    return 0;
}
