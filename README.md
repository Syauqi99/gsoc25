# GeomScale in Google Summer of Code 2025

**Student Name:** Syauqi Rahmat Perdana

**Institution:** Politecnico di Milano 
 
**Project:** Randomized Solver for Semidefinite Programs

---

This is my answer to GSOC 25 GeomScale test for a **randomized solver for semidefinite programs**.

**Answer Results:**  
- **Test1:** Run using R scripts located in the `scripts` folder.  
- **Test2:** Built as an executable named `test2_main` in the `build` folder.  
- **Test3:** Built as an executable named `test3_main` in the `build` folder.  

**Dependencies:**  
- Eigen3  
- Volesti  

---

## Test Details

- **Easy (Test1):** Compile and run Volesti. Use the R extension to visualize sampling (with any sampling algorithm) in a polytope.  
  - Run the R script located in `scripts/test1.R`.  

- **Medium (Test2):** Extend the hit-and-run algorithm to sample from the boundary of the spectrahedron (feasible region of an SDP).  
  - Build and run the executable `test2_main` in the `build` folder.  

- **Hard (Test3):** Implement an interior point algorithm for linear programming.  
  - Build and run the executable `test3_main` in the `build` folder.  

---

## Build and Run Instructions

1. Create and navigate to the build folder:
   ```
   mkdir build && cd build
   ```
2. Run CMake and make:
   ```
   cmake ..
   make
   ```
3. Run the test executables:
   - For Boundary Spectrahedron Sampling (test2_main):  
     ```
     ./test2_main
     ```
   - For Interior Point Methods (test3_main):  
     - To use the default LP data file:
       ```
       ./test3_main
       ```
     - Or specify a custom LP data file path:
       ```
       ./test3_main /path/to/your/lpdata.txt
       ```

**Notes:**  
- All necessary data is stored in the `data` folders.  
- The script to generate visualizations and to run the feasible linear programming test is available in the `GSOC_25.ipynb` notebook.