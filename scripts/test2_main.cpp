// VolEsti (volume computation and sampling library)

// Copyright (c) 2021 Vissarion Fisikopoulos
// Copyright (c) 2021 Apostolos Chalkis
// Copyright (c) 2021 Maios Papachristou

// Licensed under GNU LGPL.3, see LICENCE file


// Edited by HZ on 11.06.2020 - mute doctest.h
#include <fstream>
#include <iostream>

#include <boost/random.hpp>

#include "random_walks/random_walks.hpp"
#include "sampling/sampling.hpp"
#include "diagnostics/univariate_psrf.hpp"

// Add the spectrahedron headers
#include "convex_bodies/spectrahedra/spectrahedron.h"
#include "convex_bodies/spectrahedra/LMI.h"
#include "SDPAFormatManager.h"
// Add this include after the other spectrahedron headers
#include "convex_bodies/spectrahedra/const_spectrahedron_wrapper.h"

// Add necessary namespaces
template
<
    typename MT,
    typename WalkType,
    typename SpectrahedronType
>
MT get_samples_boundary_spectahedron(SpectrahedronType &S)
{
    typedef typename SpectrahedronType::PointType Point;
    typedef typename SpectrahedronType::NT NT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;

    unsigned int walkL = 10, numpoints = 10000, nburns = 0, d = S.dimension();
    RNGType rng(d);
    Point StartingPoint(d);
    std::list<Point> randPoints;

    std::cout << "Starting uniform sampling on the boundary of the spectrahedron." << std::endl;
    uniform_sampling_boundary<WalkType>(randPoints, S, rng, walkL, numpoints,
                                        StartingPoint, nburns);
    std::cout << "Finished sampling. Number of points: " << randPoints.size() << std::endl;

    MT samples(d, numpoints);
    unsigned int jj = 0;
    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit != randPoints.end(); rpit++, jj++)
    {
        samples.col(jj) = (*rpit).getCoefficients();
    }
    return samples;
}

// New helper function to write samples to an Excel-compatible CSV file
template<typename MT>
void write_samples_to_excel(const MT &samples, const std::string &filename) {
    std::ofstream out(filename);
    if(!out) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }
    for(int i = 0; i < samples.rows(); ++i) {
        for(int j = 0; j < samples.cols(); ++j) {
            out << samples(i, j) << (j < samples.cols()-1 ? "," : "");
        }
        out << "\n";
    }
    out.close();
    std::cout << "Written samples to " << filename << std::endl;
}

template <typename NT, typename WalkType = BRDHRWalk>
void sample_spectrahedron_boundary(){
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef ConstSpectrahedronWrapper<Point>  SpectrahedronType;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    
    std::string filepath = "/home/syauqirp/gsoc25/tp_library/volesti/test/SDP/sdp__2_8.txt";

    SpectrahedronType S;

    SdpaFormatManager<NT> sdpaManager;
    std::ifstream in(filepath);
    std::cout << "Loading spectrahedron from file: " << filepath << std::endl;
    if(!in){
        std::cerr << "Error: cannot open file " << filepath << std::endl;
        std::exit(1);
    }
    
    Point objFunction;
    sdpaManager.loadSDPAFormatFile(in, S, objFunction); // Load directly into the wrapper
    in.close();
    std::cout << "Spectrahedron loaded successfully." << std::endl;
    std::cout.flush();
    
    // Debug: check initialization of S
    int dim = S.dimension();
    std::cout << "DEBUG: Spectrahedron dimension: " << dim << std::endl;

    Point initialPoint(S.getLMI().dimension());
    S.set_interior_point(initialPoint);
    
    std::cout << "Starting boundary sampling." << std::endl;
    MT samples = get_samples_boundary_spectahedron<MT, WalkType, SpectrahedronType>(S);
    std::cout << "Boundary sampling completed." << std::endl;

    // Write the samples to an Excel-compatible CSV file
    write_samples_to_excel(samples, "/home/syauqirp/gsoc25/samples_output.csv");

    VT score = univariate_psrf<NT, VT>(samples);
    std::cout << "PSRF score: " << score.maxCoeff() << std::endl;

    if(score.maxCoeff() < 1.1) {
        std::cout << "PSRF test passed: score is below 1.1" << std::endl;
    } else {
        std::cout << "PSRF test failed: score is above 1.1" << std::endl;
    }
}

int main() {
    std::cout << "Starting spectrahedron boundary sampling program" << std::endl;
    try {
        sample_spectrahedron_boundary<double, BRDHRWalk>();
        std::cout << "Program completed successfully" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Unknown error occurred" << std::endl;
        return 1;
    }
    return 0;
}