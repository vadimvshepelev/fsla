#include "C1DProblem.h"

// Test for non-ideal (Mie-Gruneisen) EOS of Bolotova-Nigmatullin
C1DProblem prNBtest = C1DProblem("NBtestHLL", 100., 0., 1.e9, 1000., 0., 1.e5, 0., 1., 0., 100.e-6, .7, 1000, .9, "tt"); 
// 5 Toro tests
C1DProblem prToro1Idealtest = C1DProblem("Toro-1", 1.,           .75,      1.,    .125,       0.,      .1,  0., 1., 0.,   .2, .3, 100, .3, "tt");
C1DProblem prToro2Idealtest = C1DProblem("Toro-2", 1.,           -2.,      .4,      1.,       2.,      .4,  0., 1., 0.,  .15, .5, 100, .9, "tt");
C1DProblem prToro3Idealtest = C1DProblem("Toro-3", 1.,            0.,   1000.,      1.,       0.,     .01,  0., 1., 0., .012, .5, 100, .9, "tt");
C1DProblem prToro4Idealtest = C1DProblem("Toro-4", 5.99924,  19.5975, 460.894, 5.99242, -6.19633, 46.0950,  0., 1., 0., .035, .4, 100, .9, "tt");
C1DProblem prToro5Idealtest = C1DProblem("Toro-5-ex", 1.,     -19.59745,   1000.,      1., -19.59745,     .01,  0., 1., 0., .012, .8, 100, .3, "tt");

// 08.02.2020 First Riemann problem test for laser volume target problem for Al
C1DProblem prLaserVTAlIdealTest1 = C1DProblem("LaserVTAl1-0.05nm-cell", 2700., 0., 300.e9, 2., 0., 19.4872e9, -50.e-9, 50e-9, 0., 1.e-12, 0., 2000, .9, "tt");
// 08.02.2020 Second Riemann problem test for laser volume target problem for Al
C1DProblem prLaserVTAlIdealTest2 = C1DProblem("LaserVTAl2-0.05nm-cell", 2700., 0., 19.4872e9, 2700., 0., 300.e9, -100.e-9, 0., 0., 1.e-12, -50.e-9, 1000, .9, "tt");
// 3 tests for V. Denisenko article about Godunov method
C1DProblem prDenisenko1 = C1DProblem("Denisenko-1-8100", 2., 0., 2., 1., 0., 1., 0., 1., 0., .225, .5, 8100, .9, "tt");
C1DProblem prDenisenko2 = C1DProblem("Denisenko-2-8100", 1., -1., 1., 1., 1., 1., 0., 1., 0., .15, .5, 8100, .9, "tt");
C1DProblem prDenisenko3 = C1DProblem("Denisenko-3-8100", 3., 4., 2., 2., 2., 1., 0., 1., 0., .09, .5, 8100, .9, "tt");
