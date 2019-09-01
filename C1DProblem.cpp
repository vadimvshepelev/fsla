#include "C1DProblem.h"

C1DProblem prNBtest = C1DProblem("NBtest", 1000., 0., 1.e9, 1000., 0., 1.e5, 0., 1., 0., .1, .3, 100, .9, "tt");
C1DProblem prToro1Idealtest = C1DProblem("Toro-1", 1.,           .75,      1.,    .125,       0.,      .1,  0., 1., 0.,   .2, .3, 100, .9, "tt");
C1DProblem prToro2Idealtest = C1DProblem("Toro-2", 1.,           -2.,      .4,      1.,       2.,      .4,  0., 1., 0.,  .15, .5, 100, .9, "tt");
C1DProblem prToro3Idealtest = C1DProblem("Toro-3", 1.,            0.,   1000.,      1.,       0.,     .01,  0., 1., 0., .012, .5, 100, .9, "tt");
C1DProblem prToro4Idealtest = C1DProblem("Toro-4", 5.99924,  19.5975, 460.894, 5.99242, -6.19633, 46.0950,  0., 1., 0., .035, .4, 100, .9, "tt");
C1DProblem prToro5Idealtest = C1DProblem("Toro-5", 1.,     -19.59745,   1000.,      1., -1959745,     .01,  0., 1., 0., .012, .8, 100, .9, "tt");
