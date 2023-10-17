#include "F1DReconstruction.h"

#include <algorithm>
#include <array>
// #include <execution>
// #include <functional>
#include <numeric>
#include <ranges>
// #include <span>

// #include "Eigen/Eigen/Dense"

#include "_matrix4.h"
#include "_vector4.h"

// ENO2 stencils
const double _rm1r[] = {-1./2., 3./2.,  0.   };
const double _r0r[]  = {    0., 1./2.,  1./2.};
const double _rm1l[] = { 1./2., 1./2.,  0.   };
const double _r0l[]  = {    0., 3./2., -1./2.};


// Possible ENO-3 stencils
//std::array<std::array<const double, 6>, 3> _wm_nminus {{
//	{{0., 0., 0., 11./6, -7./6, 2./6}},
//	{{0., 0., 2./6, 5./6, -1./6, 0.}},
//	{{0., -1./6, 5./6, 2./6, 0., 0.}}
//}};

//std::array<std::array<const double, 6>, 3> _wm_nplus {{
//	{{2./6, -7./6, 11./6, 0., 0., 0.}},
//	{{0., -1./6, 5./6, 2./6, 0., 0.}},
//	{{0., 0., 2./6, 5./6, -1./6, 0.}}
//}};


F1DReconstruction::F1DReconstruction(C1DField& _fld) {
	// std::size_t i=0;
	double _v[] = {0., 0., 0., 0.};
	// vector<double> tempVect = vector<double>(_v, _v+3);
	for (std::size_t i = 0; i < _fld.U.size(); ++ i) {
		// ULx.push_back(tempVect);  // vector<vector<double>>
		ULx.push_back(Vector4::ZERO);  // vector<Vector4>
		// URx.push_back(tempVect);
		URx.push_back(Vector4::ZERO);
	}
}


F1DENO2Reconstruction::F1DENO2Reconstruction(C1DField& _fld)
	: F1DReconstruction(_fld),
	  rm1r(vector<double>(_rm1r, _rm1r+3)),
	  r0r(vector<double>(_r0r, _r0r+3)),
	  rm1l(vector<double>(_rm1l, _rm1l+3)),
	  r0l(vector<double>(_r0l, _r0l+3)) {
/*	rm1r = vector<double>(_rm1r, _rm1r+3);		// Для типа vector<T> при инициализация через два параметра первый это указатель на первый элемент,
	r0r = vector<double>(_r0r, _r0r+3);			// а второй -- указатель на последний элемент. Таким образом очень удобно инициализировать вектор
	rm1l = vector<double>(_rm1l, _rm1l+3);		// через статический массив, который, в свою очередь можно инициализировать списком значений {..., ..., ... etc.}
	r0l = vector<double>(_r0l, _r0l+3);			// Пока не дошли до C11, едем на этом. :) */
}


F1DENO3Reconstruction::F1DENO3Reconstruction(C1DField& _fld)
	: F1DENO2Reconstruction(_fld) {}



F1DWENO5Reconstruction::F1DWENO5Reconstruction(C1DField& _fld)
	: F1DENO3Reconstruction(_fld) {
	DISCRETE_LAMBDA5 = prediscretizeWENO5LambdaMapping(
				1000000/*00*/, 1./3.);
}


F1DCharWiseWENO5Reconstruction::F1DCharWiseWENO5Reconstruction(
		C1DField& _fld, FEOS& _eos)
	: F1DWENO5Reconstruction(_fld), eos(_eos) {}


void F1DENO2Reconstruction::calc(C1DField& fld) /*override*/ {
	int i = 0, n = 0;
	const int nComp = 4;
	int imin = fld.imin, imax = fld.imax;
	auto&& U = fld.U;
	Vector4 diffPlus = Vector4::ZERO, diffMinus=Vector4::ZERO;
	int stencilIndex[] = {0, 0, 0, 0};
	for (i = imin - 1; i < imax + 1; ++ i) {
		diffPlus = Vector4(U[i+1][0], U[i+1][1], U[i+1][2], 0.)
				- Vector4(U[i][0], U[i][1], U[i][2], 0.);

		diffMinus = Vector4(U[i][0], U[i][1], U[i][2], 0.)
				- Vector4(U[i-1][0], U[i-1][1], U[i-1][2], 0.);

		for (n = 0; n < nComp - 1; ++ n) {
		    // Choose stencil
			if (fabs(diffMinus[n]) <= fabs(diffPlus[n]))
				stencilIndex[n] = i - 1; else stencilIndex[n] = i;

			// Calculate ENO-2 approximations themselves
			if (stencilIndex[n] == i - 1) {
				URx[i][n] = rm1r[0]*U[i-1][n]
						+ rm1r[1]*U[i][n]
						+ rm1r[2]*U[i+1][n];

				ULx[i][n] = rm1l[0]*U[i-1][n]
						+ rm1l[1]*U[i][n]
						+ rm1l[2]*U[i+1][n];
			} else {
				URx[i][n] = r0r[0]*U[i-1][n]
						+ r0r[1]*U[i][n]
						+ r0r[2]*U[i+1][n];

				ULx[i][n] = r0l[0]*U[i-1][n]
						+ r0l[1]*U[i][n]
						+ r0l[2]*U[i+1][n];
			}
		}
	}
}


short chooseENO3Stencil(const std::ranges::sized_range auto&& u_stencil) {
	/* Choose the best ENO3 stencil for the reconstruction (the smoothest
	 * result with the smallest variation (finite difference).
	 *
	 * Return a number from {0, 1, 2} that represent respectively
	 * 	0 <-> [j-2, j-1, j+0          ],
	 * 	1 <-> [     j-1, j+0, j+1     ],
	 * 	2 <-> [          j+0, j+1, j+2].
	 */

	if (std::fabs(u_stencil[2] - u_stencil[1])
			<= std::fabs(u_stencil[3] - u_stencil[2])) {
		if (std::fabs(u_stencil[2] - 2. * u_stencil[1] + u_stencil[0])
				<= std::fabs(u_stencil[3] - 2. * u_stencil[2] + u_stencil[1])
				) {
			return 0;
		} else {
			return 1;
		}
	}

	if (std::fabs(u_stencil[3] - 2. * u_stencil[2] + u_stencil[1])
			<= std::fabs(u_stencil[4] - 2. * u_stencil[3] + u_stencil[2]))
		return 1;

	return 2;
}


double computeENO3ReconstructionKernel(
		const std::ranges::sized_range auto&& u_stencil,
		short which_stencil) {
	/* 3rd order ENO reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 * Choose a stencil with `which_stencil`:
	 * 	0: [j-2, j-1, j+0          ],
	 * 	1: [     j-1, j+0, j+1     ],
	 * 	2: [          j+0, j+1, j+2].
	 */
	//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}}

	switch (which_stencil) {
		case 0: return (2. * u_stencil[0]
						- 7. * u_stencil[1]
						+ 11. * u_stencil[2]) / 6.; break;
		case 1: return (-1. * u_stencil[1]
						+ 5. * u_stencil[2]
						+ 2. * u_stencil[3]) / 6.; break;
		case 2: return (2. * u_stencil[2]
						+ 5. * u_stencil[3]
						- 1. * u_stencil[4]) / 6.; break;
		default: return 0.;
	}
}


double computeENO3ReconstructionKernelRev(
		const std::ranges::sized_range auto&& u_stencil,
		short which_stencil) {
	/* 3rd order ENO reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_minus or reversed f_plus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '-'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '+').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 * Choose a stencil with `which_stencil`:
	 * 	2: [j+3, j+2, j+1          ],
	 * 	1: [     j+2, j+1, j+0     ],
	 * 	0: [          j+1, j+0, j-1].
	 */
	//		{{0., 0., 0., 11./6, -7./6, 2./6}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}}

	switch (which_stencil) {
		case 0: return (-1. * u_stencil[0]
						+ 5. * u_stencil[1]
						+ 2. * u_stencil[2]) / 6.; break;
		case 1: return (2. * u_stencil[1]
						+ 5. * u_stencil[2]
						- 1. * u_stencil[3]) / 6.; break;
		case 2: return (11. * u_stencil[2]
						- 7. * u_stencil[3]
						+ 2. * u_stencil[4]) / 6.; break;
		default: return 0.;
	}
}


void F1DENO3Reconstruction::calcComponent_(
		const std::ranges::common_range auto&& u,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells /* = 3*/) {
	/* Simple finite-volume essentially non-oscillatory 3-rd order
	 * (FV ENO-3) reconstruction for a 1-D field `u`.
	 *
	 * Fill `u_plus_rec` and `u_minus_rec` with the reconstructed
	 * values.
	 *
	 * `u` is assumed to contain precisely `n_ghost_cells` ghost nodes
	 * on the left and on the right. It should be at least 3 !!
	 */

	const std::size_t stencil_size = 5;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = _actual_stencil_size / 2;  // 3

	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{
		mini - 1,
		maxi + 1
	};

	double uhatplus = 0.;
	double uhatminus = 0.;

	auto&& j_it_p = std::ranges::begin(u);  // f_plus
	// auto j_it_m = std::ranges::begin(u);  // f_minus

	// std::vector<T> u_plus(_actual_stencil_size);
	// std::vector<T> u_minus(_actual_stencil_size);

	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size - 1);
	auto&& u_plus = std::ranges::views::counted(j_it_p, 6);
	auto&& u_minus = u_plus | std::ranges::views::reverse;
	short which_stencil = chooseENO3Stencil(
				std::ranges::views::counted(std::ranges::begin(u_plus), 5)
				);

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(u);  // u_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size - 1);
		u_plus = std::ranges::views::counted(j_it_p, 6);
		u_minus = u_plus | std::ranges::views::reverse;

		uhatplus = computeENO3ReconstructionKernel(
					std::ranges::views::counted(
						std::ranges::begin(u_plus), 5),
					which_stencil
					);

		which_stencil = chooseENO3Stencil(
					std::ranges::views::counted(
						std::ranges::begin(u_plus) + 1, 5)
					);

//		uhatminus = computeENO3ReconstructionKernelRev(
//					std::ranges::views::counted(
//						std::ranges::begin(u_plus)+1, 5),
//					which_stencil
//				);
		uhatminus = computeENO3ReconstructionKernel(
					std::ranges::subrange(
						std::begin(u_minus),
						std::end(u_minus) - 1
					),
					2 - which_stencil
					);

		u_plus_rec[j] = uhatplus;
		u_minus_rec[j] = uhatminus;

//		std::cout << (u_minus_rec[j]-u_plus_rec[j])*(u[j]-u[j+1])
//					<< "\n";
	}
}


void F1DENO3Reconstruction::calc(C1DField& fld) /*override*/ {
	/* Perform finite-volume component-wise essentially non-oscillatory
	 * 3-rd order (FV c'mp't-wise ENO-3) reconstruction for a 1-D Euler
	 * equation field, i. e. for any field 3-component field `C1DField`.
	 */

	auto&& U = std::views::all(fld.U);
	auto&& u_plus = std::views::all(URx);
	auto&& u_minus = std::views::all(ULx) | std::ranges::views::drop(1);
	auto&& components = {
		&Vector4::x
		, &Vector4::y
		, &Vector4::z
		// , &Vector4::w
	};
	std::for_each(
//			std::execution::par_unseq,
			std::ranges::begin(components),
			std::ranges::end(components),
			[&](auto&& kth_vector_component) {
		F1DENO3Reconstruction::calcComponent_(
			U       | std::ranges::views::transform(kth_vector_component),
			u_plus  | std::ranges::views::transform(kth_vector_component),
			u_minus | std::ranges::views::transform(kth_vector_component),
			3
			// n_size = fld.imax - fld.imin + 1
			);
	});
}


template <typename T1, typename T2>
std::valarray<double> vecMatDot(const T1& vec, const T2& mat) {
	/* Multiply a vector by a matrix (a sum product over
	 * the last axis of mat and vec).
	 */

	std::valarray<double> result(std::ranges::size(mat));

	for (std::size_t k = 0; k < std::ranges::size(mat)
			&& k < std::ranges::size(vec); ++ k)
		result[k] = std::inner_product(std::ranges::begin(mat[k]),
									   std::ranges::end(mat[k]),
									   std::ranges::begin(vec), 0.);

	return result;
}


//std::array<double, 3> smoothness_indicators(const T1& f_stencil) {
std::valarray<double> betaSmoothnessIndicators(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO5 smoothness indicators of Jiang and Shu (1996)
	 * for each of the 3 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 3 local interpolating
	 * polynomials in the sub-stencils of 5-node `f_stencil`.
	 *
	 * This allows for (2*3-1)=5th order accuracy from the 3rd
	 * order Eno schemes.
	 */

	// std::array<double, 3> res;
	double beta1;
	double beta2;
	double beta3;

	beta1 = ((13./12.) * std::pow(
					f_stencil[0] - 2.*f_stencil[1] + f_stencil[2], 2)
				+ (1./4.) * std::pow(
					f_stencil[0] - 4.*f_stencil[1] + 3.*f_stencil[2], 2
					));

	beta2 = ((13./12.) * std::pow(
					f_stencil[1] - 2.*f_stencil[2] + f_stencil[3], 2)
				+ (1./4.) * std::pow((f_stencil[1] - f_stencil[3]), 2));

	beta3 = ((13./12.) * std::pow(
					f_stencil[2] - 2.*f_stencil[3] + f_stencil[4], 2)
				+ (1./4.) * std::pow(
					3.*f_stencil[2] - 4.*f_stencil[3] + f_stencil[4], 2
					));

	return {beta1, beta2, beta3};
}


std::valarray<double> betaSmoothnessIndicatorsZQ(
		const std::ranges::sized_range auto& f_stencil) {
	/* Return the WENO5-ZQ smoothness indicators of  and Shu (1996)
	 * for each of the 3 substencils.
	 * That is the sum of the normalized squares of the scaled
	 * L2-norms of all the derivatives of 3 local interpolating
	 * polynomials in the sub-stencils of 5-node `f_stencil`.
	 *
	 * This allows for (2*3-1)=5th order accuracy.
	 */

	// std::array<double, 3> res;
	double beta1;
	double beta2;
	double beta3;

	beta1 = std::pow(-82. * f_stencil[1] + 11. * f_stencil[0]
			 + 82. * f_stencil[3] - 11. * f_stencil[4]
			 + 2. * f_stencil[1] - f_stencil[0]
			 - 2. * f_stencil[3] + f_stencil[4], 2) * (1./14400.)
			+ std::pow((40. * f_stencil[1] - 3. * f_stencil[0]
			  - 74. * f_stencil[2]
			  + 40. * f_stencil[3] - 3. * f_stencil[4], 2) * (1./56.)
			 + (123./455.) * (-4. * f_stencil[1] + f_stencil[2]
			  + 6. * f_stencil[2] - f_stencil[3] + f_stencil[4]) * (1./24.),
			2) * (3./13.)
			+ std::pow((2. * f_stencil[1] - f_stencil[0]
			 - 2. * f_stencil[1] + f_stencil[0]) * (1./12.), 2) * (781./20.)
			+ std::pow((-4. * f_stencil[1] + f_stencil[0]
			 + 6. * f_stencil[2]
			 - 4. * f_stencil[3] + f_stencil[4]) * (1./24.),
			2) * (1421461./2275.);

	beta2 = std::pow(f_stencil[1] - f_stencil[2], 2);

	beta3 = std::pow(f_stencil[2] - f_stencil[3], 2);

	return {beta1, beta2, beta3};
}


//template <ArithmeticWith<numeric_val> T>
//std::valarray<double> betaSmoothnessIndicatorsMat(
//		const std::ranges::common_range auto& f_stencil,
//		std::array<std::array<const T, 6>, 6> _coefs[] = plus_coefs) {
//	/* Return the smoothness indicators beta_k, k=0,1,2
//	 * for each of the 3 substencils of `f_stencil`.
//	 */

//	std::valarray<double> res(3);

//	for (std::size_t k = 0; k < 3; ++ k) {
//		// for (std::size_t k = 0; k < half_size + 1; ++ k)
//		res[k] = std::inner_product(
//			std::ranges::begin(f_stencil),
//			std::ranges::end(f_stencil),
//			std::ranges::begin(vecMatDot(f_stencil, _coefs[k])),
//			0.
//		);
//	}


//	return res;
//}


std::valarray<double> f3OrdReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 3rd order reconstructions of f(j) from all the 3 3-element
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 */
	//		{{2./6, -7./6, 11./6, 0., 0., 0.}},
	//		{{0., -1./6, 5./6, 2./6, 0., 0.}},
	//		{{0., 0., 2./6, 5./6, -1./6, 0.}}

	std::valarray<double> q_res(3);

	q_res[0] = (2. * f_stencil[0]
			  - 7. * f_stencil[1]
			 + 11. * f_stencil[2]) / 6.;

	q_res[1] = (-1. * f_stencil[1]
			   + 5. * f_stencil[2]
			   + 2. * f_stencil[3]) / 6.;

	q_res[2] = (2. * f_stencil[2]
			  + 5. * f_stencil[3]
			  - 1. * f_stencil[4]) / 6.;

	return q_res;
}


std::valarray<double> f4OrdZQReconstructionFromStencil(
		const std::ranges::sized_range auto& f_stencil) {
	/* 4th order reconstruction of f(j) from
	 * substencils of `f_stencil` (f_plus or reversed f_minus:
	 * receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *               (or [j+3, j+2, j+1, j+0, j-1, ...] for '-').
	 *                     ^    ^    ^    ^    ^    ^
	 *                     0    1    2    3    4    |
	 */

	std::valarray<double> q_res(3);

	q_res[0] = f_stencil[2] + (
				-82. * f_stencil[1] + 11. * f_stencil[0]
				+ 82. * f_stencil[3] - 11. * f_stencil[4]
			) * (1./120.) * 0.5
			+ (40. * f_stencil[1] - 3. * f_stencil[0]
				- 74. * f_stencil[2]
				+ 40. * f_stencil[3] - 3. * f_stencil[4]
			) * (1./56.) * (0.25 - (1./12.))
			+ (2. * f_stencil[1] - f_stencil[0]
				- 2. * f_stencil[3] + f_stencil[4]
			) * (1./12.) * (0.125 - (3./20.) * 0.5)
			+ (-4. * f_stencil[1] + f_stencil[0]
				+ 6. * f_stencil[2]
				-4. * f_stencil[3] + f_stencil[4]
			) * (1./24.) * (0.0625 - (3./14.) * 0.25 + (3./560.));

	q_res[1] = f_stencil[2] + (f_stencil[2] - f_stencil[1]) * 0.5;

	q_res[2] = f_stencil[2] + (f_stencil[3] - f_stencil[2]) * 0.5;

	return q_res;
}


double henrickGMappingForLambda(double lambda_weno_weight,
						   double lambda_ideal = 1./3.) {
	/* The mapping function g by Henrick modified for symmetric
	 * lambda-weights by Zheng Hong, Zhengyin Ye and Kun Ye.
	 */

	double square_ideal = lambda_ideal * lambda_ideal;

	return lambda_weno_weight
			* (lambda_ideal
			   + square_ideal
			   - 3. * lambda_ideal * lambda_weno_weight
			   + lambda_weno_weight * lambda_weno_weight)
			/ (square_ideal
			   + lambda_weno_weight
			   * (1. - 2. * lambda_ideal));
}


//double henrickGMapping(double omega_weno_weight,
//				  double d_ideal = 1./3.) {
//	/* The original mapping function g by Henrick et al. */

//	double d_square = d_ideal * d_ideal;
//	return omega_weno_weight
//			* (d_ideal
//			   + d_square
//			   - 3. * d_ideal * omega_weno_weight
//			   + omega_weno_weight * omega_weno_weight)
//			/ (d_square
//			   + omega_weno_weight
//			   * (1. - 2.*d_ideal));
//}


// const std::ranges::common_range auto&&
double alphaWENO5FMWeight(
		double beta_IS_coefficient,
		double epsilon = 1e-40,
		double p = 2.) {
	/* Compute appropriate alpha(α)-weights for the WENO5-FM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights.
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. / std::pow(epsilon + beta_IS_coefficient, p);
}


double alphaWENO5ZMWeight(
		double beta_IS_coefficient,
		double tau_5,
		double epsilon = 1e-40,
		double p = 2.) {
	/* Compute appropriate alpha(α)-weights for the WENO5-ZM scheme,
	 * by which the inverses of smoothness indicators are meant,
	 * so inverse beta(β) with the caveat of aritificially finite
	 * answers using the added epsilon-parameter to the beta-weights.
	 *
	 * `p` controls (increases) the amount of numerical dissipation
	 * (it's recommended to take it = r-1 for 2r-1 order schemes,
	 * so 2 for WENO5).
	 */

	return 1. + std::pow(tau_5 / (beta_IS_coefficient + epsilon), p);
}


void lambdaWENO5FMWeights(
		const std::ranges::common_range auto&& alpha_weights,
		std::ranges::common_range auto&& res) {
	/* FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
	 * for WENO5-FM or WENO5-ZM
	 * due to Zheng Hong, Zhengyin Ye and Kun Ye:
	 * lambda_weights = alpha_weights / alpha_weights.sum();
	 */

	// return alpha_weights / alpha_weights.sum();
	double sum = std::reduce(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				0.);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(alpha_weights),
				std::ranges::end(alpha_weights),
				std::ranges::begin(res),
				[sum](const auto alpha) {
		return alpha / sum;
	});
}


std::vector<
	double> F1DWENO5Reconstruction::prediscretizeWENO5LambdaMapping(
		std::size_t N,
		double lambda_ideal = 1./3.) {
	/* Pre-discrete mapping method for WENO5-FM
	 * (idea due to Hong et al.) to increase performace.
	 * Construct the lambda mapping `valarray` of size `N + 1`
	 * for lambdas in the range [0..1].
	 */

	std::vector<double> res_lookup_table(N + 1);

//	auto ns = std::ranges::common_view(
//			std::ranges::views::iota(std::size_t(0))
//				| std::views::take(N + 1)
//	);

//	std::transform(
//				std::execution::par_unseq,
//				std::ranges::begin(ns),
//				std::ranges::end(ns),
//				std::ranges::begin(res_lookup_table),
//				[N](std::size_t n) {
//		return henrickGMappingForLambda<T>(
//					static_cast<T>(n) / static_cast<T>(N));
//	});

	for (std::size_t n = 0; n <= N; ++ n)
		res_lookup_table[n] = henrickGMappingForLambda(
					static_cast<double>(n) / static_cast<double>(N),
					lambda_ideal);

	return res_lookup_table;
}


//std::ranges::common_range auto omegaWENO5FMWeights(
//		const std::ranges::common_range auto&& lambda_weights) {
//	/* From Henrick et al.'s mappings of g(λ_k) for the improved
//	 * symmetric normalized lambda-weights of Hong, Ye & Ye
//	 * and linear weights d_k we get the new corrected resultant
//	 * normalized WENO5-FM (WENO5-ZM) omega (ω_k-)weights for WENO5-FM
//	 * (again due to Zheng Hong, Zhengyin Ye and Kun Ye).
//	 */

//	// The ideal weights (they generate the central upstream fifth-order
//	// scheme for the 5-point stencil), which are in WENO usu. called
//	// linear weights:
//	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
//	// From them WENO5-Z and WENO-M will calculate the non-linear
//	// alpha and omega weights.

//	// In WENO5-FM, further, we have one ideal value for λ
//	// \overbar{lambda} = 1/3
//	// T lambda_ideal = 1/3;
//	// In the smooth region the smoothness indicators β_k ought to
//	// be equal for all sub-stencils, and thus the weight's ideal
//	// value must be unique.

//	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
//	// omega_weights = d_lin_weights * alpha_weights;
//	// lambda_weights

//	// And only to the λ-weights a mapping in the spirit of
//	// Henrick et al. is applied:
//	auto gMap = [](double x) -> double {
//		return henrickGMappingForLambda(x);
//	};

//	std::valarray<double> alpha_weights(3);
//	std::ranges::transform(
//				lambda_weights,
//				std::ranges::begin(alpha_weights),
//				gMap);
//	// α*-weights

//	// From α*=g(λ_k) and d_k we get the new corrected resultant
//	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
//	// omega_weights = d_lin_weights * alpha_weights;
//	std::valarray<double> omega_weights = d_lin_weights * alpha_weights;
//	omega_weights /= omega_weights.sum();

//	return omega_weights;
//}


std::ranges::common_range auto omegaWENOFMWeights(
		const std::ranges::common_range auto&& lambda_weights,
		const std::valarray<double>& d_ideal_lin_weights,
		const std::vector<double>& discrete_lambda
					/*= DISCRETE_LAMBDA5*/) {
	/* From Henrick et al.'s mappings of g(λ_k) for the improved
	 * symmetric normalized lambda-weights of Hong, Ye & Ye
	 * and linear weights d_k we get the new corrected resultant
	 * normalized WENO5-FM (WENO5-ZM) omega (ω_k-)weights for WENO5-FM
	 * (again due to Zheng Hong, Zhengyin Ye and Kun Ye).
	 */

	// From them WENO5-Z and WENO-M will calculate the non-linear
	// alpha and omega weights.

	// In WENO5-FM, further, we have one ideal value for λ
	// \overbar{lambda} = 1/3
	// T lambda_ideal = 1/3;
	// In the smooth region the smoothness indicators β_k ought to
	// be equal for all sub-stencils, and thus the weight's ideal
	// value must be unique.

	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	// lambda_weights

	// And only to the λ-weights a mapping in the spirit of
	// Henrick et al. is applied:
//	auto gMap = [](T x) -> T {
//		return henrickGMappingForLambda(x);
//	};

//	std::valarray<T> alpha_weights(3);
//	std::ranges::transform(
//				lambda_weights,
//				std::ranges::begin(alpha_weights),
//				gMap);
	const std::size_t n = std::ranges::size(discrete_lambda) - 1;
	std::valarray<double> alpha_weights(
				std::ranges::size(d_ideal_lin_weights));
	std::ranges::transform(
				lambda_weights,
				std::ranges::begin(alpha_weights),
				[n, &discrete_lambda](double lambda) -> double {
		return discrete_lambda[
				static_cast<std::size_t>(
					static_cast<double>(n) * lambda)];
	});
	// α*-weights

	// From α*=g(λ_k) and d_k we get the new corrected resultant
	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	std::valarray<double> omega_weights = d_ideal_lin_weights
			* alpha_weights;
	omega_weights /= omega_weights.sum();

	return omega_weights;
}


std::ranges::common_range auto F1DWENO5Reconstruction::omegaWENO5FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};

	return omegaWENOFMWeights(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA5);
}


std::ranges::common_range auto F1DWENO5Reconstruction::omegaWENO5ZQMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil), which are in WENO usu. called
	// (optimal) linear weights:
	std::valarray<double> d_lin_weights = {0.98, 0.01, 0.01};

	return omegaWENOFMWeights(
				std::move(lambda_weights), d_lin_weights,
				DISCRETE_LAMBDA5);
}


//std::ranges::common_range auto omegaWENO7FMWeights(
//		const std::ranges::common_range auto&& lambda_weights) {
//	// The ideal weights (they generate the central upstream seventh-order
//	// scheme for the 7-point stencil), which are in WENO usu. called
//	// (optimal) linear weights:
//	// std::valarray<T> d_lin_weights = {4./35., 18./35., 12./35., 1./35.};
//	std::valarray<double> d_lin_weights = {1./35., 12./35., 18./35., 4./35.};

//	return omegaWENOFMWeights(
//				std::move(lambda_weights), d_lin_weights,
//				DISCRETE_LAMBDA7);
//}


double computeFHatWENO5JSReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		double eps = 1e-40, double p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);
//	beta_IS_coefs[0] = (13.0/12.0)*std::pow(
//				f_stencil[0]-2.0*f_stencil[1]+f_stencil[2], 2)
//			+ 0.25*std::pow(
//				f_stencil[0]-4.0*f_stencil[1]+3.0*f_stencil[2], 2);
//	beta_IS_coefs[1] = (13.0/12.0)*std::pow(
//				f_stencil[1]-2.0*f_stencil[2]+f_stencil[3], 2)
//			+ 0.25*std::pow(f_stencil[1]-f_stencil[3], 2);
//	beta_IS_coefs[2] = (13.0/12.0)*std::pow(
//				f_stencil[2]-2.0*f_stencil[3]+f_stencil[4], 2)
//			+ 0.25*std::pow(
//				3.0*f_stencil[2]-4.0*f_stencil[3]+f_stencil[4], 2);

	// non-linear non-scaled (α-)weights
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<double> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
//	std::valarray<double> alpha_weights(3);
//	alpha_weights[0] = 1.0e-1/std::pow((eps+beta_IS_coefs[0]), 2);
//	alpha_weights[1] = 6.0e-1/std::pow((eps+beta_IS_coefs[1]), 2);
//	alpha_weights[2] = 3.0e-1/std::pow((eps+beta_IS_coefs[2]), 2);

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<double> omega_weights = alpha_weights
			/ alpha_weights.sum();
//	std::valarray<double> omega_weights(3);
//	omega_weights[0] = alpha_weights[0]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[1] = alpha_weights[1]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);
//	omega_weights[2] = alpha_weights[2]
//			/ (alpha_weights[0] + alpha_weights[1] + alpha_weights[2]);

	// vecMatDot<double>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<double> eno_reconstructed_f
			= f3OrdReconstructionFromStencil(f_stencil);

//	 eno_reconstructed_f[0] = f_stencil[0]/3.0
//			 - 7.0/6.0*f_stencil[1] + 11.0/6.0*f_stencil[2];
//	 eno_reconstructed_f[1] =-f_stencil[1]/6.0
//			 + 5.0/6.0*f_stencil[2] + f_stencil[3]/3.0;
//	 eno_reconstructed_f[2] = f_stencil[2]/3.0
//			 + 5.0/6.0*f_stencil[3] - f_stencil[4]/6.0;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


double computeFHatWENO5JSReconstructionKernelRev(
		const std::ranges::sized_range auto&& f_stencil,
		double eps = 1e-40, double p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<double> d_lin_weights = {0.3, 0.6, 0.1};
	std::valarray<double> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<double> omega_weights = alpha_weights
			/ alpha_weights.sum();

	// vecMatDot<double>(u_..., WmN...) stores a 3-rd order estimate
	// of f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<double> eno_reconstructed_f(3);
	eno_reconstructed_f[0] = (-1. * f_stencil[0]
						   + 5. * f_stencil[1]
						   + 2. * f_stencil[2]) / 6.;
	eno_reconstructed_f[1] = (2. * f_stencil[1]
						   + 5. * f_stencil[2]
						   - 1. * f_stencil[3]) / 6.;
	eno_reconstructed_f[2] = (11. * f_stencil[2]
						   - 7. * f_stencil[3]
						   + 2. * f_stencil[4]) / 6.;

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];


	return f_hat;
}


double computeFHatWENO5MReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		double eps = 1e-40, double p = 2.) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

	// non-linear non-scaled (α-)weights
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
	std::valarray<double> alpha_weights = d_lin_weights
			/ std::pow(eps + beta_IS_coefs, p);
	// — we have no need of them in WENO-FM!

	// scaled (normalized) non-linear (ω-)weights (ENO weights)
	std::valarray<double> omega_weights = alpha_weights
			/ alpha_weights.sum();

	// vecMatDot<double>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::transform(
				/*std::execution::par_unseq,*/
				std::ranges::begin(omega_weights),
				std::ranges::end(omega_weights),
				std::ranges::begin(d_lin_weights),
				std::ranges::begin(omega_weights),
				[](auto w, auto d) {
		return henrickGMappingForLambda(w, d);
	});  // we obtain a new non-normalized weight alpha*

	omega_weights = omega_weights
			/ omega_weights.sum(); // normalize it

	std::valarray<double> eno_reconstructed_f
			= f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


double F1DWENO5Reconstruction::computeFHatWENO5FMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		double eps/* = 1e-40*/, double p/* = 2.*/) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// (and nothing more in WENO and WENO-(F)M);
	// but in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

	std::array<double, 3> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[eps, p](auto beta) {
		return alphaWENO5FMWeight(beta, eps, p);
	});

	std::array<double, 3> lambda_weights;
	lambdaWENO5FMWeights(std::move(alpha_weights), lambda_weights);

	std::valarray<double> omega_weights = omegaWENO5FMWeights(
				std::move(lambda_weights));

	// vecMatDot<double>(u_..., WmN...) stores a 3-rd order estimate of
	// f_{i+1/2} via linear combinations with WmNplus coefficients
	// for each substencil which is then used to calculate
	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
	// using the nonlinear weights [ω]
//	f_hat = std::inner_product(
//		std::ranges::begin(omega_weights), std::ranges::end(omega_weights),
//		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)), 0.
//	);

	std::valarray<double> eno_reconstructed_f
			= f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


double F1DWENO5Reconstruction::computeFHatWENO5ZMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		double eps/* = 1e-40*/, double p/* = 2.*/) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 *
	 * Borges et al.'s WENO5-Z. Then in improves this by the symmetric
	 * mapping of Hong et al.
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// and in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)). For WENO7-Z(M)
	// p = 4.
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicators(f_stencil);
	double tau_5 = std::abs(beta_IS_coefs[2] - beta_IS_coefs[0]);

	std::array<double, 3> alpha_weights;
	std::ranges::transform(
				beta_IS_coefs,
				std::ranges::begin(alpha_weights),
				[tau_5, eps, p](auto beta) {
		return alphaWENO5ZMWeight(beta, tau_5, eps, p);
	});

	std::array<double, 3> lambda_weights;
	lambdaWENO5FMWeights(std::move(alpha_weights), lambda_weights);

	std::valarray<double> omega_weights = omegaWENO5FMWeights(
				std::move(lambda_weights));

	std::valarray<double> eno_reconstructed_f
			= f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


double F1DWENO5Reconstruction::computeFHatWENO5ZQMReconstructionKernel(
		const std::ranges::sized_range auto&& f_stencil,
		double eps/* = 1e-40*/, double p/* = 2.*/) {
	/* Calculate (reconstruct) one of the two split monotone numerical
	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
	 *                      ^    ^    ^    ^    ^    ^
	 *                      0    1    2    3    4    |
	 * in either case for convenience).
	 *
	 * I.e. this function implements the upwind reconstruction which
	 * should be used for positive fluxes (with information propagated
	 * from left to right) if the nodes are passed in order. However,
	 * the downwind reconstruction should obviously look the same
	 * modulo flipping the values with respect to j+0, so that it
	 * becomes downwind biased and takes one extra point to the right
	 * instead of taking one extra to the left. In other words, to get
	 * this to behave as a downwind reconstrution we need to pass
	 * the points symmetric to those of upwind reconstruction with
	 * respect to j+0:
	 * [j+3, j+2, j+1, j+0, j-1, ...]. (We reverse the points in
	 *      [j-2, j-1, j+0, j+1, j+2] j+3 and get
	 *                  |
	 * [j+3, j+2, j+1, j+0, j-1] j-2.)
	 *
	 * Borges et al.'s WENO5-Z. Then in improves this by the symmetric
	 * mapping of Hong et al.
	 */

	// `p` controls (increases) the amount of numerical dissipation
	// and in WENO-Z(M) changing the value of p alters convergence
	// rates at critical points (it's recommended to take it = r-1
	// for 2r-1 order schemes, so 2 for WENO5-Z(M)). For WENO7-Z(M)
	// p = 4.
	//
	// `eps` is a small positive parameter to avoid the denominator
	// of weights being zero
	// (though it, too, can be significant for convergence properties
	// and should ideally be tailored to the specific comp. problem,
	// as first noted and more or less fully outlined by Henrick et al.)

	// f_stencil = f_plus (or a reversed stencil for f_minus)

	std::valarray<double> beta_IS_coefs(3);

	double f_hat = 0.;

	// smoothness indicators of the stencil
	// (measure how smooth u is in the stencil)
//	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

	// The non-matrix variant seems to be faster(?)
	beta_IS_coefs = betaSmoothnessIndicatorsZQ(f_stencil);
	double tau_5 = std::pow(std::abs(beta_IS_coefs[0] - beta_IS_coefs[1])
			+ std::abs(beta_IS_coefs[0] - beta_IS_coefs[2]), 2) * 0.25;

	std::valarray<double> d_lin_weights = {0.98, 0.01, 0.01};
	std::array<double, 3> alpha_weights; p = 1.;
	std::transform(
				std::ranges::begin(beta_IS_coefs),
				std::ranges::end(beta_IS_coefs),
				std::ranges::begin(alpha_weights),
				std::ranges::begin(d_lin_weights),
				[tau_5, eps, p](auto beta, auto gamma) {
		return gamma * alphaWENO5ZMWeight(beta, tau_5, eps, p);
	});

	std::array<double, 3> lambda_weights;
	lambdaWENO5FMWeights(std::move(alpha_weights), lambda_weights);


	std::valarray<double> omega_weights = {
				lambda_weights[0], lambda_weights[1], lambda_weights[2]};
//	std::valarray<double> omega_weights = omegaWENO5ZQMWeights(
//				std::move(lambda_weights));


	std::valarray<double> eno_reconstructed_f
			= f4OrdZQReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * (eno_reconstructed_f[0] * (1./d_lin_weights[0])
				- eno_reconstructed_f[1] * (
					d_lin_weights[1] / d_lin_weights[0])
				- eno_reconstructed_f[2] * (
					d_lin_weights[2] / d_lin_weights[0]))
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


//template <ArithmeticWith<numeric_val> T>
//T computeFHatWENO5FMReconstructionKernelRev(
//											std::span<double, 5> f_stencil,
//											T eps = 1e-40, T p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical
//	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
//	 * (receives 5 values [j-2, j-1, j+0, j+1, j+2, ...] for '+'
//	 *                (or [j+3, j+2, j+1, j+0, j-1, ...] for '-')
//	 *                      ^    ^    ^    ^    ^    ^
//	 *                      0    1    2    3    4    |
//	 * in either case for convenience).
//	 *
//	 * I.e. this function implements the downwind reconstruction which
//	 * should be used for negative fluxes (with information propagated
//	 * from right to left) if the nodes are passed in order.
//	 */

//	std::valarray<double> beta_IS_coefs(3);

//	T f_hat = 0.;

//	// smoothness indicators of the stencil
//	// (measure how smooth u is in the stencil)
////	 beta_IS_coefs = betaSmoothnessIndicatorsMat(f_stencil);

//	// The non-matrix variant seems to be faster(?)
//	// beta_IS_coefs = betaSmoothnessIndicators(f_stencil);

//	T f_prev2 = f_stencil[0];
//	T f_prev1 = f_stencil[1];
//	T f_curr0 = f_stencil[2];
//	T f_next1 = f_stencil[3];
//	T f_next2 = f_stencil[4];
////	T f_prev2 = f_stencil[4];
////	T f_prev1 = f_stencil[3];
////	T f_curr0 = f_stencil[2];
////	T f_next1 = f_stencil[1];
////	T f_next2 = f_stencil[0];

//	beta_IS_coefs[0] = ((13./12.) * std::pow(
//							f_prev2 - 2.*f_prev1 + f_curr0, 2)
//				+ (1./4.) * std::pow(
//							f_prev2 - 4.*f_prev1 + 3.*f_curr0, 2));

//	beta_IS_coefs[1] = ((13./12.) * std::pow(
//							f_prev1 - 2.*f_curr0 + f_next1, 2)
//				+ (1./4.) * std::pow((f_/usr/bin/prev1 - f_next1), 2));

//	beta_IS_coefs[2] = ((13./12.) * std::pow(
//							f_curr0 - 2.*f_next1 + f_next2, 2)
//				+ (1./4.) * std::pow(
//							3.*f_curr0 - 4.*f_next1 + f_next2, 2));

//	// non-linear non-scaled (α-)weights
//	std::valarray<double> d_lin_weights = {0.3, 0.6, 0.1};
//	std::valarray<double> alpha_weights = d_lin_weights
//			/ std::pow(eps + beta_IS_coefs, p);
//	// — we have no need of them in WENO-FM!

//	// Instead we use:
////	std::valarray<double> alpha_weights = alphaWENO5FMWeights(
////		std::move(beta_IS_coefs), eps, p
////	);

//	// scaled (normalized) non-linear (ω-)weights (ENO weights)
//	 std::valarray<double> omega_weights = alpha_weights
//			 / alpha_weights.sum();

//	// FM(ZM)-improved scaled (normalized) symmetric (λ-)weights
//	// due to Zheng Hong, Zhengyin Ye and Kun Ye:
////	 lambda_weights = alpha_weights / alpha_weights.sum();
////	std::valarray<double> lambda_weights = lambdaWENO5FMWeights(
////		std::move(alpha_weights)
////	);

////	std::valarray<double> omega_weights = omegaWENO5FMWeights(
////		std::move(lambda_weights)
////	);

//	// vecMatDot<double>(u_..., WmN...) stores a 3-rd order estimate of
//	// f_{i+1/2} via linear combinations with WmNplus coefficients
//	// for each substencil which is then used to calculate
//	// f_hat = ∑ ω * q = [ω] * (WmN(+/-) * [f])
//	// using the nonlinear weights [ω]
////	f_hat = std::inner_product(
////		std::ranges::begin(omega_weights),
////		std::ranges::end(omega_weights),
////		std::ranges::begin(f3OrdReconstructionFromStencil(f_stencil)),
////		0.
////	);

//	 std::valarray<
//		 double
//	 > eno_reconstructed_f(3);
//	 eno_reconstructed_f[0] = (-1. * f_stencil[0]
//							   + 5. * f_stencil[1]
//							   + 2. * f_stencil[2]) / 6.;

//	 eno_reconstructed_f[1] = (2. * f_stencil[1]
//							   + 5. * f_stencil[2]
//							   - 1. * f_stencil[3]) / 6.;

//	 eno_reconstructed_f[2] = (11. * f_stencil[2]
//							   - 7. * f_stencil[3]
//							   + 2. * f_stencil[4]) / 6.;

//	f_hat = omega_weights[0] * eno_reconstructed_f[0]
//			+ omega_weights[1] * eno_reconstructed_f[1]
//			+ omega_weights[2] * eno_reconstructed_f[2];

//	return f_hat;
//}


//void F1DWENO5Reconstruction::calcComponent_(
//		const std::ranges::common_range auto&& u,
//		std::ranges::common_range auto&& u_plus_rec,
//		std::ranges::common_range auto&& u_minus_rec,
//		std::size_t n_ghost_cells /* = 3*/) {
//	/* Component-wise finite-volume WENO5FM (WENO5-FM) - space
//	 * reconstruction method with the global Lax-Friedrichs (LF) flux
//	 * splitting.
//	 *
//	 * Usually, componentwise reconstruction produces satisfactory
//	 * results for schemes up to third-order accuracy, while characteristic
//	 * reconstruction produces better nonoscillatory results for
//	 * higher-order accuracy, albeit with an increased computational cost.
//	 */

//	const unsigned order = 5;
//	const std::size_t stencil_size = order;
//	// const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
//	const std::size_t half_size = order / 2;  // 2

//	// r = (order + 1) / 2 = 3
//	assert(n_ghost_cells >= 3);
//	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
//	const std::size_t mini = n_ghost_cells;  // at least 3
//	// const std::size_t maxi = n_ghost_cells + n_size - 1;
//	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
//	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
//	// [g      g      g      i      i      i      i      i      i      ...]
//	// {0      1      2      3      4      5}     6      7      8      ...
//	//  |             |      |      |
//	// itr            j nGhostCells end()

//	// WENO5 stencils

//	// Coefficients WmN(+/-) before fluxes at the stencil nodes
//	// to find component stencils
//	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
//	// of the WENO interpolator.
//	// So 3rd order approximation coefficients associated with
//	// the substencils.

//	// Calculation of f_hat, the numerical flux of u (whichever
//	// is chosen), requires the approximation of u that uses at
//	// the j-th cell of u-discretization a group of cell average
//	// values on the left (`u_minus`) and on the right (`u_plus`).
//	// So left- and right-biased approximations respectively.
//	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
//	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
//	// for convenience and uniformity we represent both using the
//	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
//	// std::valarray<double> u_plus(_actual_stencil_size);   // f_plus
//	// std::valarray<double> u_minus(_actual_stencil_size);  // f_minus

//	// For the purpose of linear stability (upwinding),
//	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
//	// dfminus/du <= 0), is performed.
//	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
//	// and smooth - we need the positive and negative fluxes to have
//	// as many derivatives as the order of our finite-difference WENO)
//	// is chosen here. `f_plus` uses a biased stencil with 1 point to
//	// the left, while `f_minus` uses a biased stencil with 1 point to
//	// the right.

//	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
//	// a monotone numerical flux consistent with the physical one
//	// (f_hat(u, u) = f(u)), will be construced at the end of
//	// the loop below from `fhatminus` and `fhatplus` using
//	// `f_minus` and `f_plus`.
//	// N.B.! This can be replaced by an exact or approximate
//	// Riemann solver (see Toro, 2009). Somehow...
//	// Not every monotone flux can be writtenin the flux split form.
//	// For example, the Godunov flux cannot.
//	// std::valarray<double> numerical_flux(0., u.size());
//	// f = std::valarray<double>(u.size());

//	auto j_it_p = std::ranges::begin(u);  // f_plus

////	std::ranges::transform(
////			monotone_flux_components[0],
////			monotone_flux_components[1],
////			std::ranges::begin(f), [](const auto fp, const auto fm) {

////	})
//	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
//	auto u_plus = std::ranges::views::counted(j_it_p, 6);
//	auto u_minus = u_plus | std::ranges::views::reverse;

//	for (std::size_t j : shifted_index_range) {
//		j_it_p = std::ranges::begin(u);  // f_plus
//		std::advance(j_it_p, j + half_size + 1 - stencil_size);
//		u_plus = std::ranges::views::counted(j_it_p, 6);

//		u_plus_rec[j] = computeFHatWENO5FMReconstructionKernel(
//			std::ranges::views::counted(
//						std::ranges::begin(u_plus), 5), eps, p
//		);


//		u_minus = u_plus | std::ranges::views::reverse;
//		u_minus_rec[j] = computeFHatWENO5FMReconstructionKernel(
//			std::ranges::subrange(
//				std::ranges::begin(u_minus),
//						std::ranges::end(u_minus) - 1), eps, p
//		);
////		u_minus_rec[j] = computeFHatWENO5JSReconstructionKernelRev(
////			std::ranges::views::counted(
////						std::ranges::begin(u_minus)+1, 5), eps, p
////		);
//	}

//	// std::cout << " done!" << "\n";
//}


//double computeFHatWENO7FMReconstructionKernel(
//		const std::ranges::sized_range auto&& f_stencil,
//		double eps = 1e-40, double p = 2.) {
//	/* Calculate (reconstruct) one of the two split monotone numerical
//	 * fluxes `fhatplus`/`fhatminus` at a point j+0 for a given stencil
//	 * (receives the following 7 values
//	 *     [j-3, j-2, j-1, j+0, j+1, j+2, j+3, ...] for '+'
//	 * (or [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...] for '-')
//	 *       ^    ^    ^    ^    ^    ^    ^    ^
//	 *       0    1    2    3    4    5    6    |
//	 * in either case for convenience).
//	 *
//	 * I.e. this function implements the upwind reconstruction which
//	 * should be used for positive fluxes (with information propagated
//	 * from left to right) if the nodes are passed in order. However,
//	 * the downwind reconstruction should obviously look the same
//	 * modulo flipping the values with respect to j+0, so that it
//	 * becomes downwind biased and takes one extra point to the right
//	 * instead of taking one extra to the left. In other words, to get
//	 * this to behave as a downwind reconstrution we need to pass
//	 * the points symmetric to those of upwind reconstruction with
//	 * respect to j+0:
//	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2, ...]. (We reverse the points in
//	 *      [j-3, j-2, j-1, j+0, j+1, j+2, j+3] j+4 and get
//	 *                       |
//	 * [j+4, j+3, j+2, j+1, j+0, j-1, j-2] j-3.)
//	 */

//	// `p` controls (increases) the amount of numerical dissipation
//	// (and nothing more in WENO and WENO-(F)M);
//	// but in WENO-Z(M) changing the value of p alters convergence
//	// rates at critical points (it's recommended to take it = r-1
//	// for 2r-1 order schemes, so 2 for WENO5-Z(M)).
//	//
//	// `eps` is a small positive parameter to avoid the denominator
//	// of weights being zero
//	// (though it, too, can be significant for convergence properties
//	// and should ideally be tailored to the specific comp. problem,
//	// as first noted and more or less fully outlined by Henrick et al.)

//	// f_stencil = f_plus (or a reversed stencil for f_minus)

//	std::valarray<double> beta_IS_coefs(4);

//	double f_hat = 0.;
//	beta_IS_coefs = betaSmoothnessIndicatorsWENO7BS<T>(f_stencil);

//	std::array<double, 4> alpha_weights;
//	std::ranges::transform(
//				beta_IS_coefs,
//				std::ranges::begin(alpha_weights),
//				[eps, p](auto beta) {
//		return alphaWENO5FMWeight(beta, eps, p);
//	});

//	std::array<double, 4> lambda_weights;
//	lambdaWENO5FMWeights(std::move(alpha_weights), lambda_weights);

//	std::valarray<double> omega_weights = omegaWENO7FMWeights(
//				std::move(lambda_weights));

//	std::valarray<double> eno_reconstructed_f
//			= f4OrdReconstructionFromStencil(f_stencil);

//	f_hat = omega_weights[0] * eno_reconstructed_f[0]
//			+ omega_weights[1] * eno_reconstructed_f[1]
//			+ omega_weights[2] * eno_reconstructed_f[2]
//			+ omega_weights[3] * eno_reconstructed_f[3];

//	return f_hat;
//}


void calcHydroStageFVWENO5(
		const std::ranges::common_range auto&& u,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		auto&& computeWENOReconstructionKernel,
		std::size_t n_ghost_cells = 3,
		double eps = 1e-40,
		double p = 2.) {
	/* Component-wise finite-volume WENO5FM (FV WENO5-FM) - space
	 * reconstruction method with the global Lax-Friedrichs (LF) flux
	 * splitting.
	 *
	 * Usually, componentwise reconstruction produces satisfactory
	 * results for schemes up to third-order accuracy, while characteristic
	 * reconstruction produces better nonoscillatory results for
	 * higher-order accuracy, albeit with an increased computational cost.
	 */

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	// const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{mini - 1, maxi + 1};
	// [g      g      g      i      i      i      i      i      i      ...]
	// {0      1      2      3      4      5}     6      7      8      ...
	//  |             |      |      |
	// itr            j nGhostCells end()

	// WENO5 stencils

	// Coefficients WmN(+/-) before fluxes at the stencil nodes
	// to find component stencils
	// [q_k] = WmN(+/-) * [ f[j-2] ... f[j+2] ]
	// of the WENO interpolator.
	// So 3rd order approximation coefficients associated with
	// the substencils.

	// Calculation of f_hat, the numerical flux of u (whichever
	// is chosen), requires the approximation of u that uses at
	// the j-th cell of u-discretization a group of cell average
	// values on the left (`u_minus`) and on the right (`u_plus`).
	// So left- and right-biased approximations respectively.
	// `u_plus` represents the cells [j-2, j-1, j, j+1, j+2],
	// and `u_minus` represents the cells [j-1, j, j+1, j+2, j+3];
	// for convenience and uniformity we represent both using the
	// same combined structure of    [j-2, j-1, j, j+1, j+2, j+3].
	// std::valarray<double> u_plus(_actual_stencil_size);   // f_plus
	// std::valarray<double> u_minus(_actual_stencil_size);  // f_minus

	// For the purpose of linear stability (upwinding),
	// a flux splitting, f = fplus + fminus (dfplus/du >= 0 and
	// dfminus/du <= 0), is performed.
	// Lax-Friedrichs (LF) flux splitting (because it is the simplest
	// and smooth - we need the positive and negative fluxes to have
	// as many derivatives as the order of our finite-difference WENO)
	// is chosen here. `f_plus` uses a biased stencil with 1 point to
	// the left, while `f_minus` uses a biased stencil with 1 point to
	// the right.

	// So an LF flux	`numerical_flux`, f_hat(u_minus, u_plus),
	// a monotone numerical flux consistent with the physical one
	// (f_hat(u, u) = f(u)), will be construced at the end of
	// the loop below from `fhatminus` and `fhatplus` using
	// `f_minus` and `f_plus`.
	// N.B.! This can be replaced by an exact or approximate
	// Riemann solver (see Toro, 2009). Somehow...
	// Not every monotone flux can be writtenin the flux split form.
	// For example, the Godunov flux cannot.
	// std::valarray<double> numerical_flux(0., u.size());
	// f = std::valarray<double>(u.size());

	auto j_it_p = std::ranges::begin(u);  // f_plus

//	std::ranges::transform(
//			monotone_flux_components[0],
//			monotone_flux_components[1],
//			std::ranges::begin(f), [](const auto fp, const auto fm) {

//	})
	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = u_plus | std::ranges::views::reverse;

	for (std::size_t j : shifted_index_range) {
		j_it_p = std::ranges::begin(u);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		u_plus = std::ranges::views::counted(j_it_p, 6);

		u_plus_rec[j] = computeWENOReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);


		u_minus = u_plus | std::ranges::views::reverse;
		u_minus_rec[j] = computeWENOReconstructionKernel(
			std::ranges::subrange(
				std::ranges::begin(u_minus),
						std::ranges::end(u_minus) - 1), eps, p
		);
//		u_minus_rec[j] = computeFHatWENO5JSReconstructionKernelRev(
//			std::ranges::views::counted(
//						std::ranges::begin(u_minus)+1, 5), eps, p
//		);
	}

	// std::cout << " done!" << "\n";
}


//void calcHydroStageFVWENO5JS(
//		const std::ranges::common_range auto&& u,
//		std::ranges::common_range auto&& u_plus_rec,
//		std::ranges::common_range auto&& u_minus_rec,
//		std::size_t n_ghost_cells = 3,
//		double eps = 1e-40,
//		double p = 2.) {
//	calcHydroStageFVWENO5(
//				std::ranges::views::all(u),
//				std::ranges::views::all(u_plus_rec),
//				std::ranges::views::all(u_minus_rec),
//				[](const std::ranges::sized_range auto&& stencil,
//							double eps, double p) -> double {
//					return computeFHatWENO5JSReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
//				},
//				n_ghost_cells,
//				eps,
//				p);
//}


//void calcHydroStageFVWENO5M(
//		const std::ranges::common_range auto&& u,
//		std::ranges::common_range auto&& u_plus_rec,
//		std::ranges::common_range auto&& u_minus_rec,
//		std::size_t n_ghost_cells = 3,
//		double eps = 1e-40,
//		double p = 2.) {
//	calcHydroStageFVWENO5(
//				std::ranges::views::all(u),
//				std::ranges::views::all(u_plus_rec),
//				std::ranges::views::all(u_minus_rec),
//				[](const std::ranges::sized_range auto&& stencil,
//							double eps, double p) -> double {
//					return computeFHatWENO5MReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
//				},
//				n_ghost_cells,
//				eps,
//				p);
//}


//void calcHydroStageFVWENO5ZM(
//		const std::ranges::common_range auto&& u,
//		double t,
//		std::ranges::common_range auto&& u_plus_rec,
//		std::ranges::common_range auto&& u_minus_rec,
//		std::size_t n_ghost_cells = 3,
//		double eps = 1e-40,
//		double p = 2.) {
//	calcHydroStageFVWENO5(
//				std::ranges::views::all(u), t,
//				std::ranges::views::all(u_plus_rec),
//				std::ranges::views::all(u_minus_rec),
//				[](const std::ranges::sized_range auto&& stencil,
//							double eps, double p) -> double {
//					return computeFHatWENO5ZMReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
//				},
//				n_ghost_cells,
//				eps,
//				p);
//}


void F1DWENO5Reconstruction::calcComponent_(
		const std::ranges::common_range auto&& u,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells /* = 3*/) {
	calcHydroStageFVWENO5(
				std::ranges::views::all(u),
				std::ranges::views::all(u_plus_rec),
				std::ranges::views::all(u_minus_rec),
				[&](const std::ranges::sized_range auto&& stencil,
							double eps, double p) -> double {
//					return computeFHatWENO5JSReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
//					return computeFHatWENO5MReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
					return computeFHatWENO5ZMReconstructionKernel(
								std::ranges::views::all(stencil), eps, p);
//					return computeFHatWENO5FMReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);

//					return computeFHatWENO5ZQMReconstructionKernel(
//								std::ranges::views::all(stencil), eps, p);
				},
				n_ghost_cells,
				eps,
				p);
}


void F1DWENO5Reconstruction::calc(C1DField& fld) /*override*/ {
	/* Perform finite-volume component-wise weighted essentially
	 * non-oscillatory 5-th order (FV c'mp't-wise WENO5-JS) reconstruction
	 * for a 1-D Euler equation field, i. e. for any field 3-component field
	 * `C1DField`.
	 */

	auto&& U = std::views::all(fld.U);
	auto&& u_plus = std::views::all(URx);
	auto&& u_minus = std::views::all(ULx) | std::ranges::views::drop(1);
	auto&& components = {
		&Vector4::x
		, &Vector4::y
		, &Vector4::z
		// , &Vector4::w
	};
	std::for_each(
//			std::execution::par_unseq,
			std::ranges::begin(components),
			std::ranges::end(components),
			[&](auto&& kth_vector_component) {
		calcComponent_(
			U       | std::ranges::views::transform(kth_vector_component),
			u_plus  | std::ranges::views::transform(kth_vector_component),
			u_minus | std::ranges::views::transform(kth_vector_component),
			fld.imin
			// n_size = fld.imax - fld.imin + 1
			);
	});
}


void F1DCharWiseWENO5Reconstruction::calc_(
		const std::ranges::common_range auto&& u,
		const std::ranges::common_range auto&& q_avg,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		auto& project,
		std::size_t n_ghost_cells /* = 3*/) {

	const unsigned order = 5;
	const std::size_t stencil_size = order;
	const std::size_t _actual_stencil_size = stencil_size + 1;  // 6
	const std::size_t half_size = order / 2;  // 2

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::common_view(
				std::views::iota(mini - 1)
					| std::views::take(maxi + 1 - (mini - 1) + 1));

	auto j_it_p = std::ranges::begin(u);  // f_plus

	std::advance(j_it_p, mini - 1 + half_size + 1 - stencil_size);
	auto u_plus = std::ranges::views::counted(j_it_p, 6);
	auto u_minus = u_plus | std::ranges::views::reverse;

	auto components = {
		std::make_pair(0, &Vector4::x),
		std::make_pair(1, &Vector4::y),
		std::make_pair(2, &Vector4::z)
//		&Vector4<T>::w
	};

	std::for_each(
				//std::execution::par_unseq,
				std::ranges::begin(shifted_index_range),
				std::ranges::end(shifted_index_range),
				[&](std::size_t j) {
		j_it_p = std::ranges::begin(u);  // f_plus
		std::advance(j_it_p, j + half_size + 1 - stencil_size);
		// if ()  // maybe try a discontinuity detector/indicator
		auto proj_u_j = [&](auto u) -> decltype(u) {
			const auto q = q_avg[j];
			return project(q, u);
		};

		for (auto& comp : components) {
			u_plus = std::ranges::views::counted(j_it_p, 6);
			u_plus_rec[j][comp.first] = computeFHatWENO5FMReconstructionKernel(
				std::ranges::views::counted(
							std::ranges::begin(u_plus), 5)
						| std::ranges::views::transform(proj_u_j)
						| std::ranges::views::transform(
							comp.second), eps, p
			);

			u_minus = u_plus | std::ranges::views::reverse;
			u_minus_rec[j][comp.first] = computeFHatWENO5FMReconstructionKernel(
				std::ranges::views::counted(
							std::ranges::begin(u_minus), 5)
						| std::ranges::views::transform(proj_u_j)
						| std::ranges::views::transform(
							comp.second), eps, p
			);
		}
	});
}


/*Eigen::Matrix<double, 3, 3>*/Matrix4 EigenLeft1DEulerEigenMatrix(
		Vector4 vec, FEOS& eos) {
	double u = vec[1];
	if (vec[0] != 0.)
		u /= vec[0];

//	double phi_square = 0.5 * gamma_m * u * u;
	double e = vec[2]/vec[0] - .5*u*u;
	double p = eos.getp(vec[0],e);
	double c = eos.getc(vec[0], p);

	double c_s_square = c * c;
	double c_s = std::abs(std::sqrt(c_s_square));

	double uc = u * c_s;
	double h = (vec[3] + p) / vec[0];
//	if (vec[0] != 0.)
//		h = (vec[3] + p) / vec[0];

//	Eigen::Matrix<T, 3, 3> l_mat {
//		{1. - phi_square / c_s_square,
//					gamma_m * u / c_s_square, -gamma_m / c_s_square},
//		{phi_square - uc, +c_s - gamma_m * u, gamma_m              },
//		{phi_square + uc, -c_s - gamma_m * u, gamma_m              }
//	};
//	Eigen::Matrix<double, 3, 3> l_mat {
//		{1.,      1.,          1.     },
//		{u - c_s, u + 0.,      u + c_s},
//		{H - uc,  0.5 * u * u, H + uc }
//	};

	double b = 0.;
	if (vec[0] != 0.)
		b = eos.getdpde(vec[0], e) / vec[0];

//	Eigen::Matrix<double, 3, 3> l_mat {
//		{1.,                           1.,        1.},
//		{u - c_s,                  u + 0.,   u + c_s},
//		{h -  uc,      h - c_s_square / b,   h +  uc},
//	};
	Matrix4 l_mat (
		1.,                           1.,        1.,	0.,
		u - c_s,                  u + 0.,   u + c_s,	0.,
		h -  uc,      h - c_s_square / b,   h +  uc,	0.,
		0.,	0.,	0.,	0.
	);

	return l_mat;
}


/*Eigen::Matrix<double, 3, 3>*/Matrix4 EigenRight1DEulerEigenMatrix(
		Vector4 vec, FEOS& eos) {
	double u = vec[1];
	if (vec[0] != 0.)
		u /= vec[0];

//	double phi_square = 0.5 * gamma_m * u * u;
	double e = vec[2]/vec[0] - .5*u*u;
	double p = eos.getp(vec[0],e);
	double c = eos.getc(vec[0], p);

	double c_s_square = c * c;
	double beta = 1.;
	if (c_s_square != 0.)
		beta /= (2. * c_s_square);
//	else
//		beta = 0.;
	double c_s = std::sqrt(c_s_square);

	double uc = u * c_s;
	double h = (vec[3] + p) / vec[0];
//	if (vec[0] != 0.)
//		h = (vec[3] + p) / vec[0];

//	Eigen::Matrix<T, 3, 3> r_mat {
//		{1.,                   beta,             beta            },
//		{u,                    beta * (u + c_s), beta * (u - c_s)},
//		{phi_square / gamma_m, beta * (H + uc),  beta * (H - uc) },
//	};
//	Eigen::Matrix<double, 3, 3> r_mat {
//		{H + c_s * (u - c_s) / gamma_m,       -(u + c_s / gamma_m), 1.},
//		{-2. * H + 4. * c_s_square / gamma_m, 2. * u,              -2.},
//		{H - c_s * (u + c_s) / gamma_m,       -u + c_s / gamma_m,   1.},
//	};
//	r_mat *= gamma_m * beta;

	double dpde = eos.getdpde(vec[0], e);
	double b = 0.;
	if (vec[0] != 0.)
		b = dpde / vec[0];

	// harder to compute for MG
	// double dpdrho = eos.getdpdrho(vec[0], e);
	double dpdrho = c_s_square - p * b / vec[0];

	double theta = u * u
			- vec[3] / vec[0]
			+ vec[0] * dpdrho / dpde;

//	Eigen::Matrix<double, 3, 3> r_mat {
//		{theta  +  uc / b, -(u + c_s / b),    1.},
//		{2. * (h - u * u),         2. * u,   -2.},
//		{theta  -  uc / b,   -u + c_s / b,    1.},
//	};
	Matrix4 r_mat(
		b * beta * (theta  +  uc / b), -(u + c_s / b) * b * beta,    b * beta,	0.,
		b * beta * 2. * (h - u * u),         2. * u * b * beta,   -2. * b * beta,	0.,
		(theta  -  uc / b) * b * beta,   (-u + c_s / b) * b * beta,    1. * b * beta,	0.,
		0.,	0.,	0.,	0.
	);

	return r_mat;
}


Vector4 projectOntoCharacteristics(
		Vector4 conservative_variables, Vector4 vec, FEOS& eos) {
	return Vector4(EigenRight1DEulerEigenMatrix(
				conservative_variables, eos)
			* vec/*Eigen::Matrix<double, 3, 1>{vec[0], vec[1], vec[2]}*/);
}


Vector4 projectCharacteristicVariablesBackOntoConserved(
		Vector4 conservative_variables, Vector4 vec, FEOS& eos) {
	return Vector4(EigenLeft1DEulerEigenMatrix(
				conservative_variables, eos)
			* vec/*Eigen::Matrix<double, 3, 1>{vec[0], vec[1], vec[2]}*/);
}


template <typename T>
T average(T left, T right) {
	/* A simple arithmetic mean average of 2 values. */
	return (left + right) * 0.5;
}


void F1DCharWiseWENO5Reconstruction::calc(C1DField& fld) /*override*/ {

	std::vector<Vector4> avg(std::ranges::size(fld.U));

	auto&& U = std::ranges::views::all(fld.U);

	std::transform(
				//std::execution::par_unseq,
				std::ranges::begin(U) + 1,
				std::ranges::end(U),
				std::ranges::begin(U),
				std::ranges::begin(avg) + 1,
				[](auto q_l, auto q_r) {
		return average<Vector4>(q_l, q_r);
	});

	auto project = [&](
			Vector4 q_ast, Vector4 vec) -> Vector4 {
		return projectOntoCharacteristics(q_ast, vec, eos);
	};

	auto project_back = [&](auto q_ast, auto f) {
		return projectCharacteristicVariablesBackOntoConserved(
					q_ast, f, eos);
	};

	auto&& u_plus = std::ranges::views::all(URx);
	auto&& u_minus = std::ranges::views::all(ULx)
			| std::ranges::views::drop(1);


	calc_(
				std::move(U),
				std::ranges::views::all(avg),
				u_plus,
				u_minus,
				project,
				fld.imin);

	/*u_plus = std::ranges::views::all(URx);*/
	std::transform(
//				std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg) - 1,
				std::ranges::begin(u_plus),
				std::ranges::begin(u_plus),
				project_back);

	//u_minus = std::ranges::views::all(ULx)
	//	| std::ranges::views::drop(1);
	std::transform(
				//std::execution::par_unseq,
				std::ranges::begin(avg),
				std::ranges::end(avg) - 1,
				std::ranges::begin(u_minus),
				std::ranges::begin(u_minus),
				project_back);
}
