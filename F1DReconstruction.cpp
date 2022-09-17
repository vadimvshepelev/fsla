#include"F1DReconstruction.h"

#include <algorithm>
#include <array>
#include <execution>
#include <ranges>
// #include <span>

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
	: F1DENO3Reconstruction(_fld) {}


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
			std::execution::par_unseq,
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
			   * (1. - 2.*square_ideal));
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


std::ranges::common_range auto omegaWENO5FMWeights(
		const std::ranges::common_range auto&& lambda_weights) {
	/* From Henrick et al.'s mappings of g(λ_k) for the improved
	 * symmetric normalized lambda-weights of Hong, Ye & Ye
	 * and linear weights d_k we get the new corrected resultant
	 * normalized WENO5-FM (WENO5-ZM) omega (ω_k-)weights for WENO5-FM
	 * (again due to Zheng Hong, Zhengyin Ye and Kun Ye).
	 */

	// The ideal weights (they generate the central upstream fifth-order
	// scheme for the 5-point stencil), which are in WENO usu. called
	// linear weights:
	std::valarray<double> d_lin_weights = {0.1, 0.6, 0.3};
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
	auto gMap = [](double x) -> double {
		return henrickGMappingForLambda(x);
	};

	std::valarray<double> alpha_weights(3);
	std::ranges::transform(
				lambda_weights,
				std::ranges::begin(alpha_weights),
				gMap);
	// α*-weights

	// From α*=g(λ_k) and d_k we get the new corrected resultant
	// normalized WENO5-FM (WENO5-ZM) (ω_k-)weights:
	// omega_weights = d_lin_weights * alpha_weights;
	std::valarray<double> omega_weights = d_lin_weights * alpha_weights;
	omega_weights /= omega_weights.sum();

	return omega_weights;
}


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

	omega_weights = omega_weights / omega_weights.sum(); // normalize it

	std::valarray<double> eno_reconstructed_f
			= f3OrdReconstructionFromStencil(f_stencil);

	f_hat = omega_weights[0] * eno_reconstructed_f[0]
			+ omega_weights[1] * eno_reconstructed_f[1]
			+ omega_weights[2] * eno_reconstructed_f[2];

	return f_hat;
}


double computeFHatWENO5FMReconstructionKernel(
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
	beta_IS_coefs = betaSmoothnessIndicators<double>(f_stencil);

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


void F1DWENO5Reconstruction::calcComponent_(
		const std::ranges::common_range auto&& u,
		std::ranges::common_range auto&& u_plus_rec,
		std::ranges::common_range auto&& u_minus_rec,
		std::size_t n_ghost_cells /* = 3*/) {
	/* Component-wise finite-volume WENO5FM (WENO5-FM) - space
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

		u_plus_rec[j] = computeFHatWENO5JSReconstructionKernel(
			std::ranges::views::counted(
						std::ranges::begin(u_plus), 5), eps, p
		);


		u_minus = u_plus | std::ranges::views::reverse;
		u_minus_rec[j] = computeFHatWENO5JSReconstructionKernel(
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
			std::execution::par_unseq,
			std::ranges::begin(components),
			std::ranges::end(components),
			[&](auto&& kth_vector_component) {
		F1DWENO5Reconstruction::calcComponent_(
			U       | std::ranges::views::transform(kth_vector_component),
			u_plus  | std::ranges::views::transform(kth_vector_component),
			u_minus | std::ranges::views::transform(kth_vector_component),
			3
			// n_size = fld.imax - fld.imin + 1
			);
	});
}

