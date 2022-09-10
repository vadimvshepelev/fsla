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
		// ULx.push_back(tempVect);
		ULx.push_back(Vector4::ZERO);
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


void F1DENO2Reconstruction::calc(C1DField& fld) {
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

	if (std::abs(u_stencil[2] - u_stencil[1])
			<= std::abs(u_stencil[3] - u_stencil[2])) {
		if (std::abs(u_stencil[2] - 2. * u_stencil[1] + u_stencil[0])
				<= std::abs(u_stencil[3] - 2. * u_stencil[2] + u_stencil[1])
				) {
			return 0;
		} else {
			return 1;
		}
	}

	if (std::abs(u_stencil[3] - 2. * u_stencil[2] + u_stencil[1])
			<= std::abs(u_stencil[4] - 2. * u_stencil[3] + u_stencil[2]))
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

	// r = (order + 1) / 2 = 3
	assert(n_ghost_cells >= 3);
	// const std::size_t n_ghost_cells = (stencil_size + 1) / 2;  // 3
	const std::size_t mini = n_ghost_cells;  // at least 3
	// const std::size_t maxi = n_ghost_cells + n_size - 1;
	const std::size_t maxi = std::ranges::size(u) - n_ghost_cells - 1;
	auto shifted_index_range = std::ranges::iota_view{
		mini - 1,
		maxi + 1
	};

	double uhatminus = 0.;
	double uhatplus = 0.;

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


void F1DENO3Reconstruction::calc(C1DField& fld) {
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
			// std::execution::par_unseq,
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
