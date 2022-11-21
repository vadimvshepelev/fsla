#include "odesolver.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <execution>
#include <numeric>
#include <ranges>

void SSPERK3_3::solve() /* override */ {
	/* Optimal 3rd Order 3 Stage Explicit Total Variation Diming
	 * / Diminishing (Strong Stability Preserving)
	 * Runge-Kutta Scheme (TVD RK3 / SSPRK(3,3))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 * This is generally known as the Shu–Osher method.
	 * It has 3 stages. SSP coefficient = 1.
	 * See Shu and Osher (1988).
	 *
	 * (It is important to note that this method is
	 * not only strong stability preserving but
	 * also internally stable.)
	 */

	pr.setbcs(fld.U);
	double dt = fld.dt;
	double dx = fld.dx;
	auto&& U = fld.U;
	auto&& F = fld.F;
	auto&& U1 = inter_fld1.U;
	auto&& U2 = inter_fld2.U;
	auto&& F1 = inter_fld1.F;
	auto&& F2 = inter_fld2.F;
//	const std::size_t n_size = pr.nx;
	const std::size_t full_size = std::ranges::size(U);

	inter_fld1.dt = fld.dt;
	inter_fld2.dt = fld.dt;
	inter_fld1.t = fld.t;
	inter_fld2.t = fld.t;

	auto&& iv = std::ranges::common_view(
				std::ranges::views::iota(std::size_t(0))
					| std::views::take(full_size - 1)
	);

	// ------------------------First Stage----------------------------
	mtd.calcFluxField(pr, eos, fld);  // L1 = L[u^n]
	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U1, &U, &F](std::size_t k) {
		U1[k] = U[k] + dt / dx * (F[k] - F[k+1]);
	});
	// u(1) = u^n + Δt L[u^n]

	pr.setbcs(U1);


	// ------------------------Second Stage---------------------------
	mtd.calcFluxField(pr, eos, inter_fld1);  // L2 = L[u(1)]

	iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(full_size - 1)
	);
	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U2, &U, &F1, &U1](std::size_t k) {
		U2[k] = (3. * U[k] + dt / dx * (F1[k] - F1[k+1]) + U1[k]) * 0.25;
	});
	// u(2) = 0.75 * u^n + 0.25 * u(1) + 0.25 * Δt L[u(1)]

	pr.setbcs(U2);


	// ------------------------Third Stage----------------------------
	mtd.calcFluxField(pr, eos, inter_fld2);  // L3 = L[u(2)]


	iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(full_size - 1)
	);
	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[&U2, &F2, dx, dt](std::size_t k) {
		U2[k] = 2. * (U2[k] + dt / dx * (F2[k]-F2[k+1]));
	});
	std::transform(
				std::execution::par_unseq,
				std::ranges::begin(U/* | interior_view*/),
				std::ranges::end(U/* | interior_view*/),
				std::ranges::begin(U2/* | interior_view*/),
				std::ranges::begin(U/* | interior_view*/),
				[dt](const auto u, const auto y) {
		return (u + y) * (1./3.);
	});
	// u^(n+1) = (1/3) * u^n + (2/3) * u(2) + (2/3) * Δt L[u(2)]

	pr.setbcs(U);
}


void ERK6_5::solve() /* override */ {
	/* 5th Order 6 Stage Explicit Runge-Kutta Scheme (ERK(6,5))
	 * to discretize a method-of-lines (MOL) ODE
	 * du/dt = L[u], where L is some spatial operator.
	 *
	 * This is an embedded RK(6, 5) / SSP(3, 3) method introduced
	 * in Colin Barr Macdonald's thesis:
	 * see 'Constructing High-Order Runge-Kutta Methods with
	 * Embedded Strong-Stability-Preserving Pairs' (2001)
	 * and successfully employed by Henrick et al. with WENO5-M
	 * to resolve detonation waves.
	 */

	double dt = fld.dt;
	double dx = fld.dx;
	auto&& U = fld.U;
	auto&& F = fld.F;
	auto&& Y0 = U;
	auto&& L0 = F;
	auto&& Y1 = inter_fld1.U; inter_fld1.dt = dt;
	auto&& L1 = inter_fld1.F;
	auto&& Y2 = inter_fld2.U; inter_fld2.dt = dt;
	auto&& L2 = inter_fld2.F;
	auto&& Y3 = inter_fld3.U; inter_fld3.dt = dt;
	auto&& L3 = inter_fld3.F;
	auto&& Y4 = inter_fld4.U; inter_fld4.dt = dt;
	auto&& L4 = inter_fld4.F;
	auto&& Y5 = inter_fld5.U; inter_fld5.dt = dt;
	auto&& L5 = inter_fld5.F;

	// std::slice Nint(3, nSize, 1);
	//	const std::size_t n_size = pr.nx;
	const std::size_t full_size = std::ranges::size(U);


	auto iv = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(0))
				| std::views::take(std::ranges::size(U))
	);

//	auto interior_view = std::views::drop(n_ghost_points)
//			| std::views::take(n_size)
//			| std::ranges::views::common;
//	auto interior_view = std::ranges::views::common;

	// ------------------------First Stage----------------------------
	pr.setbcs(fld.U);
	mtd.calcFluxField(pr, eos, fld);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[&](std::size_t k) {
		Y1[k] = U[k] + (L0[k] - L0[k+1]) * dt / dx;
	});

	pr.setbcs(Y1);


	// ------------------------Second Stage-----------------------------
	mtd.calcFluxField(pr, eos, inter_fld1);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U, &L0, &L1, &Y2](std::size_t k) {
		Y2[k] =  U[k]
				+ (L0[k] - L0[k+1] + L1[k] - L1[k+1]) * 0.25 * dt / dx;
	});

	pr.setbcs(Y2);


	// ------------------------Third Stage------------------------------
	mtd.calcFluxField(pr, eos, inter_fld2);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U, &L0, &L1, &L2, &Y3](std::size_t k) {
		Y3[k] =  U[k] + (2046. * (L0[k] - L0[k+1])
					- 454. * (L1[k] - L1[k+1])
					+ 1533. * (L2[k] - L2[k+1])) * dt / dx / 15625.;
	});

	pr.setbcs(Y3);


	// ------------------------Fourth Stage-----------------------------
	mtd.calcFluxField(pr, eos, inter_fld3);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U, &L0, &L1, &L2, &L3, &Y4](std::size_t k) {
		Y4[k] =  U[k] + (-2217. * (L0[k] - L0[k+1])
					+ 1533. * (L1[k] - L1[k+1])
					- 566. * (L2[k] - L2[k+1])
					+ 12500. * (L3[k] - L3[k+1])) * dt / dx / 16875.;
	});

	pr.setbcs(Y4);


	// ------------------------Fifth Stage------------------------------
	mtd.calcFluxField(pr, eos, inter_fld4);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U, &L0, &L1, &L2, &L3, &L4, &Y5](std::size_t k) {
		Y5[k] =  U[k]
				+ (11822. * (L0[k] - L0[k+1])
					- 6928. * (L1[k] - L1[k+1])
					- 4269. * (L2[k] - L2[k+1])
					- 12500. * (L3[k] - L3[k+1])
					+ 33750. * (L4[k] - L4[k+1])) * dt / dx / 21875.;
	});

	pr.setbcs(Y5);


	// ------------------------Final Stage------------------------------
	mtd.calcFluxField(pr, eos, inter_fld5);

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(iv/* | interior_view*/),
				std::ranges::end(iv/* | interior_view*/),
				[dt, dx, &U, &L0, &L1, &L3, &L4, &L5](std::size_t k) {
		U[k] += (14. * (L0[k] - L0[k+1])
					+ 125. * (L3[k] - L3[k+1])
					+ 162. * (L4[k] - L4[k+1])
					+ 35. * (L5[k] - L5[k+1])) * dt / dx / 336.;
	});

	pr.setbcs(fld.U);
}


void eBDF5::solve() /* override */ {
	/* 5th Order Explicit 5 Step extrapolated
	 * Backward Differentiation Formula (eBDF5).
	 * At least it only has one stage.
	 *
	 * Stable with WENO with CFL <= 0.2.
	 */

	const double dt = fld.dt;
	const double dx = fld.dx;
	std::ranges::common_range auto&& U = fld.U;
	std::ranges::common_range auto&& F = fld.F;
	if (!prepared) {  // Makes more sense as a separate method...
		double t = fld.t;
		fld.dt = std::pow(dt, get_order() / starting_solver.get_order());
		// super small starting time steps
		// ~ dx^(5. / starting_solver.order)

		(steps[(step_number - 1) - 4])->U = U;  // copy
		/*starting_solver.*/mtd.calcFluxField(pr, eos, fld);  // calc. F
		(steps[(step_number - 1) - 4])->F = F;  // copy
		(steps[(step_number - 1) - 4])->t = t;
		(steps[(step_number - 1) - 4])->dt = fld.dt;

		for (auto& prev_step_field : std::ranges::views::all(steps)
				| std::ranges::views::drop(1)) {  // 4 iterations
			fld.t += fld.dt;
			fld.dt = std::pow(/*starting_solver.*/mtd.calcdt(pr, eos, fld),
							  get_order() / starting_solver.get_order());
			starting_solver.solve();  // solving with some 1-step solver

			prev_step_field->U = U;
			prev_step_field->F = F;
			// copy — subject to optimization and major refactoring?
			prev_step_field->dt = fld.dt;
			prev_step_field->t = fld.t;
		}

		fld.t += (steps[(step_number - 1) - 1])->dt;
		fld.dt = mtd.calcdt(pr, eos, fld);

		prepared = true;
	}

	std::vector<Vector4>& Y1
			= (steps[(step_number - 1) - 1])->U;
	std::vector<Vector4>& L1
			= (steps[(step_number - 1) - 1])->F;
	std::vector<Vector4>& Y2
			= (steps[(step_number - 1) - 2])->U;
	std::vector<Vector4>& L2
			= (steps[(step_number - 1) - 2])->F;
	std::vector<Vector4>& Y3
			= (steps[(step_number - 1) - 3])->U;
	std::vector<Vector4>& L3
			= (steps[(step_number - 1) - 3])->F;
	std::vector<Vector4>& Y4
			= (steps[(step_number - 1) - 4])->U;
	std::vector<Vector4>& L4
			= (steps[(step_number - 1) - 4])->F;
	std::vector<Vector4>& Ytemp
			= (steps[(step_number - 1)])->U;
	std::vector<Vector4>& L0
			= (steps[(step_number - 1)])->F;
	Ytemp = U;  // copy (!)

	const std::size_t n_size = pr.nx;
	const std::size_t n_ghost_points = pr.get_order();

	auto interior_view_indices = std::ranges::common_view(
			std::ranges::views::iota(std::size_t(n_ghost_points/*0*/))
				| std::views::take(n_size/*std::ranges::size(U)*/)
	);

	(steps[(step_number - 1)])->dt = dt;
	(steps[(step_number - 1)])->t = fld.t;
	mtd.calcFluxField(pr, eos, *(steps[(step_number - 1)]));  // fill L0

	std::for_each(
				std::execution::par_unseq,
				std::ranges::begin(interior_view_indices),
				std::ranges::end(interior_view_indices),
				[&U, &Y1, &Y2, &Y3, &Y4,
					&L0, &L1, &L2, &L3, &L4, dt, dx](std::size_t k) {
		U[k] = ((300. * (U[k] - Y1[k])
					+ 200. * Y2[k]
					- 75. * Y3[k]
					+ 12. * Y4[k])
				+ dt / dx * 10. * (30. * (L0[k] - L0[k + 1])
					- 60. * (L1[k] - L1[k + 1])
					+ 60. * (L2[k] - L2[k + 1])
					- 30. * (L3[k] - L3[k + 1])
					+ 6. * (L4[k] - L4[k + 1]))) * (1./137.);
	});

	pr.setbcs(U);

	// Shuffle pointers to replace the steps
	steps[(step_number - 1) - 4].swap(
				steps[(step_number - 1) - 3]);  // U4 <- U3, L4 <- L3
	steps[(step_number - 1) - 3].swap(
				steps[(step_number - 1) - 2]);  // U3 <- U2, L3 <- L2
	steps[(step_number - 1) - 2].swap(
				steps[(step_number - 1) - 1]);  // U2 <- U1, L2 <- L1
	steps[(step_number - 1) - 1].swap(
				steps[(step_number - 1) - 0]);  // U1 <- U(old), L1 <- L0
}
