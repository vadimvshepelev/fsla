#ifndef _F1DRECONSTRUCTION_H_
#define _F1DRECONSTRUCTION_H_

#include "_vector4.h"
#include "C1DProblem.h"
#include "C1DField.h"

#include <vector>

using namespace std;

class F1DReconstruction {
public:
	// vector<vector<double>> ULx, URx;
	vector<Vector4> ULx, URx;
	F1DReconstruction(C1DField& fld);
	virtual void calc(C1DField& fld) = 0;
};


class F1DENO2Reconstruction : public F1DReconstruction {
	vector<double> rm1r, r0r, rm1l, r0l;
public:
	F1DENO2Reconstruction(C1DField& fld);
	void calc(C1DField& fld) override;
};


class F1DENO3Reconstruction : public F1DENO2Reconstruction {
public:
	F1DENO3Reconstruction(C1DField& fld);
	void calc(C1DField& fld) override;
private:
	void calcComponent_(const std::ranges::common_range auto&& u,
						std::ranges::common_range auto&& u_plus_rec,
						std::ranges::common_range auto&& u_minus_rec,
						std::size_t n_ghost_cells = 3);
};


class F1DWENO5Reconstruction : public F1DENO3Reconstruction {
public:
	F1DWENO5Reconstruction(C1DField& fld);
	void calc(C1DField& fld) override;
	// std::string type = "JS";

	std::string type = "FM";
	const double eps = 1e-40;
	const double p = 2;
private:
	void calcComponent_(const std::ranges::common_range auto&& u,
						std::ranges::common_range auto&& u_plus_rec,
						std::ranges::common_range auto&& u_minus_rec,
						std::size_t n_ghost_cells = 3);
protected:
	std::vector<double> DISCRETE_LAMBDA5;

	static std::vector<double> prediscretizeWENO5LambdaMapping(
			std::size_t, double);

	std::ranges::common_range auto omegaWENO5FMWeights(
			const std::ranges::common_range auto&& lambda_weights);

	std::ranges::common_range auto omegaWENO5ZQMWeights(
			const std::ranges::common_range auto&& lambda_weights);

	double computeFHatWENO5FMReconstructionKernel(
			const std::ranges::sized_range auto&& f_stencil,
			double eps = 1e-40, double p = 2.);

	double computeFHatWENO5ZMReconstructionKernel(
			const std::ranges::sized_range auto&& f_stencil,
			double eps = 1e-40, double p = 2.);

	double computeFHatWENO5ZQMReconstructionKernel(
			const std::ranges::sized_range auto&& f_stencil,
			double eps = 1e-40, double p = 2.);
};


class F1DCharWiseWENO5Reconstruction : public F1DWENO5Reconstruction {
public:
	F1DCharWiseWENO5Reconstruction(C1DField& fld, FEOS& eos);
	void calc(C1DField& fld) override;
	std::string type = "FM";
	FEOS& eos;
private:
	void calc_(
			const std::ranges::common_range auto&& u,
			const std::ranges::common_range auto&& q_avg,
			std::ranges::common_range auto&& u_plus_rec,
			std::ranges::common_range auto&& u_minus_rec,
			auto& project,
			std::size_t n_ghost_cells = 3);
};


#endif
