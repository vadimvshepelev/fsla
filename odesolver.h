#ifndef ODESOLVER_H
#define ODESOLVER_H

#include <memory>

#include "C1DProblem.h"
#include "FEOS.h"
#include "C1DField.h"
#include "C1DMethod.h"

class ODESolver {
protected:
	C1DProblem& pr;
	FEOS& eos;
	C1DField& fld;
	C1DMethod& mtd;

	static constexpr const std::size_t order = 0;

public:
	ODESolver(
			C1DProblem& _pr,
			FEOS& _eos,
			C1DField& _fld,
			C1DMethod& _mtd)
		: pr(_pr), eos(_eos), fld(_fld), mtd(_mtd) {}
	virtual void solve() = 0;
	virtual constexpr std::size_t get_order() = 0;
};


class ForwardEuler : public ODESolver {
protected:
	static constexpr const std::size_t order = 1;

public:
	ForwardEuler(
			C1DProblem& _pr,
			FEOS& _eos,
			C1DField& _fld,
			C1DMethod& _mtd)
		: ODESolver(_pr, _eos, _fld, _mtd) {}

	void solve() override { mtd.calc(pr, eos, fld); }

	constexpr std::size_t get_order() override { return order; };
};


class SSPERK3_3 : public ODESolver {
	C1DField inter_fld1;
	C1DField inter_fld2;

protected:
	static constexpr const std::size_t order = 3;

public:
	SSPERK3_3(
			C1DProblem& _pr,
			FEOS& _eos,
			C1DField& _fld,
			C1DMethod& _mtd)
		: ODESolver(_pr, _eos, _fld, _mtd),
		  inter_fld1(C1DField(_pr)),
		  inter_fld2(C1DField(_pr)) {}

	void solve() override;

	constexpr std::size_t get_order() override { return order; };
};


class ERK6_5 : public ODESolver {
	C1DField inter_fld1;
	C1DField inter_fld2;
	C1DField inter_fld3;
	C1DField inter_fld4;
	C1DField inter_fld5;

protected:
	static constexpr const std::size_t order = 5;

public:
	ERK6_5(
			C1DProblem& _pr,
			FEOS& _eos,
			C1DField& _fld,
			C1DMethod& _mtd)
		: ODESolver(_pr, _eos, _fld, _mtd),
		  inter_fld1(C1DField(_pr)),
		  inter_fld2(C1DField(_pr)),
		  inter_fld3(C1DField(_pr)),
		  inter_fld4(C1DField(_pr)),
		  inter_fld5(C1DField(_pr)) {}

	void solve() override;

	constexpr std::size_t get_order() override { return order; };
};


template <std::size_t N>
class MultistepSolver : public ODESolver {
public:
	static constexpr const size_t step_number = N;

protected:
	std::array<std::unique_ptr<C1DField>, step_number/* - 1*/> steps;
	static constexpr const std::size_t order = 0;

public:
	ODESolver& starting_solver;

	MultistepSolver(
				C1DProblem& _pr,
				FEOS& _eos,
				C1DField& _fld,
				C1DMethod& _mtd,
				ODESolver& _startup)
			: ODESolver(_pr, _eos, _fld, _mtd),
			  starting_solver(_startup) {
		for (std::size_t k = 0; k < step_number/* - 1*/; ++ k)
			steps[k] = std::make_unique<C1DField>(_pr);
	}

	virtual void solve() override = 0;
	virtual constexpr std::size_t get_order() override = 0;
};


class eBDF5 : public MultistepSolver<5> {
private:
	bool prepared = false;

protected:
	static constexpr const std::size_t order = 5;

public:
	eBDF5(
				C1DProblem& _pr,
				FEOS& _eos,
				C1DField& _fld,
				C1DMethod& _mtd,
				ODESolver& _startup)
			: MultistepSolver(_pr, _eos, _fld, _mtd, _startup) {}

	void solve() override;
	constexpr std::size_t get_order() override { return order; };
};

#endif // ODESOLVER_H
