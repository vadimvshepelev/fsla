#include<assert.h>
#include"C1DField.h"

C1DField::C1DField(C1DProblem& pr)
	: imin(pr.get_order()),
	  imax(pr.get_order() + pr.nx),
	  x(vector<double>(pr.get_order() + pr.nx + pr.get_order() + 1)),
	  dx((pr.xmax - pr.xmin) / pr.nx),
	  t(0.), dt(0.),
	  U(vector<vector<double>>()),
	  newU(vector<vector<double>>()),
	  F(vector<vector<double>>()) {
	const std::size_t len = pr.get_order() + pr.nx + pr.get_order();

	U.resize(len);
	newU.resize(len);
	F.resize(len);

	const std::size_t variables_vector_size = 3;

	for (std::size_t i = 0; i < len/*U.size()*/; ++ i) {
		U[i].resize(variables_vector_size);
		newU[i].resize(variables_vector_size);
		F[i].resize(variables_vector_size);
	}
}

C1DField::~C1DField() {
	x.clear();
	U.clear();
	newU.clear();
	F.clear();
}
