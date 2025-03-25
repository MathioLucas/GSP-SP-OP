// for debugging purposes, this recreates a failed randomly generated test case

#include "gsp-sp-op.hxx"
#include "GraphGenerator.hxx"
#include <cstdlib>
#include <iostream>

int main(int argc, char * argv[]) {
	if (argc < 7) return 0;

	graph g = generate_graph(atol(argv[1]), atol(argv[2]), atol(argv[3]), atol(argv[4]), atol(argv[5]), atol(argv[6]));

	gsp_sp_op_result r = GSP_SP_OP(g);
	bool win = r.authenticate(g);
	if (!win) {
		std::cout << "UH OH YOU'VE GONE AND MESSED IT UP YOU ABSOLUTE BUFFOON\n";
	} else {
		std::cout << "YAAAY\n";
	}
}