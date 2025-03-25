// this tests on a bunch of randomly generated graphs

#include "gsp-sp-op.hxx"
#include "GraphGenerator.hxx"
#include <iostream>
#include <random>
#include <limits>

int main(int argc, char * argv[]) {
	std::random_device rd;
	std::default_random_engine re{rd()};

	std::uniform_int_distribution<long> seed{std::numeric_limits<long>::min(), std::numeric_limits<long>::max()};
	std::uniform_int_distribution<long> nC{2, 100};
	std::uniform_int_distribution<long> lC{3, 10};
	std::uniform_int_distribution<long> nK{0, 1};
	std::uniform_int_distribution<long> lK{4, 10};
	std::uniform_int_distribution<long> three_edges{0, 1};

	for (int i = 0; i < 100000; i++) { // testing testing 123
		long nC_ = nC(rd);
		long lC_ = lC(rd);
		long nK_ = nK(rd);
		long lK_ = lK(rd);
		long three_edges_ = three_edges(rd);
		long seed_ = seed(rd);
		graph g = generate_graph(nC_, lC_, nK_, lK_, three_edges_, seed_);

		std::cout << "\n====== TESTING GRAPH " << i << " ======\n";
		gsp_sp_op_result r = GSP_SP_OP(g);
		bool win = r.authenticate(g);
		if (!win) {
			std::cout << "UH OH YOU'VE GONE AND MESSED IT UP YOU ABSOLUTE BUFFOON\n";
			std::cout << "graph: " << i << ", params: " << nC_ << " " << lC_ << " " << nK_ << " " << lK_ << " " << three_edges_ << " " << seed_ << "\n";
			break;
		}
	}

	/*gsp_sp_op_result tmp{};

	std::vector<int> cut_verts;
	cut_verts.reserve(g.n);
	for (int i = 0; i < g.n; i++) {
		cut_verts.emplace_back(-1);
	}
	auto bicomps = get_bicomps(g, cut_verts, tmp);*/

	/*std::fstream fout{"test cases/massive/K1_1000000.txt"};

	for (int i = 1; i <= 1000000; i++) {
		fout << "0 " << i << "\n";
	}*/

	/*std::random_device rd;
	std::default_random_engine re{rd()};
	std::uniform_int_distribution<int> dist{1, 1000000000};
	int n = 10000000;

	std::vector<int> v;
	v.reserve(n);
	std::generate_n(std::back_inserter(v), n, [&](){ return dist(re); });

	radix_sort(v);

	for (int i = 0; i < n; i++) {
		std::cout << v[i] << " ";
	}


	for (int bit = 28; bit >= 0; bit -= 4) {
		std::cout << "========= BIT " << bit << ": ===========\n";
		int mask = 15 << bit;
		std::cout << "mask " << mask << "\n";
		for (int i = 0; i < 256; i++) {
			int bucket = (i & mask) >> bit;
			std::cout << "num " << i << ": bucket " << bucket << "\n";
		}
	}*/
}