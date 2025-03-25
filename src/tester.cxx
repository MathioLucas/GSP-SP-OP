// this tests on every .txt file in a subdirectory, given as a command-line argument

#include "gsp-sp-op.hxx"
#include <fstream>
#include <iostream>
#include <filesystem>

int main(int argc, char * argv[]) {
	if (argc < 2) return 1;

	auto dir = std::filesystem::path{argv[1]};

	if (!std::filesystem::is_directory(std::filesystem::status(dir))) {
		std::cout << "nondir\n";
		return 1;
	}
	std::cout << "dir " << dir << "\n\n";
	for (auto file : std::filesystem::directory_iterator{dir}) {
		if (file.path().extension() == ".txt") {
			std::cout << "\n======= FILE " << file << " =======\n";

			std::fstream fin{file.path()};

			if (!fin.is_open()) {
				std::cout << "opening error\n";
				continue;
			}

			graph g;
			fin >> g;
			V_LOG(g)

			gsp_sp_op_result r = GSP_SP_OP(g);
			bool win = r.authenticate(g);
			if (!win) {
				std::cout << "UH OH YOU'VE GONE AND MESSED IT UP YOU ABSOLUTE BUFFOON\n";
				break;
			}
		}
	}
}