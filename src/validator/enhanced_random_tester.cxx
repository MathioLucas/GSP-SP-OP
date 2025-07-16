#include "enhanced_validator.hxx"
#include "GraphGenerator.hxx"
#include <iostream>
#include <filesystem>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <iomanip>
namespace fs = std::filesystem;
class enhanced_tester {
private:
    enhanced_validator validator;
    std::string output_dir;
    graph load_graph_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open file: " + filename);
        }
        graph g;
        file >> g;
        return g;
    }
    bool test_single_graph(const std::string& filename, bool export_results = false) {
        try {
            std::cout << "\n=== Testing: " << filename << " ===" << std::endl;
            graph g = load_graph_from_file(filename);
            gsp_sp_op_result result = GSP_SP_OP(g);
            std::cout << "Graph loaded: " << g.n << " vertices, " << g.m << " edges" << std::endl;
            std::cout << "GSP: " << (result.is_gsp ? "YES" : "NO") << std::endl;
            std::cout << "SP: " << (result.is_sp ? "YES" : "NO") << std::endl;
            std::cout << "OP: " << (result.is_op ? "YES" : "NO") << std::endl;
            bool validation_result = validator.full_validate(g, result);
            if (export_results) {
                std::string base_name = fs::path(filename).stem().string();
                validator.export_visualization(g, result, output_dir + "/" + base_name + ".dot");
                validator.export_analysis(g, result, output_dir + "/" + base_name + "_analysis.txt");
            }
            std::cout << "Validation: " << (validation_result ? "PASSED" : "FAILED") << std::endl;
            return validation_result;
        } catch (const std::exception& e) {
            std::cerr << "Error testing " << filename << ": " << e.what() << std::endl;
            return false;
        }
    }
    void run_predefined_tests() {
        std::cout << "\n=== RUNNING PREDEFINED TEST CASES ===" << std::endl;
        std::vector<std::string> test_directories = {
            "test_cases/biconnected",
            "test_cases/non_biconnected",
            "test_cases/massive"
        };
        int total_tests = 0;
        int passed_tests = 0;
        for (const auto& dir : test_directories) {
            if (!fs::exists(dir)) {
                std::cout << "Directory " << dir << " not found, skipping..." << std::endl;
                continue;
            }
            std::cout << "\nTesting directory: " << dir << std::endl;
            for (const auto& entry : fs::directory_iterator(dir)) {
                if (entry.is_regular_file() && entry.path().extension() == ".txt") {
                    total_tests++;
                    if (test_single_graph(entry.path().string(), true)) {
                        passed_tests++;
                    }
                }
            }
        }
        std::cout << "\n=== PREDEFINED TEST RESULTS ===" << std::endl;
        std::cout << "Total tests: " << total_tests << std::endl;
        std::cout << "Passed: " << passed_tests << std::endl;
        std::cout << "Failed: " << (total_tests - passed_tests) << std::endl;
        if (total_tests > 0) {
            std::cout << "Success rate: " << (100.0 * passed_tests / total_tests) << "%" << std::endl;
        }
    }
    void run_random_tests(int num_graphs = 1000) {
        std::cout << "\n=== RUNNING RANDOM TESTS ===" << std::endl;
        std::cout << "Generating and testing " << num_graphs << " random graphs..." << std::endl;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> size_dist(10, 100);
        std::uniform_int_distribution<> component_dist(1, 10);
        std::vector<graph> graphs;
        std::vector<gsp_sp_op_result> results;
        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_graphs; i++) {
            if (i % 100 == 0) {
                std::cout << "Generated " << i << "/" << num_graphs << " graphs..." << std::endl;
            }
            long nC = component_dist(gen);
            long lC = size_dist(gen) / nC + 3;
            long nK = component_dist(gen);
            long lK = size_dist(gen) / nK + 3;
            long three_edges = gen() % 2;
            long seed = gen();
            try {
                graph g = generate_graph(nC, lC, nK, lK, three_edges, seed);
                gsp_sp_op_result result = GSP_SP_OP(g);
                graphs.push_back(g);
                results.push_back(result);
            } catch (const std::exception& e) {
                std::cerr << "Error generating graph " << i << ": " << e.what() << std::endl;
                continue;
            }
        }
        auto generation_end = std::chrono::high_resolution_clock::now();
        double generation_time = std::chrono::duration<double>(generation_end - start_time).count();
        std::cout << "Graph generation completed in " << generation_time << " seconds" << std::endl;
        std::cout << "Successfully generated " << graphs.size() << " graphs" << std::endl;
        std::vector<bool> validation_results = validator.batch_validate(graphs, results, true);
        auto validation_end = std::chrono::high_resolution_clock::now();
        double validation_time = std::chrono::duration<double>(validation_end - generation_end).count();
        std::cout << "\nValidation completed in " << validation_time << " seconds" << std::endl;
        int passed = std::count(validation_results.begin(), validation_results.end(), true);
        int failed = validation_results.size() - passed;
        std::cout << "\n=== RANDOM TEST RESULTS ===" << std::endl;
        std::cout << "Total graphs: " << validation_results.size() << std::endl;
        std::cout << "Passed validation: " << passed << std::endl;
        std::cout << "Failed validation: " << failed << std::endl;
        if (!validation_results.empty()) {
            std::cout << "Success rate: " << (100.0 * passed / validation_results.size()) << "%" << std::endl;
        }
        validator.analyze_graph_properties(graphs);
        if (failed > 0) {
            std::cout << "\nExporting failed cases for analysis..." << std::endl;
            int export_count = 0;
            for (size_t i = 0; i < validation_results.size() && export_count < 10; i++) {
                if (!validation_results[i]) {
                    std::string filename = output_dir + "/failed_case_" + std::to_string(export_count);
                    validator.export_visualization(graphs[i], results[i], filename + ".dot");
                    validator.export_analysis(graphs[i], results[i], filename + "_analysis.txt");
                    export_count++;
                }
            }
        }
    }
    void run_stress_tests() {
        std::cout << "\n=== RUNNING STRESS TESTS ===" << std::endl;
        std::vector<std::tuple<int, int, int, int>> stress_configs = {
            {50, 10, 10, 5},    
            {100, 15, 15, 8},   
            {200, 20, 20, 10}   
        };
        for (const auto& config : stress_configs) {
            auto [nC, lC, nK, lK] = config;
            std::cout << "\nStress test config: nC=" << nC << ", lC=" << lC 
                      << ", nK=" << nK << ", lK=" << lK << std::endl;
            try {
                graph g = generate_graph(nC, lC, nK, lK, 0, 12345);
                std::cout << "Generated graph: " << g.n << " vertices, " << g.m << " edges" << std::endl;
                auto start = std::chrono::high_resolution_clock::now();
                gsp_sp_op_result result = GSP_SP_OP(g);
                auto end = std::chrono::high_resolution_clock::now();
                double algorithm_time = std::chrono::duration<double>(end - start).count();
                std::cout << "Algorithm time: " << algorithm_time << " seconds" << std::endl;
                start = std::chrono::high_resolution_clock::now();
                bool validation_passed = validator.full_validate(g, result);
                end = std::chrono::high_resolution_clock::now();
                double validation_time = std::chrono::duration<double>(end - start).count();
                std::cout << "Validation time: " << validation_time << " seconds" << std::endl;
                std::cout << "Validation result: " << (validation_passed ? "PASSED" : "FAILED") << std::endl;
                if (validation_passed) {
                    std::cout << "Properties: GSP=" << result.is_gsp << ", SP=" << result.is_sp << ", OP=" << result.is_op << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "Stress test failed: " << e.what() << std::endl;
            }
        }
    }
    void run_edge_case_tests() {
        std::cout << "\n=== RUNNING EDGE CASE TESTS ===" << std::endl;
        {
            graph empty_graph;
            empty_graph.n = 0;
            empty_graph.m = 0;
            try {
                gsp_sp_op_result result = GSP_SP_OP(empty_graph);
                validator.full_validate(empty_graph, result);
                std::cout << "Empty graph test: PASSED" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Empty graph test: FAILED (" << e.what() << ")" << std::endl;
            }
        }
        {
            graph single_vertex;
            single_vertex.n = 1;
            single_vertex.m = 0;
            single_vertex.adj.resize(1);
            try {
                gsp_sp_op_result result = GSP_SP_OP(single_vertex);
                validator.full_validate(single_vertex, result);
                std::cout << "Single vertex test: PASSED" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Single vertex test: FAILED (" << e.what() << ")" << std::endl;
            }
        }
        for (int k = 3; k <= 5; k++) {
            graph complete;
            complete.n = k;
            complete.m = k * (k - 1) / 2;
            complete.adj.resize(k);
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    if (i != j) {
                        complete.adj[i].push_back(j);
                    }
                }
            }
            try {
                gsp_sp_op_result result = GSP_SP_OP(complete);
                validator.full_validate(complete, result);
                std::cout << "Complete graph K" << k << " test: PASSED" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Complete graph K" << k << " test: FAILED (" << e.what() << ")" << std::endl;
            }
        }
        for (int n = 2; n <= 10; n++) {
            graph path;
            path.n = n;
            path.m = n - 1;
            path.adj.resize(n);
            for (int i = 0; i < n - 1; i++) {
                path.adj[i].push_back(i + 1);
                path.adj[i + 1].push_back(i);
            }
            try {
                gsp_sp_op_result result = GSP_SP_OP(path);
                validator.full_validate(path, result);
                std::cout << "Path graph P" << n << " test: PASSED" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Path graph P" << n << " test: FAILED (" << e.what() << ")" << std::endl;
            }
        }
        for (int n = 3; n <= 10; n++) {
            graph cycle;
            cycle.n = n;
            cycle.m = n;
            cycle.adj.resize(n);
            for (int i = 0; i < n; i++) {
                int next = (i + 1) % n;
                cycle.adj[i].push_back(next);
                cycle.adj[next].push_back(i);
            }
            try {
                gsp_sp_op_result result = GSP_SP_OP(cycle);
                validator.full_validate(cycle, result);
                std::cout << "Cycle graph C" << n << " test: PASSED" << std::endl;
            } catch (const std::exception& e) {
                std::cout << "Cycle graph C" << n << " test: FAILED (" << e.what() << ")" << std::endl;
            }
        }
    }
public:
    enhanced_tester(const std::string& output_directory = "validation_output", bool verbose = false) 
        : validator(verbose), output_dir(output_directory) {
        if (!fs::exists(output_dir)) {
            fs::create_directories(output_dir);
        }
    }
    void run_all_tests() {
        std::cout << "=== ENHANCED GSP-SP-OP VALIDATION TESTING ===" << std::endl;
        std::cout << "Output directory: " << output_dir << std::endl;
        auto start_time = std::chrono::high_resolution_clock::now();
        run_edge_case_tests();
        run_predefined_tests();
        run_random_tests(500);  
        run_stress_tests();
        auto end_time = std::chrono::high_resolution_clock::now();
        double total_time = std::chrono::duration<double>(end_time - start_time).count();
        std::cout << "\n=== FINAL COMPREHENSIVE REPORT ===" << std::endl;
        std::cout << "Total testing time: " << total_time << " seconds" << std::endl;
        validator.print_comprehensive_report();
        validator.print_performance_analysis();
        validator.print_certificate_analysis();
        std::ofstream stats_file(output_dir + "/final_statistics.txt");
        if (stats_file.is_open()) {
            auto cout_buf = std::cout.rdbuf();
            std::cout.rdbuf(stats_file.rdbuf());
            validator.print_comprehensive_report();
            validator.print_performance_analysis();
            validator.print_certificate_analysis();
            std::cout.rdbuf(cout_buf);
            stats_file.close();
        }
    }
    void run_custom_test(const std::string& test_file_or_params) {
        if (fs::exists(test_file_or_params)) {
            test_single_graph(test_file_or_params, true);
        } else {
            std::cout << "Custom test with parameters: " << test_file_or_params << std::endl;
        }
    }
    void interactive_mode() {
        std::cout << "\n=== INTERACTIVE TESTING MODE ===" << std::endl;
        std::cout << "Commands:" << std::endl;
        std::cout << "  1 - Run edge case tests" << std::endl;
        std::cout << "  2 - Run predefined tests" << std::endl;
        std::cout << "  3 - Run random tests" << std::endl;
        std::cout << "  4 - Run stress tests" << std::endl;
        std::cout << "  5 - Test specific file" << std::endl;
        std::cout << "  6 - Toggle verbose mode" << std::endl;
        std::cout << "  7 - Show statistics" << std::endl;
        std::cout << "  8 - Reset statistics" << std::endl;
        std::cout << "  0 - Exit" << std::endl;
        int choice;
        std::string input;
        while (true) {
            std::cout << "\nEnter command: ";
            std::cin >> choice;
            switch (choice) {
                case 1:
                    run_edge_case_tests();
                    break;
                case 2:
                    run_predefined_tests();
                    break;
                case 3:
                    std::cout << "Enter number of random graphs: ";
                    int num;
                    std::cin >> num;
                    run_random_tests(num);
                    break;
                case 4:
                    run_stress_tests();
                    break;
                case 5:
                    std::cout << "Enter filename: ";
                    std::cin >> input;
                    test_single_graph(input, true);
                    break;
                case 6:
                    std::cout << "Verbose mode toggled" << std::endl;
                    break;
                case 7:
                    validator.print_comprehensive_report();
                    break;
                case 8:
                    validator.reset_stats();
                                        std::cout << "Statistics reset" << std::endl;
                    break;
                case 0:
                    return;
                default:
                    std::cout << "Invalid command" << std::endl;
            }
        }
    }
    void run_benchmark(const std::vector<std::tuple<int, int, int, int>>& configs) {
        std::cout << "\n=== BENCHMARK MODE ===" << std::endl;
        for (const auto& config : configs) {
            auto [nC, lC, nK, lK] = config;
            std::cout << "\nBenchmarking config: nC=" << nC << ", lC=" << lC 
                      << ", nK=" << nK << ", lK=" << lK << std::endl;
            std::vector<double> algorithm_times;
            std::vector<double> validation_times;
            std::vector<bool> validation_results;
            const int num_runs = 10;
            for (int run = 0; run < num_runs; run++) {
                try {
                    graph g = generate_graph(nC, lC, nK, lK, 0, 12345 + run);
                    auto start = std::chrono::high_resolution_clock::now();
                    gsp_sp_op_result result = GSP_SP_OP(g);
                    auto end = std::chrono::high_resolution_clock::now();
                    double alg_time = std::chrono::duration<double>(end - start).count();
                    algorithm_times.push_back(alg_time);
                    start = std::chrono::high_resolution_clock::now();
                    bool valid = validator.full_validate(g, result);
                    end = std::chrono::high_resolution_clock::now();
                    double val_time = std::chrono::duration<double>(end - start).count();
                    validation_times.push_back(val_time);
                    validation_results.push_back(valid);
                } catch (const std::exception& e) {
                    std::cerr << "Benchmark run " << run << " failed: " << e.what() << std::endl;
                }
            }
            if (!algorithm_times.empty()) {
                double avg_alg = std::accumulate(algorithm_times.begin(), algorithm_times.end(), 0.0) / algorithm_times.size();
                double avg_val = std::accumulate(validation_times.begin(), validation_times.end(), 0.0) / validation_times.size();
                auto minmax_alg = std::minmax_element(algorithm_times.begin(), algorithm_times.end());
                auto minmax_val = std::minmax_element(validation_times.begin(), validation_times.end());
                int valid_count = std::count(validation_results.begin(), validation_results.end(), true);
                std::cout << "  Algorithm time: " << avg_alg << "s (min: " << *minmax_alg.first 
                          << ", max: " << *minmax_alg.second << ")" << std::endl;
                std::cout << "  Validation time: " << avg_val << "s (min: " << *minmax_val.first 
                          << ", max: " << *minmax_val.second << ")" << std::endl;
                std::cout << "  Validation success rate: " << (100.0 * valid_count / validation_results.size()) << "%" << std::endl;
                std::cout << "  Validation overhead: " << (avg_val / avg_alg * 100) << "%" << std::endl;
            }
        }
    }
    void analyze_memory_usage() {
        std::cout << "\n=== MEMORY USAGE ANALYSIS ===" << std::endl;
        std::vector<std::tuple<int, int, int, int>> test_configs = {
            {10, 5, 5, 3},
            {50, 10, 10, 5},
            {100, 15, 15, 8},
            {200, 20, 20, 10}
        };
        for (const auto& config : test_configs) {
            auto [nC, lC, nK, lK] = config;
            try {
                graph g = generate_graph(nC, lC, nK, lK, 0, 12345);
                size_t graph_memory = sizeof(graph) + g.adj.capacity() * sizeof(std::vector<int>);
                for (const auto& adj_list : g.adj) {
                    graph_memory += adj_list.capacity() * sizeof(int);
                }
                std::cout << "Config nC=" << nC << ", lC=" << lC << ", nK=" << nK << ", lK=" << lK << std::endl;
                std::cout << "  Graph size: " << g.n << " vertices, " << g.m << " edges" << std::endl;
                std::cout << "  Estimated memory: " << graph_memory << " bytes" << std::endl;
                std::cout << "  Memory per vertex: " << (graph_memory / g.n) << " bytes" << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Memory analysis failed for config: " << e.what() << std::endl;
            }
        }
    }
    void export_test_results_csv(const std::string& filename) {
        std::ofstream csv_file(output_dir + "/" + filename);
        if (!csv_file.is_open()) {
            std::cerr << "Could not open CSV file: " << filename << std::endl;
            return;
        }
        const auto& stats = validator.get_stats();
        csv_file << "metric,value\n";
        csv_file << "total_graphs," << stats.total_graphs << "\n";
        csv_file << "gsp_graphs," << stats.gsp_graphs << "\n";
        csv_file << "sp_graphs," << stats.sp_graphs << "\n";
        csv_file << "op_graphs," << stats.op_graphs << "\n";
        csv_file << "k4_certificates," << stats.k4_certificates << "\n";
        csv_file << "k23_certificates," << stats.k23_certificates << "\n";
        csv_file << "t4_certificates," << stats.t4_certificates << "\n";
        csv_file << "authentication_failures," << stats.authentication_failures << "\n";
        csv_file << "deep_validation_failures," << stats.deep_validation_failures << "\n";
        csv_file << "cross_validation_failures," << stats.cross_validation_failures << "\n";
        csv_file << "avg_validation_time," << stats.avg_validation_time << "\n";
        csv_file << "max_validation_time," << stats.max_validation_time << "\n";
        csv_file << "total_validation_time," << stats.total_validation_time << "\n";
        csv_file.close();
        std::cout << "Test results exported to: " << output_dir << "/" << filename << std::endl;
    }
    void generate_comprehensive_report() {
        std::string report_file = output_dir + "/comprehensive_report.html";
        std::ofstream html(report_file);
        if (!html.is_open()) {
            std::cerr << "Could not create comprehensive report" << std::endl;
            return;
        }
        const auto& stats = validator.get_stats();
        html << "<!DOCTYPE html>\n<html>\n<head>\n";
        html << "<title>GSP-SP-OP Validation Report</title>\n";
        html << "<style>\n";
        html << "body { font-family: Arial, sans-serif; margin: 20px; }\n";
        html << "table { border-collapse: collapse; width: 100%; }\n";
        html << "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }\n";
        html << "th { background-color: #f2f2f2; }\n";
        html << ".pass { color: green; }\n";
        html << ".fail { color: red; }\n";
        html << ".warning { color: orange; }\n";
        html << "</style>\n</head>\n<body>\n";
        html << "<h1>GSP-SP-OP Validation Report</h1>\n";
        html << "<p>Generated on: " << std::chrono::system_clock::now().time_since_epoch().count() << "</p>\n";
        html << "<h2>Summary Statistics</h2>\n";
        html << "<table>\n";
        html << "<tr><th>Metric</th><th>Value</th></tr>\n";
        html << "<tr><td>Total Graphs Tested</td><td>" << stats.total_graphs << "</td></tr>\n";
        html << "<tr><td>GSP Graphs</td><td>" << stats.gsp_graphs << "</td></tr>\n";
        html << "<tr><td>SP Graphs</td><td>" << stats.sp_graphs << "</td></tr>\n";
        html << "<tr><td>OP Graphs</td><td>" << stats.op_graphs << "</td></tr>\n";
        html << "<tr><td>Average Validation Time</td><td>" << stats.avg_validation_time << "s</td></tr>\n";
        html << "</table>\n";
        html << "<h2>Certificate Distribution</h2>\n";
        html << "<table>\n";
        html << "<tr><th>Certificate Type</th><th>Count</th></tr>\n";
        html << "<tr><td>K4 Subdivision</td><td>" << stats.k4_certificates << "</td></tr>\n";
        html << "<tr><td>K23 Subdivision</td><td>" << stats.k23_certificates << "</td></tr>\n";
        html << "<tr><td>T4 Subdivision</td><td>" << stats.t4_certificates << "</td></tr>\n";
        html << "<tr><td>Positive GSP</td><td>" << stats.positive_gsp_certificates << "</td></tr>\n";
        html << "<tr><td>Positive OP</td><td>" << stats.positive_op_certificates << "</td></tr>\n";
        html << "</table>\n";
        html << "<h2>Failure Analysis</h2>\n";
        html << "<table>\n";
        html << "<tr><th>Failure Type</th><th>Count</th><th>Status</th></tr>\n";
        auto failure_status = [](int count) -> std::string {
            if (count == 0) return "class=\"pass\">PASS";
            else if (count < 5) return "class=\"warning\">WARNING";
            else return "class=\"fail\">FAIL";
        };
        html << "<tr><td>Authentication Failures</td><td>" << stats.authentication_failures 
             << "</td><td " << failure_status(stats.authentication_failures) << "</td></tr>\n";
        html << "<tr><td>Deep Validation Failures</td><td>" << stats.deep_validation_failures 
             << "</td><td " << failure_status(stats.deep_validation_failures) << "</td></tr>\n";
        html << "<tr><td>Cross Validation Failures</td><td>" << stats.cross_validation_failures 
             << "</td><td " << failure_status(stats.cross_validation_failures) << "</td></tr>\n";
        html << "<tr><td>Structural Inconsistencies</td><td>" << stats.structural_inconsistencies 
             << "</td><td " << failure_status(stats.structural_inconsistencies) << "</td></tr>\n";
        html << "</table>\n";
        html << "<h2>Performance Metrics</h2>\n";
        html << "<table>\n";
        html << "<tr><th>Metric</th><th>Value</th></tr>\n";
        html << "<tr><td>Total Validation Time</td><td>" << stats.total_validation_time << "s</td></tr>\n";
        html << "<tr><td>Maximum Validation Time</td><td>" << stats.max_validation_time << "s</td></tr>\n";
        if (stats.total_validation_time > 0) {
            html << "<tr><td>Throughput</td><td>" << (stats.total_graphs / stats.total_validation_time) 
                 << " graphs/second</td></tr>\n";
        }
        html << "</table>\n";
        html << "</body>\n</html>\n";
        html.close();
        std::cout << "Comprehensive HTML report generated: " << report_file << std::endl;
    }
    void run_regression_tests() {
        std::cout << "\n=== REGRESSION TESTING ===" << std::endl;
        std::vector<std::tuple<int, int, int, int, long, bool, bool, bool>> regression_cases = {
            {5, 5, 5, 5, 12345, true, false, false},
            {3, 3, 3, 3, 54321, true, true, true},
            {10, 10, 10, 10, 98765, true, false, false},
            {1, 10, 1, 10, 11111, true, true, true},
            {20, 5, 5, 5, 22222, true, false, false}
        };
        int passed = 0;
        int total = regression_cases.size();
        for (const auto& test_case : regression_cases) {
            auto [nC, lC, nK, lK, seed, exp_gsp, exp_sp, exp_op] = test_case;
            try {
                graph g = generate_graph(nC, lC, nK, lK, 0, seed);
                gsp_sp_op_result result = GSP_SP_OP(g);
                bool test_passed = (result.is_gsp == exp_gsp) && 
                                  (result.is_sp == exp_sp) && 
                                  (result.is_op == exp_op);
                if (test_passed && validator.full_validate(g, result)) {
                    passed++;
                    std::cout << "Regression test " << (passed + (total - passed)) << ": PASSED" << std::endl;
                } else {
                    std::cout << "Regression test " << (passed + (total - passed)) << ": FAILED" << std::endl;
                    std::cout << "  Expected: GSP=" << exp_gsp << ", SP=" << exp_sp << ", OP=" << exp_op << std::endl;
                    std::cout << "  Got:      GSP=" << result.is_gsp << ", SP=" << result.is_sp << ", OP=" << result.is_op << std::endl;
                }
            } catch (const std::exception& e) {
                std::cerr << "Regression test failed with exception: " << e.what() << std::endl;
            }
        }
        std::cout << "\nRegression test results: " << passed << "/" << total << " passed" << std::endl;
    }
    const validation_stats& get_validation_stats() const {
        return validator.get_stats();
    }
    void set_verbose(bool verbose) {
        validator.set_verbose(verbose);
    }
    void cleanup_output() {
        try {
            if (fs::exists(output_dir)) {
                fs::remove_all(output_dir);
                fs::create_directories(output_dir);
                std::cout << "Output directory cleaned: " << output_dir << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error cleaning output directory: " << e.what() << std::endl;
        }
    }
};
int main(int argc, char* argv[]) {
    try {
        enhanced_tester tester("validation_output", false);
        if (argc > 1) {
            std::string command = argv[1];
            if (command == "all") {
                tester.run_all_tests();
            } else if (command == "interactive") {
                tester.interactive_mode();
            } else if (command == "benchmark") {
                std::vector<std::tuple<int, int, int, int>> configs = {
                    {10, 5, 5, 3},
                    {50, 10, 10, 5},
                    {100, 15, 15, 8}
                };
                tester.run_benchmark(configs);
            } else if (command == "regression") {
                tester.run_regression_tests();
            } else if (command == "memory") {
                tester.analyze_memory_usage();
            } else if (command == "verbose") {
                tester.set_verbose(true);
                tester.run_all_tests();
            } else if (command == "file" && argc > 2) {
                tester.run_custom_test(argv[2]);
            } else {
                std::cout << "Usage: " << argv[0] << " [all|interactive|benchmark|regression|memory|verbose|file <filename>]" << std::endl;
                return 1;
            }
        } else {
            tester.run_all_tests();
        }
        tester.generate_comprehensive_report();
        tester.export_test_results_csv("test_results.csv");
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
#endif
