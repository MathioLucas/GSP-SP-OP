#include "gsp-sp-op.hxx"
#include "enhanced_validator.hxx"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <map>
int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <test_directory>" << std::endl;
        std::cerr << "Example: " << argv[0] << " \"test cases/biconnected\"" << std::endl;
        return 1;
    }
    std::string test_dir = argv[1];
    std::cout << "=== ENHANCED TESTER ===" << std::endl;
    std::cout << "Testing directory: " << test_dir << std::endl;
    std::cout << "Enhanced with deep validation and cross-checking" << std::endl;
    enhanced_validator validator;
    int total_files = 0;
    int successful_tests = 0;
    int failed_tests = 0;
    int auth_failures = 0;
    int deep_validation_failures = 0;
    int cross_validation_failures = 0;
    std::map<std::string, int> error_categories;
    auto start_time = std::chrono::high_resolution_clock::now();
    try {
        std::filesystem::create_directories("test_results");
        for (const auto& entry : std::filesystem::directory_iterator(test_dir)) {
            if (entry.path().extension() == ".txt") {
                total_files++;
                std::string filename = entry.path().string();
                std::string basename = entry.path().stem().string();
                std::cout << "\n====== TESTING FILE: " << filename << " ======" << std::endl;
                std::ifstream file(filename);
                if (!file.is_open()) {
                    std::cerr << "Could not open file: " << filename << std::endl;
                    failed_tests++;
                    error_categories["file_read_error"]++;
                    continue;
                }
                graph g;
                try {
                    file >> g;
                    file.close();
                    std::cout << "Graph loaded: " << g.n << " vertices, " << g.m << " edges" << std::endl;
                    if (g.n <= 0 || g.m < 0) {
                        std::cout << "✗ INVALID GRAPH PARAMETERS!" << std::endl;
                        failed_tests++;
                        error_categories["invalid_graph"]++;
                        continue;
                    }
                    gsp_sp_op_result result = GSP_SP_OP(g);
                    std::cout << "=== ALGORITHM RESULTS ===" << std::endl;
                    std::cout << "GSP (Generalized Series-Parallel): " << (result.is_gsp ? "YES" : "NO") << std::endl;
                    std::cout << "SP (Series-Parallel): " << (result.is_sp ? "YES" : "NO") << std::endl;
                    std::cout << "OP (Outerplanar): " << (result.is_op ? "YES" : "NO") << std::endl;
                    std::cout << "\n=== CERTIFICATE ANALYSIS ===" << std::endl;
                    if (result.gsp_reason) {
                        if (dynamic_cast<const positive_cert_gsp*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Positive (GSP decomposition tree)" << std::endl;
                        } else if (dynamic_cast<const negative_cert_K4*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Negative (K4 subdivision)" << std::endl;
                        } else if (dynamic_cast<const negative_cert_K23*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Negative (K23 subdivision)" << std::endl;
                        } else if (dynamic_cast<const negative_cert_T4*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Negative (T4 subdivision)" << std::endl;
                        } else if (dynamic_cast<const negative_cert_tri_comp_cut*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Negative (tri-comp-cut)" << std::endl;
                        } else if (dynamic_cast<const negative_cert_tri_cut_comp*>(result.gsp_reason)) {
                            std::cout << "GSP Certificate: Negative (tri-cut-comp)" << std::endl;
                        }
                    }
                    if (result.op_reason) {
                        if (dynamic_cast<const positive_cert_op*>(result.op_reason)) {
                            std::cout << "OP Certificate: Positive (outerplanar embedding)" << std::endl;
                        } else {
                            std::cout << "OP Certificate: Negative (non-outerplanar)" << std::endl;
                        }
                    }
                    std::cout << "\n=== STANDARD AUTHENTICATION ===" << std::endl;
                    bool standard_auth = result.authenticate(g);
                    if (!standard_auth) {
                        std::cout << "✗ STANDARD AUTHENTICATION FAILED!" << std::endl;
                        failed_tests++;
                        auth_failures++;
                        error_categories["auth_failure"]++;
                        std::string failure_viz = "test_results/" + basename + "_auth_failure.dot";
                        validator.export_visualization(g, result, failure_viz);
                        continue;
                    }
                    std::cout << "✓ Standard authentication successful" << std::endl;
                    std::cout << "\n=== DEEP VALIDATION ===" << std::endl;
                    bool deep_valid = validator.deep_validate(g, result);
                    if (!deep_valid) {
                        std::cout << "✗ DEEP VALIDATION FAILED!" << std::endl;
                        deep_validation_failures++;
                        error_categories["deep_validation_failure"]++;
                        std::string deep_failure_viz = "test_results/" + basename + "_deep_failure.dot";
                        validator.export_visualization(g, result, deep_failure_viz);
                    } else {
                        std::cout << "✓ Deep validation successful" << std::endl;
                    }
                    std::cout << "\n=== CROSS VALIDATION ===" << std::endl;
                    bool cross_valid = validator.cross_validate(g, result);
                    if (!cross_valid) {
                        std::cout << "✗ CROSS VALIDATION FAILED!" << std::endl;
                        cross_validation_failures++;
                        error_categories["cross_validation_failure"]++;
                        std::string cross_failure_viz = "test_results/" + basename + "_cross_failure.dot";
                        validator.export_visualization(g, result, cross_failure_viz);
                    } else {
                        std::cout << "✓ Cross validation successful" << std::endl;
                    }
                    std::string viz_filename = "test_results/" + basename + "_visualization.dot";
                    validator.export_visualization(g, result, viz_filename);
                    std::cout << "Visualization exported to: " << viz_filename << std::endl;
                    std::string report_filename = "test_results/" + basename + "_report.txt";
                    std::ofstream report(report_filename);
                    report << "=== DETAILED REPORT FOR " << filename << " ===" << std::endl;
                    report << "Graph: " << g.n << " vertices, " << g.m << " edges" << std::endl;
                    report << "GSP: " << (result.is_gsp ? "YES" : "NO") << std::endl;
                    report << "SP: " << (result.is_sp ? "YES" : "NO") << std::endl;
                    report << "OP: " << (result.is_op ? "YES" : "NO") << std::endl;
                    report << "Standard Authentication: " << (standard_auth ? "PASS" : "FAIL") << std::endl;
                    report << "Deep Validation: " << (deep_valid ? "PASS" : "FAIL") << std::endl;
                    report << "Cross Validation: " << (cross_valid ? "PASS" : "FAIL") << std::endl;
                    report << "\n=== ADJACENCY LIST ===" << std::endl;
                    for (int i = 0; i < g.n; i++) {
                        report << "Vertex " << i << ": ";
                        for (int j : g.adj[i]) {
                            report << j << " ";
                        }
                        report << std::endl;
                    }
                    report << "\n=== THEORETICAL ANALYSIS ===" << std::endl;
                    report << "Density: " << (double)g.m / (g.n * (g.n - 1) / 2) << std::endl;
                    report << "Average degree: " << (double)(2 * g.m) / g.n << std::endl;
                    if (g.n >= 3) {
                        int max_planar_edges = 3 * g.n - 6;
                        int max_outerplanar_edges = 2 * g.n - 3;
                        report << "Max planar edges: " << max_planar_edges << std::endl;
                        report << "Max outerplanar edges: " << max_outerplanar_edges << std::endl;
                        report << "Planar feasible: " << (g.m <= max_planar_edges ? "YES" : "NO") << std::endl;
                        report << "Outerplanar feasible: " << (g.m <= max_outerplanar_edges ? "YES" : "NO") << std::endl;
                    }
                    report.close();
                    if (standard_auth && deep_valid && cross_valid) {
                        std::cout << "✓ ALL VALIDATIONS PASSED" << std::endl;
                        successful_tests++;
                    } else {
                        std::cout << "⚠ SOME VALIDATIONS FAILED" << std::endl;
                    }
                } catch (const std::exception& e) {
                    std::cout << "✗ EXCEPTION DURING PROCESSING: " << e.what() << std::endl;
                    failed_tests++;
                    error_categories["exception"]++;
                    continue;
                }
            }
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
        std::cout << "\n========================================" << std::endl;
        std::cout << "=== COMPREHENSIVE TEST RESULTS ===" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Total files processed: " << total_files << std::endl;
        std::cout << "Successful tests: " << successful_tests << std::endl;
        std::cout << "Failed tests: " << failed_tests << std::endl;
        std::cout << "Authentication failures: " << auth_failures << std::endl;
        std::cout << "Deep validation failures: " << deep_validation_failures << std::endl;
        std::cout << "Cross validation failures: " << cross_validation_failures << std::endl;
        std::cout << "Processing time: " << duration.count() << " ms" << std::endl;
        std::cout << "\n=== ERROR CATEGORY BREAKDOWN ===" << std::endl;
        for (const auto& category : error_categories) {
            std::cout << category.first << ": " << category.second << std::endl;
        }
        double success_rate = (double)successful_tests / total_files * 100.0;
        std::cout << "\n=== SUCCESS METRICS ===" << std::endl;
        std::cout << "Overall success rate: " << success_rate << "%" << std::endl;
        if (success_rate >= 95.0) {
            std::cout << "✓ EXCELLENT - System is highly reliable" << std::endl;
        } else if (success_rate >= 90.0) {
            std::cout << "✓ GOOD - System is reliable with minor issues" << std::endl;
        } else if (success_rate >= 80.0) {
            std::cout << "⚠ FAIR - System has some reliability issues" << std::endl;
        } else {
            std::cout << "✗ POOR - System has significant reliability issues" << std::endl;
        }
        std::cout << "\n";
        validator.print_comprehensive_report();
        std::ofstream summary("test_results/summary_report.txt");
        summary << "=== COMPREHENSIVE TEST SUMMARY ===" << std::endl;
        summary << "Test Directory: " << test_dir << std::endl;
        summary << "Total Files: " << total_files << std::endl;
        summary << "Successful Tests: " << successful_tests << std::endl;
        summary << "Failed Tests: " << failed_tests << std::endl;
        summary << "Success Rate: " << success_rate << "%" << std::endl;
        summary << "Processing Time: " << duration.count() << " ms" << std::endl;
        summary << "\nError Categories:" << std::endl;
        for (const auto& category : error_categories) {
            summary << "  " << category.first << ": " << category.second << std::endl;
        }
        summary.close();
        std::cout << "\n=== OUTPUT FILES ===" << std::endl;
        std::cout << "All results exported to test_results/ directory" << std::endl;
        std::cout << "- Individual graph reports: *_report.txt" << std::endl;
        std::cout << "- Visualizations: *_visualization.dot" << std::endl;
        std::cout << "- Failure cases: *_failure.dot" << std::endl;
        std::cout << "- Summary report: summary_report.txt" << std::endl;
        std::cout << "\nUse 'dot -Tpng filename.dot -o filename.png' to view visualizations" << std::endl;
        if (failed_tests > 0) {
            std::cout << "\n✗ CRITICAL FAILURES DETECTED - Review failed cases" << std::endl;
            return 1;
        } else if (deep_validation_failures > 0 || cross_validation_failures > 0) {
            std::cout << "\n⚠ MINOR VALIDATION ISSUES - Review validation failures" << std::endl;
            return 2;
        } else {
            std::cout << "\n✓ ALL TESTS PASSED SUCCESSFULLY" << std::endl;
            return 0;
        }
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        return 1;
    }
}
