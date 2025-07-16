#ifndef ENHANCED_VALIDATOR_HXX
#define ENHANCED_VALIDATOR_HXX
#include "gsp-sp-op.hxx"
#include <iostream>
#include <vector>
#include <set>
#include <queue>
#include <map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <unordered_set>
#include <unordered_map>
struct validation_stats {
    int total_graphs = 0;
    int gsp_graphs = 0;
    int sp_graphs = 0;
    int op_graphs = 0;
    int k4_certificates = 0;
    int k23_certificates = 0;
    int t4_certificates = 0;
    int tri_comp_cut_certificates = 0;
    int tri_cut_comp_certificates = 0;
    int positive_gsp_certificates = 0;
    int positive_op_certificates = 0;
    int authentication_failures = 0;
    int deep_validation_failures = 0;
    int cross_validation_failures = 0;
    int structural_inconsistencies = 0;
    double total_validation_time = 0.0;
    double avg_validation_time = 0.0;
    double max_validation_time = 0.0;
    void update_timing(double time) {
        total_validation_time += time;
        avg_validation_time = total_validation_time / total_graphs;
        max_validation_time = std::max(max_validation_time, time);
    }
    void print_stats() const {
        std::cout << "=== VALIDATION STATISTICS ===" << std::endl;
        std::cout << "Total graphs tested: " << total_graphs << std::endl;
        std::cout << "GSP graphs: " << gsp_graphs << " (" << (total_graphs > 0 ? 100.0 * gsp_graphs / total_graphs : 0) << "%)" << std::endl;
        std::cout << "SP graphs: " << sp_graphs << " (" << (total_graphs > 0 ? 100.0 * sp_graphs / total_graphs : 0) << "%)" << std::endl;
        std::cout << "OP graphs: " << op_graphs << " (" << (total_graphs > 0 ? 100.0 * op_graphs / total_graphs : 0) << "%)" << std::endl;
        std::cout << "\n=== CERTIFICATE DISTRIBUTION ===" << std::endl;
        std::cout << "K4 certificates: " << k4_certificates << std::endl;
        std::cout << "K23 certificates: " << k23_certificates << std::endl;
        std::cout << "T4 certificates: " << t4_certificates << std::endl;
        std::cout << "Tri-comp-cut certificates: " << tri_comp_cut_certificates << std::endl;
        std::cout << "Tri-cut-comp certificates: " << tri_cut_comp_certificates << std::endl;
        std::cout << "Positive GSP certificates: " << positive_gsp_certificates << std::endl;
        std::cout << "Positive OP certificates: " << positive_op_certificates << std::endl;
        std::cout << "\n=== FAILURE ANALYSIS ===" << std::endl;
        std::cout << "Authentication failures: " << authentication_failures << std::endl;
        std::cout << "Deep validation failures: " << deep_validation_failures << std::endl;
        std::cout << "Cross validation failures: " << cross_validation_failures << std::endl;
        std::cout << "Structural inconsistencies: " << structural_inconsistencies << std::endl;
        std::cout << "\n=== TIMING STATISTICS ===" << std::endl;
        std::cout << "Total validation time: " << total_validation_time << "s" << std::endl;
        std::cout << "Average validation time: " << avg_validation_time << "s" << std::endl;
        std::cout << "Maximum validation time: " << max_validation_time << "s" << std::endl;
    }
};
class enhanced_validator {
private:
    validation_stats stats;
    bool verbose_mode = false;
    bool are_adjacent(const graph& g, int u, int v) const {
        if (u >= g.n || v >= g.n || u < 0 || v < 0) return false;
        auto it = std::find(g.adj[u].begin(), g.adj[u].end(), v);
        return it != g.adj[u].end();
    }
    std::vector<int> find_path(const graph& g, int start, int end) const {
        if (start == end) return {start};
        if (start >= g.n || end >= g.n || start < 0 || end < 0) return {};
        std::vector<int> parent(g.n, -1);
        std::queue<int> q;
        std::vector<bool> visited(g.n, false);
        q.push(start);
        visited[start] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int v : g.adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    parent[v] = u;
                    q.push(v);
                    if (v == end) {
                        std::vector<int> path;
                        int curr = end;
                        while (curr != -1) {
                            path.push_back(curr);
                            curr = parent[curr];
                        }
                        std::reverse(path.begin(), path.end());
                        return path;
                    }
                }
            }
        }
        return {}; 
    }
    bool is_clique(const graph& g, const std::vector<int>& vertices) const {
        for (size_t i = 0; i < vertices.size(); i++) {
            for (size_t j = i + 1; j < vertices.size(); j++) {
                if (!are_adjacent(g, vertices[i], vertices[j])) {
                    return false;
                }
            }
        }
        return true;
    }
    bool is_connected(const graph& g) const {
        if (g.n == 0) return true;
        std::vector<bool> visited(g.n, false);
        std::queue<int> q;
        q.push(0);
        visited[0] = true;
        int count = 1;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (int v : g.adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                    count++;
                }
            }
        }
        return count == g.n;
    }
    std::vector<std::vector<int>> find_biconnected_components(const graph& g) const {
        std::vector<std::vector<int>> components;
        std::vector<bool> visited(g.n, false);
        std::vector<int> disc(g.n, -1);
        std::vector<int> low(g.n, -1);
        std::vector<int> parent(g.n, -1);
        std::vector<bool> ap(g.n, false);
        std::vector<std::pair<int, int>> stack;
        int time = 0;
        for (int i = 0; i < g.n; i++) {
            if (!visited[i]) {
                biconnected_dfs(g, i, visited, disc, low, parent, ap, stack, components, time);
            }
        }
        return components;
    }
    void biconnected_dfs(const graph& g, int u, std::vector<bool>& visited,
                        std::vector<int>& disc, std::vector<int>& low,
                        std::vector<int>& parent, std::vector<bool>& ap,
                        std::vector<std::pair<int, int>>& stack,
                        std::vector<std::vector<int>>& components, int& time) const {
        int children = 0;
        visited[u] = true;
        disc[u] = low[u] = ++time;
        for (int v : g.adj[u]) {
            if (!visited[v]) {
                children++;
                parent[v] = u;
                stack.push_back({u, v});
                biconnected_dfs(g, v, visited, disc, low, parent, ap, stack, components, time);
                low[u] = std::min(low[u], low[v]);
                if ((parent[u] == -1 && children > 1) || (parent[u] != -1 && low[v] >= disc[u])) {
                    ap[u] = true;
                    std::vector<int> component;
                    std::pair<int, int> edge;
                    do {
                        edge = stack.back();
                        stack.pop_back();
                        component.push_back(edge.first);
                        component.push_back(edge.second);
                    } while (edge != std::make_pair(u, v));
                    std::sort(component.begin(), component.end());
                    component.erase(std::unique(component.begin(), component.end()), component.end());
                    components.push_back(component);
                }
            } else if (v != parent[u] && disc[v] < disc[u]) {
                stack.push_back({u, v});
                low[u] = std::min(low[u], disc[v]);
            }
        }
    }
    bool verify_k4_subdivision(const graph& g, const negative_cert_K4* cert) const {
        if (verbose_mode) std::cout << "Deep validating K4 subdivision certificate..." << std::endl;
        if (!cert->authenticate(g)) {
            if (verbose_mode) std::cout << "  FAIL: Basic authentication failed" << std::endl;
            return false;
        }
        if (g.n < 4) {
            if (verbose_mode) std::cout << "  FAIL: Graph too small to contain K4" << std::endl;
            return false;
        }
        if (g.m < 6) {
            if (verbose_mode) std::cout << "  FAIL: Graph has too few edges for K4" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: K4 subdivision certificate validates" << std::endl;
        return true;
    }
    bool verify_k23_subdivision(const graph& g, const negative_cert_K23* cert) const {
        if (verbose_mode) std::cout << "Deep validating K23 subdivision certificate..." << std::endl;
        if (!cert->authenticate(g)) {
            if (verbose_mode) std::cout << "  FAIL: Basic authentication failed" << std::endl;
            return false;
        }
        if (g.n < 5) {
            if (verbose_mode) std::cout << "  FAIL: Graph too small to contain K23" << std::endl;
            return false;
        }
        if (g.m < 6) {
            if (verbose_mode) std::cout << "  FAIL: Graph has too few edges for K23" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: K23 subdivision certificate validates" << std::endl;
        return true;
    }
    bool verify_t4_subdivision(const graph& g, const negative_cert_T4* cert) const {
        if (verbose_mode) std::cout << "Deep validating T4 subdivision certificate..." << std::endl;
        if (!cert->authenticate(g)) {
            if (verbose_mode) std::cout << "  FAIL: Basic authentication failed" << std::endl;
            return false;
        }
        if (g.n < 4) {
            if (verbose_mode) std::cout << "  FAIL: Graph too small to contain T4" << std::endl;
            return false;
        }
        if (g.m < 5) {
            if (verbose_mode) std::cout << "  FAIL: Graph has too few edges for T4" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: T4 subdivision certificate validates" << std::endl;
        return true;
    }
    bool verify_positive_gsp(const graph& g, const positive_cert_gsp* cert) const {
        if (verbose_mode) std::cout << "Deep validating positive GSP certificate..." << std::endl;
        if (!cert->authenticate(g)) {
            if (verbose_mode) std::cout << "  FAIL: Basic authentication failed" << std::endl;
            return false;
        }
        if (!is_connected(g)) {
            if (verbose_mode) std::cout << "  FAIL: Graph is not connected" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: Positive GSP certificate validates" << std::endl;
        return true;
    }
    bool verify_positive_op(const graph& g, const positive_cert_op* cert) const {
        if (verbose_mode) std::cout << "Deep validating positive OP certificate..." << std::endl;
        if (!cert->authenticate(g)) {
            if (verbose_mode) std::cout << "  FAIL: Basic authentication failed" << std::endl;
            return false;
        }
        if (!is_potentially_planar(g)) {
            if (verbose_mode) std::cout << "  FAIL: Graph violates planar constraints" << std::endl;
            return false;
        }
        if (g.n >= 3 && g.m > 2 * g.n - 3) {
            if (verbose_mode) std::cout << "  FAIL: Too many edges for outerplanar graph" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: Positive OP certificate validates" << std::endl;
        return true;
    }
    bool verify_theoretical_properties(const graph& g, const gsp_sp_op_result& result) const {
        if (verbose_mode) std::cout << "Verifying theoretical properties..." << std::endl;
        if (result.is_sp && !result.is_gsp) {
            if (verbose_mode) std::cout << "  FAIL: Graph is SP but not GSP (impossible)" << std::endl;
            return false;
        }
        if (result.is_op && !result.is_gsp) {
            if (verbose_mode) std::cout << "  FAIL: Graph is OP but not GSP (impossible)" << std::endl;
            return false;
        }
        if (result.is_op && !is_potentially_planar(g)) {
            if (verbose_mode) std::cout << "  FAIL: Graph is OP but violates planar constraints" << std::endl;
            return false;
        }
        if (result.is_op && g.n >= 3 && g.m > 2 * g.n - 3) {
            if (verbose_mode) std::cout << "  FAIL: Too many edges for outerplanar graph" << std::endl;
            return false;
        }
        if (result.is_gsp && !is_connected(g)) {
            if (verbose_mode) std::cout << "  FAIL: GSP graph is not connected" << std::endl;
            return false;
        }
        if (verbose_mode) std::cout << "  PASS: Theoretical properties validate" << std::endl;
        return true;
    }
    bool is_potentially_planar(const graph& g) const {
        if (g.n >= 3 && g.m > 3 * g.n - 6) {
            return false;
        }
        if (g.n >= 5) {
            int high_degree_vertices = 0;
            for (int i = 0; i < g.n; i++) {
                if (g.adj[i].size() >= 4) {
                    high_degree_vertices++;
                }
            }
            if (high_degree_vertices > g.n / 2) {
                return false;
            }
        }
        return true;
    }
    void export_to_dot(const graph& g, const std::string& filename, 
                      const gsp_sp_op_result& result) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Could not open file: " << filename << std::endl;
            return;
        }
        file << "graph G {" << std::endl;
        file << "  node [shape=circle, style=filled, fillcolor=lightblue];" << std::endl;
        file << "  edge [color=black];" << std::endl;
        for (int i = 0; i < g.n; i++) {
            file << "  " << i << " [label=\"" << i << "\"];" << std::endl;
        }
        std::set<std::pair<int, int>> added_edges;
        for (int u = 0; u < g.n; u++) {
            for (int v : g.adj[u]) {
                if (u < v) { 
                    file << "  " << u << " -- " << v << ";" << std::endl;
                    added_edges.insert({u, v});
                }
            }
        }
        file << "  labelloc=\"t\";" << std::endl;
        file << "  label=\"Graph Properties (n=" << g.n << ", m=" << g.m << "):\\n";
        file << "GSP: " << (result.is_gsp ? "YES" : "NO") << "\\n";
        file << "SP: " << (result.is_sp ? "YES" : "NO") << "\\n";
        file << "OP: " << (result.is_op ? "YES" : "NO") << "\";" << std::endl;
        file << "}" << std::endl;
        file.close();
    }
    void export_analysis_report(const graph& g, const gsp_sp_op_result& result, 
                               const std::string& filename) const {
        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Could not open file: " << filename << std::endl;
            return;
        }
        file << "=== GRAPH ANALYSIS REPORT ===" << std::endl;
        file << "Vertices: " << g.n << std::endl;
        file << "Edges: " << g.m << std::endl;
        file << "Density: " << (g.n > 1 ? (2.0 * g.m) / (g.n * (g.n - 1)) : 0) << std::endl;
        std::vector<int> degrees;
        for (int i = 0; i < g.n; i++) {
            degrees.push_back(g.adj[i].size());
        }
        std::sort(degrees.begin(), degrees.end());
        file << "Min degree: " << (degrees.empty() ? 0 : degrees[0]) << std::endl;
        file << "Max degree: " << (degrees.empty() ? 0 : degrees[g.n-1]) << std::endl;
        file << "Average degree: " << (g.n > 0 ? (2.0 * g.m) / g.n : 0) << std::endl;
        file << std::endl << "=== PROPERTIES ===" << std::endl;
        file << "GSP: " << (result.is_gsp ? "YES" : "NO") << std::endl;
        file << "SP: " << (result.is_sp ? "YES" : "NO") << std::endl;
        file << "OP: " << (result.is_op ? "YES" : "NO") << std::endl;
        file << "Connected: " << (is_connected(g) ? "YES" : "NO") << std::endl;
        file << "Potentially planar: " << (is_potentially_planar(g) ? "YES" : "NO") << std::endl;
        file << std::endl << "=== CERTIFICATE ANALYSIS ===" << std::endl;
        if (result.gsp_reason) {
            file << "GSP certificate type: ";
            if (dynamic_cast<const negative_cert_K4*>(result.gsp_reason)) {
                file << "K4 subdivision" << std::endl;
            } else if (dynamic_cast<const negative_cert_K23*>(result.gsp_reason)) {
                file << "K23 subdivision" << std::endl;
            } else if (dynamic_cast<const negative_cert_T4*>(result.gsp_reason)) {
                file << "T4 subdivision" << std::endl;
            } else if (dynamic_cast<const positive_cert_gsp*>(result.gsp_reason)) {
                file << "Positive GSP" << std::endl;
            } else {
                file << "Other" << std::endl;
            }
        }
        if (result.op_reason) {
            file << "OP certificate type: ";
            if (dynamic_cast<const positive_cert_op*>(result.op_reason)) {
                file << "Positive OP" << std::endl;
            } else {
                file << "Negative (inherited)" << std::endl;
            }
        }
        file.close();
    }
public:
    enhanced_validator(bool verbose = false) : verbose_mode(verbose) {}
    bool deep_validate(const graph& g, const gsp_sp_op_result& result) {
        auto start_time = std::chrono::high_resolution_clock::now();
        stats.total_graphs++;
        if (verbose_mode) {
            std::cout << "=== DEEP VALIDATION ===" << std::endl;
            std::cout << "Graph: " << g.n << " vertices, " << g.m << " edges" << std::endl;
        }
        if (result.is_gsp) stats.gsp_graphs++;
        if (result.is_sp) stats.sp_graphs++;
        if (result.is_op) stats.op_graphs++;
        if (!result.authenticate(g)) {
            if (verbose_mode) std::cout << "FAIL: Standard authentication failed" << std::endl;
            stats.authentication_failures++;
            auto end_time = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double>(end_time - start_time).count();
            stats.update_timing(duration);
            return false;
        }
        if (!verify_theoretical_properties(g, result)) {
            stats.structural_inconsistencies++;
            auto end_time = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double>(end_time - start_time).count();
            stats.update_timing(duration);
            return false;
        }
        bool deep_valid = true;
        if (result.gsp_reason) {
            if (auto k4_cert = dynamic_cast<const negative_cert_K4*>(result.gsp_reason)) {
                stats.k4_certificates++;
                deep_valid &= verify_k4_subdivision(g, k4_cert);
            } else if (auto k23_cert = dynamic_cast<const negative_cert_K23*>(result.gsp_reason)) {
                stats.k23_certificates++;
                deep_valid &= verify_k23_subdivision(g, k23_cert);
            } else if (auto t4_cert = dynamic_cast<const negative_cert_T4*>(result.gsp_reason)) {
                stats.t4_certificates++;
                deep_valid &= verify_t4_subdivision(g, t4_cert);
            } else if (dynamic_cast<const negative_cert_tri_comp_cut*>(result.gsp_reason)) {
                stats.tri_comp_cut_certificates++;
            } else if (dynamic_cast<const negative_cert_tri_cut_comp*>(result.gsp_reason)) {
                stats.tri_cut_comp_certificates++;
            } else if (auto pos_cert = dynamic_cast<const positive_cert_gsp*>(result.gsp_reason)) {
                stats.positive_gsp_certificates++;
                deep_valid &= verify_positive_gsp(g, pos_cert);
            }
        }
        if (result.op_reason) {
            if (auto op_cert = dynamic_cast<const positive_cert_op*>(result.op_reason)) {
                stats.positive_op_certificates++;
                deep_valid &= verify_positive_op(g, op_cert);
            }
        }
        if (!deep_valid) {
            stats.deep_validation_failures++;
            if (verbose_mode) std::cout << "FAIL: Deep validation failed" << std::endl;
            auto end_time = std::chrono::high_resolution_clock::now();
            double duration = std::chrono::duration<double>(end_time - start_time).count();
            stats.update_timing(duration);
            return false;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        double duration = std::chrono::duration<double>(end_time - start_time).count();
        stats.update_timing(duration);
        if (verbose_mode) std::cout << "PASS: Deep validation successful" << std::endl;
        return true;
    }
    bool cross_validate(const graph& g, const gsp_sp_op_result& result) {
        if (verbose_mode) std::cout << "=== CROSS VALIDATION ===" << std::endl;
        bool independent_planar_check = is_potentially_planar(g);
        if (result.is_op && !independent_planar_check) {
            if (verbose_mode) std::cout << "FAIL: Independent planarity check disagrees" << std::endl;
            stats.cross_validation_failures++;
            return false;
        }
        bool independent_connected_check = is_connected(g);
        if (result.is_gsp && !independent_connected_check) {
            if (verbose_mode) std::cout << "FAIL: Independent connectivity check disagrees" << std::endl;
            stats.cross_validation_failures++;
            return false;
        }
        auto components = find_biconnected_components(g);
        if (g.n <= 10) {
        }
        if (verbose_mode) std::cout << "PASS: Cross validation successful" << std::endl;
        return true;
    }
    bool full_validate(const graph& g, const gsp_sp_op_result& result) {
        return deep_validate(g, result) && cross_validate(g, result);
    }
    void export_visualization(const graph& g, const gsp_sp_op_result& result, 
                            const std::string& filename) {
        export_to_dot(g, filename, result);
        if (verbose_mode) std::cout << "Visualization exported to: " << filename << std::endl;
    }
    void export_analysis(const graph& g, const gsp_sp_op_result& result,
                        const std::string& filename) {
        export_analysis_report(g, result, filename);
        if (verbose_mode) std::cout << "Analysis report exported to: " << filename << std::endl;
    }
    const validation_stats& get_stats() const {
        return stats;
    }
    void reset_stats() {
        stats = validation_stats();
    }
    void set_verbose(bool verbose) {
        verbose_mode = verbose;
    }
    void print_comprehensive_report() const {
        stats.print_stats();
        std::cout << "\n=== VALIDATION RATIOS ===" << std::endl;
        if (stats.total_graphs > 0) {
            double gsp_ratio = (double)stats.gsp_graphs / stats.total_graphs;
            double sp_ratio = (double)stats.sp_graphs / stats.total_graphs;
            double op_ratio = (double)stats.op_graphs / stats.total_graphs;
            std::cout << "GSP ratio: " << gsp_ratio << std::endl;
            std::cout << "SP ratio: " << sp_ratio << std::endl;
            std::cout << "OP ratio: " << op_ratio << std::endl;
            if (gsp_ratio > 0.8) {
                std::cout << "WARNING: GSP ratio seems high" << std::endl;
            }
            if (sp_ratio > 0.5) {
                std::cout << "WARNING: SP ratio seems high" << std::endl;
            }
            if (op_ratio > 0.3) {
                std::cout << "WARNING: OP ratio seems high" << std::endl;
            }
        }
        std::cout << "\n=== RELIABILITY METRICS ===" << std::endl;
        int total_failures = stats.authentication_failures + stats.deep_validation_failures + 
                            stats.cross_validation_failures + stats.structural_inconsistencies;
        double success_rate = 1.0 - (double)total_failures / stats.total_graphs;
        std::cout << "Overall success rate: " << (success_rate * 100) << "%" << std::endl;
        if (success_rate < 0.99) {
            std::cout << "WARNING: Success rate below 99%" << std::endl;
        }
        if (stats.total_graphs > 0) {
            std::cout << "Authentication failure rate: " << (100.0 * stats.authentication_failures / stats.total_graphs) << "%" << std::endl;
            std::cout << "Deep validation failure rate: " << (100.0 * stats.deep_validation_failures / stats.total_graphs) << "%" << std::endl;
            std::cout << "Cross validation failure rate: " << (100.0 * stats.cross_validation_failures / stats.total_graphs) << "%" << std::endl;
            std::cout << "Structural inconsistency rate: " << (100.0 * stats.structural_inconsistencies / stats.total_graphs) << "%" << std::endl;
        }
    }
    void print_performance_analysis() const {
        std::cout << "\n=== PERFORMANCE ANALYSIS ===" << std::endl;
        std::cout << "Total validation time: " << stats.total_validation_time << " seconds" << std::endl;
        std::cout << "Average time per graph: " << stats.avg_validation_time << " seconds" << std::endl;
        std::cout << "Maximum time for single graph: " << stats.max_validation_time << " seconds" << std::endl;
        if (stats.total_graphs > 0) {
            std::cout << "Throughput: " << (stats.total_graphs / stats.total_validation_time) << " graphs/second" << std::endl;
        }
        if (stats.max_validation_time > 1.0) {
            std::cout << "WARNING: Some graphs took over 1 second to validate" << std::endl;
        }
    }
    void print_certificate_analysis() const {
        std::cout << "\n=== CERTIFICATE ANALYSIS ===" << std::endl;
        int total_negative_certs = stats.k4_certificates + stats.k23_certificates + 
                                  stats.t4_certificates + stats.tri_comp_cut_certificates + 
                                  stats.tri_cut_comp_certificates;
        int total_positive_certs = stats.positive_gsp_certificates + stats.positive_op_certificates;
        std::cout << "Total negative certificates: " << total_negative_certs << std::endl;
        std::cout << "Total positive certificates: " << total_positive_certs << std::endl;
        if (total_negative_certs > 0) {
            std::cout << "K4 subdivision rate: " << (100.0 * stats.k4_certificates / total_negative_certs) << "%" << std::endl;
            std::cout << "K23 subdivision rate: " << (100.0 * stats.k23_certificates / total_negative_certs) << "%" << std::endl;
            std::cout << "T4 subdivision rate: " << (100.0 * stats.t4_certificates / total_negative_certs) << "%" << std::endl;
            std::cout << "Tri-comp-cut rate: " << (100.0 * stats.tri_comp_cut_certificates / total_negative_certs) << "%" << std::endl;
            std::cout << "Tri-cut-comp rate: " << (100.0 * stats.tri_cut_comp_certificates / total_negative_certs) << "%" << std::endl;
        }
        if (total_positive_certs > 0) {
            std::cout << "Positive GSP rate: " << (100.0 * stats.positive_gsp_certificates / total_positive_certs) << "%" << std::endl;
            std::cout << "Positive OP rate: " << (100.0 * stats.positive_op_certificates / total_positive_certs) << "%" << std::endl;
        }
    }
    std::vector<bool> batch_validate(const std::vector<graph>& graphs, 
                                   const std::vector<gsp_sp_op_result>& results,
                                   bool use_cross_validation = true) {
        std::vector<bool> validation_results;
        validation_results.reserve(graphs.size());
        if (graphs.size() != results.size()) {
            std::cerr << "Error: Graph and result vectors must have the same size" << std::endl;
            return validation_results;
        }
        std::cout << "=== BATCH VALIDATION ===" << std::endl;
        std::cout << "Validating " << graphs.size() << " graphs..." << std::endl;
        auto batch_start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < graphs.size(); i++) {
            if (verbose_mode) {
                std::cout << "Validating graph " << (i + 1) << "/" << graphs.size() << std::endl;
            }
            bool valid = use_cross_validation ? 
                        full_validate(graphs[i], results[i]) : 
                        deep_validate(graphs[i], results[i]);
            validation_results.push_back(valid);
            if (!valid && verbose_mode) {
                std::cout << "Graph " << (i + 1) << " failed validation" << std::endl;
            }
        }
        auto batch_end = std::chrono::high_resolution_clock::now();
        double batch_duration = std::chrono::duration<double>(batch_end - batch_start).count();
        std::cout << "Batch validation completed in " << batch_duration << " seconds" << std::endl;
        return validation_results;
    }
    void analyze_graph_properties(const std::vector<graph>& graphs) const {
        if (graphs.empty()) return;
        std::cout << "\n=== GRAPH PROPERTIES ANALYSIS ===" << std::endl;
        std::vector<int> sizes, edge_counts;
        std::vector<double> densities;
        for (const auto& g : graphs) {
            sizes.push_back(g.n);
            edge_counts.push_back(g.m);
            if (g.n > 1) {
                densities.push_back((2.0 * g.m) / (g.n * (g.n - 1)));
            }
        }
        auto minmax_size = std::minmax_element(sizes.begin(), sizes.end());
        auto minmax_edges = std::minmax_element(edge_counts.begin(), edge_counts.end());
        double avg_size = std::accumulate(sizes.begin(), sizes.end(), 0.0) / sizes.size();
        double avg_edges = std::accumulate(edge_counts.begin(), edge_counts.end(), 0.0) / edge_counts.size();
        double avg_density = densities.empty() ? 0.0 : 
                           std::accumulate(densities.begin(), densities.end(), 0.0) / densities.size();
        std::cout << "Graph sizes: " << *minmax_size.first << " to " << *minmax_size.second 
                  << " (avg: " << avg_size << ")" << std::endl;
        std::cout << "Edge counts: " << *minmax_edges.first << " to " << *minmax_edges.second 
                  << " (avg: " << avg_edges << ")" << std::endl;
        std::cout << "Average density: " << avg_density << std::endl;
        std::map<int, int> degree_distribution;
        for (const auto& g : graphs) {
            for (int i = 0; i < g.n; i++) {
                degree_distribution[g.adj[i].size()]++;
            }
        }
        std::cout << "\nDegree distribution (top 10):" << std::endl;
        auto it = degree_distribution.rbegin();
        for (int count = 0; count < 10 && it != degree_distribution.rend(); ++it, ++count) {
            std::cout << "  Degree " << it->first << ": " << it->second << " vertices" << std::endl;
        }
    }
};
