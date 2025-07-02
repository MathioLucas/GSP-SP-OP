#include "graph.hxx"
#include "gsp-sp-op.hxx"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void export_to_json(const graph& g, const std::string& filename, const gsp_sp_op_result* result = nullptr) {
    std::ofstream file(filename);
    file << "{\n";
    file << "  \"vertices\": " << g.adjacency_list.size() << ",\n";
    
    int edge_count = 0;
    for (const auto& adj : g.adjacency_list) {
        edge_count += adj.size();
    }
    edge_count /= 2; 
    
    file << "  \"edges\": " << edge_count << ",\n";
    file << "  \"adjacency_list\": {\n";
    
    for (size_t i = 0; i < g.adjacency_list.size(); ++i) {
        file << "    \"" << i << "\": [";
        for (size_t j = 0; j < g.adjacency_list[i].size(); ++j) {
            file << g.adjacency_list[i][j];
            if (j < g.adjacency_list[i].size() - 1) file << ", ";
        }
        file << "]";
        if (i < g.adjacency_list.size() - 1) file << ",";
        file << "\n";
    }
    file << "  }";
    
    //  classification results is posssible
    if (result) {
        file << ",\n  \"classification\": {\n";
        file << "    \"is_gsp\": " << (result->is_gsp ? "true" : "false") << ",\n";
        file << "    \"is_sp\": " << (result->is_sp ? "true" : "false") << ",\n";
        file << "    \"is_op\": " << (result->is_op ? "true" : "false") << "\n";
        file << "  }";
    }
    
    file << "\n}\n";
    std::cout << "Exported to JSON: " << filename << std::endl;
}

void export_to_graphml(const graph& g, const std::string& filename, const gsp_sp_op_result* result = nullptr) {
    std::ofstream file(filename);
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    file << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">\n";
    file << "  <graph id=\"G\" edgedefault=\"undirected\">\n";
    
    //  vertices
    for (size_t i = 0; i < g.adjacency_list.size(); ++i) {
        file << "    <node id=\"" << i << "\"/>\n";
    }
    
    //  edges (TO avoid duplicates)
    std::set<std::pair<int, int>> added_edges;
    for (size_t i = 0; i < g.adjacency_list.size(); ++i) {
        for (int neighbor : g.adjacency_list[i]) {
            int u = static_cast<int>(i);
            int v = neighbor;
            if (u < v && added_edges.find({u, v}) == added_edges.end()) {
                file << "    <edge source=\"" << u << "\" target=\"" << v << "\"/>\n";
                added_edges.insert({u, v});
            }
        }
    }
    
    file << "  </graph>\n";
    file << "</graphml>\n";
    std::cout << "Exported to GraphML: " << filename << std::endl;
}

void export_adjacency_matrix(const graph& g, const std::string& filename) {
    std::ofstream file(filename);
    int n = g.adjacency_list.size();
    
    //  adjacency matrix
    std::vector<std::vector<int>> matrix(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int neighbor : g.adjacency_list[i]) {
            matrix[i][neighbor] = 1;
        }
    }
    

    file << "# Adjacency Matrix (" << n << "x" << n << ")\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            file << matrix[i][j];
            if (j < n - 1) file << " ";
        }
        file << "\n";
    }
    std::cout << "Exported adjacency matrix: " << filename << std::endl;
}

void export_statistics(const graph& g, const std::string& filename, const gsp_sp_op_result* result = nullptr) {
    std::ofstream file(filename);
    int n = g.adjacency_list.size();
    
    
    int edge_count = 0;
    for (const auto& adj : g.adjacency_list) {
        edge_count += adj.size();
    }
    edge_count /= 2;
    
    
    std::vector<int> degrees;
    for (const auto& adj : g.adjacency_list) {
        degrees.push_back(adj.size());
    }
    
    int max_degree = *std::max_element(degrees.begin(), degrees.end());
    int min_degree = *std::min_element(degrees.begin(), degrees.end());
    double avg_degree = 2.0 * edge_count / n;
    
    file << "Graph Statistics\n";
    file << "================\n";
    file << "Vertices: " << n << "\n";
    file << "Edges: " << edge_count << "\n";
    file << "Density: " << std::fixed << std::setprecision(4) << (2.0 * edge_count) / (n * (n - 1)) << "\n";
    file << "Average Degree: " << std::fixed << std::setprecision(2) << avg_degree << "\n";
    file << "Max Degree: " << max_degree << "\n";
    file << "Min Degree: " << min_degree << "\n";
    
    if (result) {
        file << "\nClassification Results\n";
        file << "=====================\n";
        file << "GSP (Generalized Series-Parallel): " << (result->is_gsp ? "YES" : "NO") << "\n";
        file << "SP (Series-Parallel): " << (result->is_sp ? "YES" : "NO") << "\n";
        file << "OP (Outerplanar): " << (result->is_op ? "YES" : "NO") << "\n";
        
        file << "\nSubclass: ";
        if (result->is_sp && result->is_op) {
            file << "Both SP and OP\n";
        } else if (result->is_sp && !result->is_op) {
            file << "SP but not OP\n";
        } else if (!result->is_sp && result->is_op) {
            file << "OP but not SP\n";
        } else if (result->is_gsp) {
            file << "GSP but neither SP nor OP\n";
        } else {
            file << "Not GSP, SP, or OP\n";
        }
    }
    
    std::cout << "Exported statistics: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <input.txt> <format>\n";
        std::cout << "Formats: json, graphml, adjmatrix, stats\n";
        std::cout << "Example: " << argv[0] << " testcases/biconnected/example.txt json\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string format = argv[2];
    
    
    std::ifstream file(input_file);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << input_file << std::endl;
        return 1;
    }
    
    graph g;
    file >> g;
    file.close();
    
    std::cout << "Loaded graph with " << g.adjacency_list.size() << " vertices\n";
    
   
    gsp_sp_op_result result = GSP_SP_OP(g);
    bool authenticated = result.authenticate(g);
    
    if (!authenticated) {
        std::cerr << "Warning: Certificate authentication failed!\n";
    }
    
  
    std::string base_name = input_file.substr(0, input_file.find_last_of('.'));
    std::string output_file;
    
    if (format == "json") {
        output_file = base_name + ".json";
        export_to_json(g, output_file, &result);
    } else if (format == "graphml") {
        output_file = base_name + ".graphml";
        export_to_graphml(g, output_file, &result);
    } else if (format == "adjmatrix") {
        output_file = base_name + "_matrix.txt";
        export_adjacency_matrix(g, output_file);
    } else if (format == "stats") {
        output_file = base_name + "_stats.txt";
        export_statistics(g, output_file, &result);
    } else {
        std::cerr << "Unknown format: " << format << std::endl;
        std::cerr << "Supported formats: json, graphml, adjmatrix, stats\n";
        return 1;
    }
    
    return 0;
}
