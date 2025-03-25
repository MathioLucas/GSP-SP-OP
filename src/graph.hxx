// this file contains the definition of a graph and simple graph operations

#ifndef __GRAPH_HXX__
#define __GRAPH_HXX__

#include "logging.hxx"
#include <vector>
#include <istream>
#include <ostream>

using edge_t = std::pair<int, int>;

struct graph {
	int n; // graph order
	int e; // graph size
	std::vector<std::vector<int> > adjLists; // graph adjacency lists; adjLists[i] is a vector of all the vertices vertex i is adjacent to

	bool adjacent(int e1, int e2) const { // check if two vertices are adjacent (O|deg(v)| time)
		for (int v : adjLists[e1]) {
			if (v == e2) return true;
		}

		return false;
	}

	void add_edge(int e1, int e2) {
		adjLists[e1].push_back(e2);
		adjLists[e2].push_back(e1);
	}

	void reserve(graph const& other) {
		for (int i = 0; i < other.n; i++) {
			adjLists.emplace_back();
			adjLists[i].reserve(other.adjLists[i].size());
		}
	}

	void output_adj_list(int v, std::ostream& os) const {
		os << "vertex " << v << " adjacencies: ";

		for (int v2 : adjLists[v]) {
			os << v2 << " ";
		}

		os << "\n";
	}
};

std::istream& operator>>(std::istream& is, graph& g) { // read a graph from an input stream (e.g. file)
	g = graph{};
	is >> g.n >> g.e;

	g.adjLists.reserve(g.n);

	for (int i = 0; i < g.n; i++) { // construct empty adjacency lists for each vertex
		g.adjLists.emplace_back();
	}

	for (int i = 0; i < g.e; i++) { // insert the edges into the graph
		int endpoint1, endpoint2;
		is >> endpoint1 >> endpoint2;
		g.add_edge(endpoint1, endpoint2);
	}


	for (std::vector<int> list : g.adjLists) { // free up the extra memory allocated by the vectors from being not at capacity when we inserted the edges (O(|E|) time, according to cppreference)
		list.shrink_to_fit();
	}

	return is;
}

std::ostream& operator<<(std::ostream& os, graph const& g) { // output a graph (for debugging purposes)
	os << "Graph with " << g.n << " vertices and " << g.e << " edges:\n";

	for (int i = 0; i < g.n; i++) {
		g.output_adj_list(i, os);
	}

	return os;
}

#endif