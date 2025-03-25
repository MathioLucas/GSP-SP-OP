// this file contains definitions of the possible certificates the implmentation can output, structures to store the results of its execution, and the relevant authentication algorithms

#ifndef __GSP_SP_OP_CERTS_HXX__
#define __GSP_SP_OP_CERTS_HXX__

#include "graph.hxx"
#include "logging.hxx"
#include "sp-tree.hxx"
#include "radix_sort.hxx"
#include <memory>
#include <vector>
#include <stack>

// ---------------- auxiliary functions ----------------

bool trace_path(int end1, int end2, std::vector<edge_t> const& path, graph const& g, std::vector<bool>& seen) { // trace a path, to ensure it is between end1 and end2 (either direction), all its edges are in g, and none of its internal vertices are in seen
	#ifdef __VERBOSE_LOGGING__
	for (edge_t edge : path) {
		V_LOG("(" << edge.first << ", " << edge.second << ") ")
	}
	#endif
	
	N_LOG("\n")
	if (path.size() == 0) { // check trivial path
		L_LOG("====== AUTH FAILED: no edges in path ======\n")
		return false;
	}

	if (path[0].first == end2) { // check if path is reversed, swap ends if so
		int tmp = end2;
		end2 = end1;
		end1 = tmp;
	}

	if (path[0].first != end1) { // neither endpoint matches start
		L_LOG("====== AUTH FAILED: start of path does not match either endpoint ======\n")
		return false;
	}

	if (path.back().second != end2) { // one endpoint matches, but the other doesn't
		L_LOG("====== AUTH FAILED: end of path does not match second endpoint ======\n")
		return false;
	}

	seen[end1] = true;
	int prev_v = end1;
	for (edge_t edge : path) {
		if (!g.adjacent(edge.first, edge.second)) { // check edges are in graph
			L_LOG("====== AUTH FAILED: edge (" << edge.first << ", " << edge.second << ") does not exist in graph ======\n")
			return false;
		}

		if (prev_v != edge.first) { // check edges join together to form path
				L_LOG("====== AUTH FAILED: edge (" << edge.first << ", " << edge.second << ") is not incident on the previous edge ======\n")
				return false;
		}

		prev_v = edge.second;

		if (seen[edge.second]) { // check path internally disjoint with prior paths and itself
			L_LOG("====== AUTH FAILED: duplicated vertex " << edge.second << " ======\n")
			return false;
		}

		seen[edge.second] = true;
	}

	N_LOG("path good\n")
	seen[end1] = false; // un-seen the terminating vertices (they are allowed to repeat, just not within a path)
	seen[end2] = false;
	return true;
}

int num_comps_after_removal(graph const& g, int v) { // get number of components of g after v is removed
	int retval = 0;

	std::vector<bool> seen((size_t)(g.n), false);

	for (int i = 0; i < g.n; i++) { // start DFS from every vertex, unless we've already seen it, and the number of DFS we start is the number of connected components. the correctness/time complexity are obvious by the nature of DFS
		if (seen[i] || i == v) continue;
		retval++;

		std::stack<int> dfs;
		dfs.emplace(i);

		while (!dfs.empty()) {
			int w = dfs.top();
			dfs.pop();
			seen[w] = true;

			for (int u : g.adjLists[w]) {
				if (!seen[u] && u != v) {
					dfs.emplace(u);
				}
			}
		}
	}

	return retval;
}

bool is_cut_vertex(graph const& g, int v) {
	if (num_comps_after_removal(g, v) <= 1) {
		L_LOG("\n====== AUTH FAILED: " << v << " not a cut vertex ======\n\n")
		return false;
	}

	N_LOG("yes\n")
	return true;
	
}

// ---------------- cert definitions ----------------

struct certificate {
	bool verified = false;
	virtual bool authenticate(graph const& g) = 0;
	virtual ~certificate() {}
};

struct negative_cert_K4 : certificate { // K4 subdivision, showing a graph is not GSP, OP, or SP
	int a, b, c, d; // four vertices
	std::vector<edge_t> ab, ac, ad, bc, bd, cd; // six pairwise internally disjoint paths between them

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate

		L_LOG("====== AUTHENTICATE K4: terminating vertices a: " << a << ", b: " << b << ", c: " << c << ", d: " << d << " ======\n")
		if (a == b || b == c || c == d || d == a || a == c || b == d) { // vertices must be distinct
			L_LOG("====== AUTH FAILED: terminating vertices non-distinct ======\n\n")
			return false;
		}
		std::vector<bool> seen((size_t)(g.n), false);

		N_LOG("verify ab: ") // paths must be pairwise internally disjoint, exist in graph, and terminate at endpoints
		if (!trace_path(a, b, ab, g, seen)) return false;
		N_LOG("verify ac: ")
		if (!trace_path(a, c, ac, g, seen)) return false;
		N_LOG("verify ad: ")
		if (!trace_path(a, d, ad, g, seen)) return false;
		N_LOG("verify bc: ")
		if (!trace_path(b, c, bc, g, seen)) return false;
		N_LOG("verify bd: ")
		if (!trace_path(b, d, bd, g, seen)) return false;
		N_LOG("verify cd: ")
		if (!trace_path(c, d, cd, g, seen)) return false;

		L_LOG("====== AUTH SUCCESS ======\n\n")
		verified = true;
		return true;
	}
};

struct negative_cert_K23 : certificate { // K_(2,3) subdivision, showing a graph is not OP
	int a, b; // two vertices (making up the independent set of size 2 in K_(2,3))
	std::vector<edge_t> one, two, three; // three pairwise internally disjoint paths of length at least two between them

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate

		L_LOG("====== AUTHENTICATE K23: terminating vertices a: " << a << ", b: " << b << " ======\n")

		if (a == b) { // vertices must be distinct
			L_LOG("====== AUTH FAILED: terminating vertices non-distinct ======\n\n")
			return false;
		}

		std::vector<bool> seen((size_t)(g.n), false);

		N_LOG("verify path one: ") // paths must be pairwise internally disjoint, exist in graph, and terminate at endpoints
		if (!trace_path(a, b, one, g, seen)) return false;
		if (one.size() < 2) { // paths need to have at least one internal vertex to complete K23
			L_LOG("\n====== AUTH FAILED: path one has no internal vertex ======\n\n")
			return false;
		}

		N_LOG("verify path two: ")
		if (!trace_path(a, b, two, g, seen)) return false;
		if (two.size() < 2) {
			L_LOG("\n====== AUTH FAILED: path two has no internal vertex ======\n\n")
			return false;
		}

		N_LOG("verify path three: ")
		if (!trace_path(a, b, three, g, seen)) return false;
		if (three.size() < 2) {
			L_LOG("\n====== AUTH FAILED: path three has no internal vertex ======\n\n")
			return false;
		}

		L_LOG("====== AUTH SUCCESS ======\n\n")

		verified = true;
		return true;
	}
};

struct negative_cert_T4 : certificate { // a theta_(4) subdivision with the top and bottom of the "theta" being two cut vertices of a graph, showing the graph is not SP (but may still be GSP)
	int c1, c2, a, b; // four vertices, with c1 and c2 cut vertices
	std::vector<edge_t> c1a, c1b, c2a, c2b, ab; // five pairwise internally disjoint paths between them (the path between c1 and c2 is not part of theta_(4))

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate
		L_LOG("====== AUTHENTICATE T4: terminating vertices a: " << a << ", b: " << b << ", c1: " << c1 << ", c2: " << c2 << " ======\n")

		if (a == b || a == c1 || a == c2 || b == c1 || b == c2 || c1 == c2) { // vertices must be distinct
			L_LOG("====== AUTH FAILED: terminating vertices non-distinct ======\n\n")
			return false;
		}

		N_LOG("verify c1 cut vertex: ")
		if (!is_cut_vertex(g, c1)) return false;
		N_LOG("verify c2 cut vertex: ")
		if (!is_cut_vertex(g, c2)) return false;


		std::vector<bool> seen((size_t)(g.n), false);
		N_LOG("verify path c1a: ") // paths must be pairwise internally disjoint, exist in graph, and terminate at endpoints
		if (!trace_path(c1, a, c1a, g, seen)) return false;
		N_LOG("verify path c2a: ")
		if (!trace_path(c2, a, c2a, g, seen)) return false;
		N_LOG("verify path ab: ")
		if (!trace_path(a, b, ab, g, seen)) return false;
		N_LOG("verify path c1b: ")
		if (!trace_path(c1, b, c1b, g, seen)) return false;
		N_LOG("verify path c2b: ")
		if (!trace_path(c2, b, c2b, g, seen)) return false;

		L_LOG("====== AUTH SUCCESS ======\n\n")

		verified = true;
		return true;
	}
};

struct negative_cert_tri_comp_cut : certificate { // a cut vertex that divides a graph into three or more components, showing the graph is not SP (but may still be GSP)
	int v; // the cut vertex

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate
		L_LOG("====== AUTHENTICATE THREE-COMPONENT CUT VERTEX: " << v << " ======\n")

		int comps = num_comps_after_removal(g, v);

		if (comps < 3) { // needs to have at least 3 components after v is removed
			L_LOG("====== AUTH FAILED: vertex " << v << " only splits graph into " << comps << " components ======\n\n")
			return false;
		}

		N_LOG(comps << " comps after removal\n")
		L_LOG("====== AUTH SUCCESS ======\n\n")

		verified = true;
		return true;
	}
};

struct negative_cert_tri_cut_comp : certificate { // three cut vertices which are part of a single biconnected component of a graph, showing the graph is not SP (but may still be GSP)
	int c1, c2, c3; // the three cut vertices

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate
		L_LOG("====== AUTHENTICATE BICOMP WITH THREE CUT VERTICES: cut vertices " << c1 << ", " << c2 << ", " << c3 << " ======\n")
		N_LOG("verify c1 cut vertex: ") // they all need to be cut vertices
		if (!is_cut_vertex(g, c1)) return false;
		N_LOG("verify c2 cut vertex: ")
		if (!is_cut_vertex(g, c2)) return false;
		N_LOG("verify c3 cut vertex: ")
		if (!is_cut_vertex(g, c2)) return false;

		std::vector<int> dfs_no((size_t)(g.n), 0); // this is a carbon copy of the algorithm for finding bicomps; I don't see an easier way of doing this so this will have to suffice for authentication
		std::vector<int> parent((size_t)(g.n)); 
		std::vector<int> low((size_t)(g.n));
		int cut_verts[3] = {c1, c2, c3}; // for convenience, I really should've stored all these as an array to begin with but oh well

		std::stack<edge_t> comp_edges;

		std::stack<std::pair<int, int>> dfs;

		dfs.emplace(0, 0);

		dfs_no[0] = 1;
		low[0] = 1;
		parent[0] = -1;
		int curr_dfs = 2;

		while (!dfs.empty()) {
			std::pair<int, int> p = dfs.top();
			int w = p.first;
			int u = g.adjLists[p.first][p.second];

			if (dfs_no[u] == 0) { // recurse if we haven't seen u
				dfs.push(std::pair{u, 0});
				comp_edges.emplace(w, u);
				parent[u] = w;
				dfs_no[u] = curr_dfs++;
				low[u] = dfs_no[u]; // initialize low[u]; this will never actually cause a bicomp to be detected (since dfs_no[w] < dfs_no[u] by the nature of DFS and dfs_no[u] = low[u] by this statement implies low[u] >/= dfs_no[w] if there are no outgoing back edges of a descendant of u) 
									// we do it anyway, though, so that low[u] is always defined
				continue;
			} // otherwise, we have returned from the recursive call where we visited u

			if (parent[u] == w) { // w->u is a tree edge
				if (low[u] >= dfs_no[w]) { // detect a biconnected component; if the lowest back edge coming out of a descendant of the child isn't a proper ancestor of the parent then no back edge coming out of a descendant of the child is and w must be a cut vertex
					bool seen[3] = {false, false, false};
					edge_t e;
					do {
						e = comp_edges.top();

						for (int i = 0; i < 3; i++) {
							if (e.first == cut_verts[i] || e.second == cut_verts[i]) seen[i] = true;
						}
						comp_edges.pop();
					} while (e != edge_t{w, u});

					if (seen[0] && seen[1] && seen[2]) { // if this bicomp contains all three of the vertices, we're good
						N_LOG("vertices belong to one biconnected component...\n")	
						L_LOG("====== AUTH SUCCESS ======\n\n")

						verified = true;
						return true;
					}
				}

				if (low[u] < low[w]) low[w] = low[u]; // update low, if this child has a lower value then it follows that w also has at most that value, since the descendants of a child are also descendants of its parent
													  // after processing all the children of w, we will have correctly computed low[w] (by induction on the height of w in the DFS tree)
			} else if (dfs_no[u] < dfs_no[w] && u != parent[w]) { // u ↶ w is a back edge
				comp_edges.emplace(w, u);
				if (dfs_no[u] < low[w]) low[w] = dfs_no[u]; // update low, if the outgoing back edge has a lower DFS number at its sink then it follows that w has at most that value as its low, since w is trivially a descendant of w
			}

			if ((size_t)(++dfs.top().second) >= g.adjLists[p.first].size()) { // we are done processing the adjacency list of w and are backing up to the parent
				dfs.pop(); // end the recursive call
			}
		}

		L_LOG("====== AUTH FAILED: bicomp does not contain the three cut vertices ======\n\n")
		return false;
	}
};

struct positive_cert_gsp : certificate { // SP decomposition tree showing a graph is GSP or SP
	sp_tree decomposition;
	bool is_sp; // mark true if this is an SP tree (no dangling compositions allowed)

	bool authenticate(graph const& g) override {
		if (verified) return true; // don't redo the work if we've already verified this certificate

		std::vector<int> n_src((size_t)(g.n), 0); 	  	 // n_src[i] is the number of SP subgraphs we have with source at vertex i
		std::vector<int> n_sink((size_t)(g.n), 0); 	  	 // n_sink[i] is the number of SP subgraphs we have with sink at vertex i	
		std::vector<bool> no_edge((size_t)(g.n), false); // no_edge[i] is whether we've already merged vertex i into an SP subgraph, preventing any other edges incident on it
		bool swap = false; 	  // keeps track of whether source and sink of an SP subgraph should be swapped (in an antiparallel composition)

		graph g2{}; // graph we will construct
		g2.reserve(g); // reserve edge space
		std::stack<std::pair<sp_tree_node *, int>> hist; // my implementation performs an iterative traversal of the SP tree to prevent stack overflow on exceptionally large graphs
														 // first element of each pair is the node in the recursive stack, second is our progress in pushing its children onto the stack (needed so I can restore swap after backtracking into an antiparallel node from its right child and perform the traversal in post-order)
		L_LOG("====== AUTHENTICATE " << (is_sp ? "SP" : "GSP") << " DECOMPOSITION TREE ======\n")
		if (!decomposition.root) { // the tree can't just be all a dream
			L_LOG("====== AUTH FAILED: decomposition tree does not exist ======\n\n")
			return false;
		}

		hist.emplace(decomposition.root, 0); // initialize traversal

		while (!hist.empty()) {
			sp_tree_node * curr = hist.top().first;
			V_LOG("traversal: " << *curr << ", phase: " << hist.top().second << "\n")
			int source = (swap ? curr->sink : curr->source);
			int sink = (swap ? curr->source : curr->sink);

			if (hist.top().second == 0) { // first time looking at this node
				if (!(curr->l) || !(curr->r)) { // if this is a leaf node (has one or more null children), visit it and pop it off the stack
					if (curr->l || curr->r) { // can't have one child on a decomp tree
						L_LOG("====== AUTH FAILED: node " << *curr << " malformed (one child) ======\n\n")
						return false;
					}

					if (curr->comp != c_type::edge) { // can't have non-edge leaf
						L_LOG("====== AUTH FAILED: node " << *curr << " malformed (leaf, but not an edge) ======\n\n")
						return false;
					}

					if (no_edge[source] || no_edge[sink]) { // can't have this edge incident to a vertex which has already been merged into an SP subgraph
						L_LOG("====== AUTH FAILED: edge node " << *curr << " is incident on an vertex already merged into a series SP subgraph ======\n\n")
						return false;
					}

					g2.add_edge(source, sink); // add this edge to our adj lists
					n_src[source]++; // update SP subgraph count
					n_sink[sink]++;
					hist.pop(); // finish visit
				} else { // internal node, don't visit it yet and push its children so the traversal is post-order and mark it so we can visit it on the way back
					if (curr->comp == c_type::antiparallel) swap = !swap; // if we're entering a marked (antiparallel) node's right child, flip the swap switch so we know to swap the source and sink of all descendants
																		  // the order of the trees is important; the next one we'll look at is the right subtree, which is the swapped one, the left one isn't swapped
					hist.top().second++; // increment progress in processing this node
					hist.emplace(curr->r, 0);
				}
			} else if (hist.top().second == 1) { // push left child the second time we see a node
				if (curr->comp == c_type::antiparallel) swap = !swap; // if we entered an antiparallel node from its right child, restore the swap switch to its original value
				hist.top().second++; // increment progress in processing this node
				hist.emplace(curr->l, 0);
			} else { // once processed both children, we can process the parent node
				int lsource = (swap ? curr->r->sink : curr->l->source); // swap the necessary sources and sinks if required
				int lsink = (swap ? curr->r->source : curr->l->sink);   // note we also need to swap left and right children iff swap = true, otherwise e.g {aXb}aSc{bXc} would become {bXa}cSa{cXb}, which is malformed
				int rsource = (swap ? curr->l->sink : curr->r->source);
				int rsink = (swap ? curr->l->source : curr->r->sink);

				switch (curr->comp) {
					case c_type::edge: // if it's an edge node, it's malformed; internal nodes can't be edges
						L_LOG("====== AUTH FAILED: node " << *curr << " malformed (edge, but internal) ======\n\n")
						return false;
					case c_type::series:
						if (lsource != source || rsink != sink || lsink != rsource) {
							L_LOG("====== AUTH FAILED: node " << *curr << " malformed (series children source/sink mismatch) ======\n\n")
							return false;
						}

						if (n_src[lsink] != 1 || n_sink[lsink] != 1) { // can't have pending SP subgraphs at this vertex, otherwise we won't be able to merge them
							L_LOG("====== AUTH FAILED: series node " << *curr << " has incident edges on its middle vertex " << lsink << " which cannot be merged into it ======\n\n")
							return false;
						}

						V_LOG("BLOCKING: " << lsink << "\n")
						no_edge[lsink] = true; // the middle vertex can't have any edges merged into it anymore
						n_src[lsink]--; // update SP subgraph counts
						n_sink[lsink]--;
						break;
					case c_type::parallel:
						if (lsource != source || rsource != source || lsink != sink || rsink != sink) {
							L_LOG("====== AUTH FAILED: node " << *curr << " malformed (parallel children source/sink mismatch) ======\n\n")
							return false;
						}

						n_src[source]--; // update SP subgraph counts
						n_sink[sink]--;
						break;
					case c_type::antiparallel:
						if (swap) { // if we're swapped it's the left child who is inverted instead of the right one
							if (lsource != sink || rsource != source || lsink != source || rsink != sink) {
								L_LOG("====== AUTH FAILED: node " << *curr << " malformed (antiparallel children source/sink mismatch) ======\n\n")
								return false;
							}
						} else {
							if (lsource != source || rsource != sink || lsink != sink || rsink != source) {
								L_LOG("====== AUTH FAILED: node " << *curr << " malformed (antiparallel children source/sink mismatch) ======\n\n")
								return false;
							}
						}

						n_src[source]--; // update SP subgraph counts (the right child was swapped into having the correct source/sink orientation, so this step is identical to the parallel one)
						n_sink[sink]--;
						break;
					case c_type::dangling:
						if (is_sp) { // can't have dangling compositions in an SP decomposition, just a GSP one
							L_LOG("====== AUTH FAILED: illegal dangling composition in SP decomposition tree ======\n\n")
							return false;
						}

						if (swap) { // if we're swapped it's the left child who is dangling instead of the right one
									// the swap switch sort of messes with dangling compositions; if we're swapped and a dangling happens the dangling will be connected to the wrong edge 
									// for example, if there is (v, w)S((w, u)D(w, w')), it will become after a swap ((w', w)D(u, w))S(w, v), after which the dangling should be attached to (w, v) and not be swapped but it remains attached to (u, w) and swapped)
							if (rsource != source || rsink != sink || lsink != sink) {
								L_LOG("====== AUTH FAILED: node " << *curr << " malformed (dangling children source/sink mismatch) ======\n\n")
								return false;
							}
						} else {
							if (lsource != source || lsink != sink || rsource != source) {
								L_LOG("====== AUTH FAILED: node " << *curr << " malformed (dangling children source/sink mismatch) ======\n\n")
								return false;
							}
						}

						if (swap) {
							n_src[lsource]--; // update SP subgraph counts (we need special handling if the swap switch is active per the reasoning above; the left child is dangling, and its source and sink will be erroneously switched around)
							n_sink[sink]--;
						} else {
							n_src[source]--; // update SP subgraph counts
							n_sink[rsink]--; // this instruction is not in the paper; the paper only suggests decrementing the source count for the source of the parent node
											 // the paper also suggests keeping track of whether every vertex is a cut vertex or a sink of an SP subgraph of a biconnected component, so the expected final sink counts can be determined (the sink count should be 1 iff the vertex is a sink and not a cut vertex)
											 // doing all this is a lot of other stuff to keep track of in our certificate, but this single instruction in fact avoids all of this work, and makes it so you just need to check, given a root node of the decomp tree rXt, if the n_src for r is 1 and the n_sink for t is 1 (and of course that all the other source/sink numbers are 0)
											 // the correctness should be obvious enough by induction on the number of biconnected components processed in the queue produced by GSP-SP-OP
						}

						break;
				}

				hist.pop(); // we're done with this node (pushed both children and processed parent)
			}
		}

		N_LOG("decomposition tree well-formed...\n")
		n_src[decomposition.root->source]--; // uncount the source and sink of the final SP subgraph at the root of the tree; after this point n_src and n_sink should be 0 everywhere
		n_sink[decomposition.root->sink]--;

		bool failed = false;
		for (int i = 0; i < g.n; i++) {
			if (n_src[i] != 0) { // can't have other SP subgraphs than the root one which are part of the decomposition tree
				N_LOG("OH NO: disconnected SP subgraph sourced at vertex " << i << "\n")
				failed = true;
			}
			
			if (n_sink[i] != 0) {
				N_LOG("OH NO: disconnected SP subgraph sinked at vertex " << i << "\n")
				failed = true;
			}
		}

		if (failed) {
			L_LOG("====== AUTH FAILED: additional disconnected SP subgraphs are part of the decomposition tree ======\n\n")
			return false;
		}

		N_LOG("decomposition tree connected...\n")

		for (size_t i = 0; i < g2.adjLists.size(); i++) { // radix sort adjacency lists for comparison
			std::vector<int> l1 = g.adjLists[i];
			radix_sort(l1);
			radix_sort(g2.adjLists[i]);
			if (l1 != g2.adjLists[i]) { // the adjacency lists must be the same
				L_LOG("====== AUTH FAILED: vertex " << i << " of G does not have the same adjacency list as the one produced by the decomposition tree ======\n\n")

				#ifdef __LOGGING__
				N_LOG("ORIGINAL GRAPH: ")
				g2.output_adj_list(i, std::cout);
				N_LOG("PRODUCED GRAPH: ")
				g2.output_adj_list(i, std::cout);
				#endif

				L_LOG("======================================================================\n\n")
				return false;				
			}
		}

		N_LOG("decomposition tree produces graph identical to G...\n")
		L_LOG("====== AUTH SUCCESS ======\n\n")

		verified = true;
		return true;
	}
};

struct positive_cert_op : certificate { // external boundary of an outerplanar embedding showing a graph is OP
	std::vector<std::vector<edge_t>> boundaries; // boundaries[i] is the external boundary of the ith biconnected component (in the order they were produced, forming a rooted tree of biconnected components), and the union of all the boundaries is the exterior boundary of G

	bool authenticate(graph const& g) override { // NOTE: the paper gives no guidance on how to authenticate the exterior boundaries for non-biconnected outerplanar graphs, so I have filled in the details a bit here
												 // I authenticate each biconnected component's exterior boundary separately using the method given in the paper, and keep track of a few additional things to ensure they are actually biconnected components (i.e. there are no edges for which there is no component both endpoints belong to, and there are no cycles of "biconnected components")
												 // the algorithm I use is a bit complicated, and I have not proven its correctness or time complexity, but I can write up a proof of correctness if you want (it should be relatively obvious, at least around as obvious as the other things marked as "obvious" and left unproven in the paper)

		if (verified) return true; // don't redo the work if we've already verified this certificate

		std::vector<int> comp_parent(boundaries.size(), -1);   // comp_parent[i] is the parent of the ith biconnected component in the tree of biconnected components
									   					   	   // I use this to ensure there are no errant edges between one biconnected component and another
		std::vector<int> root_vertices(boundaries.size(), -1); // root_vertices[i] is the root vertex of the ith biconnected component (as defined in the paper)
		std::vector<int> component((size_t)(g.n), -1);	  	   // component[i] is the biconnected component of greatest height a vertex belongs to
		int cycle_adjs[g.n][2]; 						   	   // cycle_adjs[i] are the two adjacencies of vertex i in the cycle created by the exterior boundary of a bicomponent
														       // we need this to check that the exterior boundary forms a cycle and walk the DFS path in O(|E|) time
		int vert_count = 0;								       // we use vert_count to count the number of vertices, to confirm there are no vertices in G left untouched

		std::vector<std::stack<edge_t>> vertex_stacks((size_t)(g.n)); // we do a (greatly, GREATLY simplified since there's no need to build decomposition trees and every ear other than the big cycle is trivial) version of the algorithm for K4 detection
		std::vector<int> dfs_no((size_t)(g.n + 1), 0); 				  // we don't need to worry about detecting K23, since all the ears are trivial and the only outerplanar violation left is a 3.5(c) one, which is actually just a 3.4(b) one
													   				  // all we need to store on the vertex stacks is the .tail (since the .SP is a single back edge with sink equal to the stack vertex)
													   				  // we can keep track of seq with just one int for its source, since we don't need to generate a decomp sequence and we know where the sink of seq is
													   				  // we don't need to keep track of any ear[w], since the back edge corresponding to the ear of every tree edge is just the back edge of the big cycle


		for (int i = 0; i < g.n; i++) {
			cycle_adjs[i][0] = -1;
			cycle_adjs[i][1] = -1;
		}

		dfs_no[g.n] = g.n;

		L_LOG("====== AUTHENTICATE EXT BOUNDARY: " << boundaries.size() << " bicomp" << (boundaries.size() != 1 ? "s" : "") << " ======\n")

		for (int cnum = 0; (size_t)(cnum) < boundaries.size(); cnum++) {
			if (boundaries[cnum].size() == 0) { // component needs to exist
				L_LOG("====== AUTH FAILED: bicomp " << cnum << " has no edges ======\n\n")
				return false;
			}

			N_LOG("bicomp " << cnum << " (" << boundaries[cnum].size() << " edges in boundary)")
			V_LOG(":")
			N_LOG("\n")
			for (edge_t e : boundaries[cnum]) { // step 1: compute component[i] for every vertex, verify that there are no cycles of "biconnected components" (which would result in those components not actually being maximal)
				V_LOG("(" << e.first << ", " << e.second << ") ")

				if (component[e.first] != -1 && component[e.first] != cnum) { // if this vertex already has a component it's in that isn't this component
					if (root_vertices[component[e.first]] == -1) { // if the pre-existing component does not have a root vertex, we have found the root vertex of that component
						root_vertices[component[e.first]] = e.first;
					} else if (root_vertices[component[e.first]] != e.first) { // if the pre-existing component already has a root vertex, it needs to be this vertex or there's an illegal overlap between components
					 														   // note that this excludes a cycle of biconnected components, since if one existed there would be no order in which we could process the components that doesn't result in an illegal overlap (the last one to be processed in a prospective cycle would always screw us over)
						L_LOG("\n====== AUTH FAILED: vertex " << e.first << " overlaps between bicomps " << cnum << " and " << component[e.first] << ", but not at root ======\n\n")
						return false;
					}
				}

				if (component[e.second] != -1 && component[e.second] != cnum) { // the same, but with the other endpoint
					if (root_vertices[component[e.second]] == -1) {
						root_vertices[component[e.second]] = e.second;
					} else if (root_vertices[component[e.second]] != e.second) {
						L_LOG("\n====== AUTH FAILED: vertex " << e.second << " overlaps between bicomps " << cnum << " and " << component[e.second] << ", but not at root ======\n\n")
						return false;
					}
				}

				component[e.first] = cnum; // the endpoints of this edge are in the component
				component[e.second] = cnum;
			}

			V_LOG("\n")
		}

		#ifdef __LOGGING__
			if (boundaries.size() > 1) {
				N_LOG("bicomps have no cycles...\n")
			}
		#endif

		for (int cnum = 0; (size_t)(cnum) < boundaries.size() - 1; cnum++) { // step 2: compute component parents, verify that the components form a rooted tree (i.e. every component except for the last one has a parent, and there are no cycles which we've already checked)
			if (root_vertices[cnum] == -1) { // component needs to have a root (if it does have a root, that root is guaranteed to belong to a different component at that point, which is the parent component)
				L_LOG("====== AUTH FAILED: bicomp " << cnum << " has no root, but is not the final bicomp ======\n\n")
				return false;
			}

			comp_parent[cnum] = component[root_vertices[cnum]]; // this component's parent is the last component that overwrote its root
		}

		root_vertices[boundaries.size() - 1] = boundaries.back()[0].first; // arbitrarily choose the root of the last component

		#ifdef __LOGGING__
			if (boundaries.size() > 1) {
				N_LOG("bicomps form rooted tree...\n")
			}
		#endif

		for (int cnum = 0; (size_t)(cnum) < boundaries.size(); cnum++) { // step 3 is basically the same as what's in the paper for verifying the ext boundary of a single biconnected component, with some minor adjustments to handle multiple components or trivial components with one edge
			N_LOG("verify bicomp " << cnum << " boundary:\n")
			int root = root_vertices[cnum];

			if (boundaries[cnum].size() == 1) { // check case of trivial bicomp (2 vertices with 1 edge between them, which we need to handle specially since the exterior boundary does not actually form a cycle in this case)
				edge_t e = boundaries[cnum][0];
				V_LOG("trivial bicomp\n")

				if (e.first != root) { // swap edge endpoints if first isn't at root
					int temp = e.first;
					e.first = e.second;
					e.second = temp;
				}
				bool edge_in_g = false;
				for (int u : g.adjLists[e.second]) { // ensure that all the outgoing edges of the single non-root vertex are either in this component or going to an immediate child component whose root vertex is this vertex
													 // if there are any edges which don't do this, then those edges do not have any single biconnected component they belong to and are illegal
					if (u == root) {
						edge_in_g = true;
					} else if (comp_parent[component[u]] != cnum || e.second != root_vertices[component[u]]) {
						L_LOG("====== AUTH FAILED: edge (" << u << ", " << e.second << ") does not belong to any bicomp ======\n\n")
						return false;
					}
				}

				if (!edge_in_g) { // ensure the edge of this trivial bicomp exists in G
					L_LOG("====== AUTH FAILED: edge in trivial bicomp " << cnum << " does not exist in graph ======\n\n")
					return false;
				}

				vert_count += 2; // trivial bicomp has two vertices

				N_LOG("bicomp " << cnum << " boundary good\n")
				continue;
			} // end handling of one-edge bicomp

			for (edge_t e : boundaries[cnum]) { // set up the cycle adjacencies
				if (cycle_adjs[e.first][0] == -1) {
					cycle_adjs[e.first][0] = e.second;
				} else if (cycle_adjs[e.first][1] == -1) {
					cycle_adjs[e.first][1] = e.second;
				} else { // this vertex already has two adjacencies, so there can be no more on the cycle and this edge is illegal
					L_LOG("====== AUTH FAILED: vertex " << e.first << " of bicomp " << cnum << " has three edges in its bicomp incident on it ======\n\n")
					return false;
				}

				if (cycle_adjs[e.second][0] == -1) {
					cycle_adjs[e.second][0] = e.first;
				} else if (cycle_adjs[e.second][1] == -1) {
					cycle_adjs[e.second][1] = e.first;
				} else {
					L_LOG("====== AUTH FAILED: vertex " << e.second << " of bicomp " << cnum << " has three edges in its bicomp incident on it ======\n\n")
					return false;
				}
			}

			std::stack<int> dfs_path; // construct the dfs path according to the cycle of the exterior boundary
									  // every vertex in this dfs has either 0 or 1 tree children, so we do nothing but back up after this and there's no need to keep track of adj list position
			int curr = root;
			int prev = -1;

			do {
				V_LOG("add to dfs path " << curr << "\n")
				dfs_path.emplace(curr);
				dfs_no[curr] = dfs_path.size();

				if (cycle_adjs[curr][1] == -1) { // every vertex in our cycle needs two adjacencies
					L_LOG("====== AUTH FAILED: vertex " << curr << " in the ext boundary of bicomp " << cnum << " lacks two incident edges in the boundary ======\n\n")
					return false;					
				}

				if (cycle_adjs[curr][0] != prev) { // go to next vertex
					prev = curr;
					int temp = cycle_adjs[curr][1]; // swap em around so the parent in the DFS path is always adjacency 0
					cycle_adjs[curr][1] = cycle_adjs[curr][0];
					cycle_adjs[curr][0] = temp;
					curr = cycle_adjs[curr][1];
				} else {
					prev = curr;
					curr = cycle_adjs[curr][1];
				}
			} while (curr != root);

			if (dfs_path.size() != boundaries[cnum].size()) { // we need to have one vertex in the DFS path for every edge in the boundary, or else the boundary is multiple cycles instead of just one
				L_LOG("====== AUTH FAILED: exterior boundary of bicomp " << cnum << " forms more than one cycle ======\n\n")
				return false;	
			}

			int comp_verts = dfs_path.size();
			vert_count += comp_verts;
			N_LOG("boundary forms single cycle with " << dfs_path.size() << " vertices\n")

			int seq_source = dfs_path.top();

			int tree_bottom = dfs_path.top(); 					 // this is used to generate a K4, we need to loop back around from the back edge of the big cycle
			std::vector<bool> seen((size_t)(comp_verts), false); // keep track of whether every edge in the cycle is in G

			while (!dfs_path.empty()) {
				int w = dfs_path.top();
				V_LOG("tree edge w: " << w << " v: " << cycle_adjs[w][0] << "\n")
				V_LOG("seq: " << seq_source << "\n")
				int earliest_outgoing = g.n; // keep track of the lexicographically earliest outgoing back edge of w, s_w, as in the complicated version of the algorithm

				while (!vertex_stacks[w].empty()) { // extend seq
					if (seq_source == vertex_stacks[w].top().first) {  // if the next link in the SP chain matches with the source of seq (there is nothing in between), pop that entry to extend seq
						if (vertex_stacks[w].top().second != -1) seq_source = vertex_stacks[w].top().second; // ensure the stack entry has a tail before we update seq source, this is handled automatically by my sp-tree implementation in the complicated version
						vertex_stacks[w].pop();
					} else { // otherwise, that link and the one at the source of seq form interlacing ears, which is a a 3.4(b) (and subsequently a 3.5(c)) violation
						L_LOG("OOPS, 3.5(c) violation in authentication process (resulting in auth failure) seq src " << seq_source << ", stk end " << vertex_stacks[w].top().first << "\n")		
						// note: since we're relying on a simplified version of the same algorithm to invalidate this certificate as we do to generate certificates to begin with, this sort of defeats the purpose of certifying algorithm (to be sure that an implementation of your algorithm is working by verifying a certificate it spits out using a different algorithm)
						// fortunately, to remedy this we can generate a K4 here and then verify the K4, turning the authentication algorithm itself into a certifying algorithm
						// the generation of the K4 is almost exactly the same as the K4 generation in the SP-OP stack pop violation, but I'm too lazy to extract that code into a function

						negative_cert_K4 k4{};
						k4.a = vertex_stacks[w].top().first;
						k4.b = seq_source;
						k4.c = w;

						for (int a = k4.a; a != k4.b; a = cycle_adjs[a][0]) k4.ab.emplace_back(a, cycle_adjs[a][0]);
						for (int b = k4.b; b != k4.c; b = cycle_adjs[b][0]) k4.bc.emplace_back(b, cycle_adjs[b][0]);

						k4.d = -1;
						int c = k4.c;
						while (k4.d == -1) {
							k4.cd.emplace_back(c, cycle_adjs[c][0]);
							c = cycle_adjs[c][0];

							for (; !vertex_stacks[c].empty(); vertex_stacks[c].pop()) {
								if (vertex_stacks[c].top().first == k4.b) {
									k4.d = c;
									break; 
								}
							}
						}

						for (int d = k4.d; d != root; d = cycle_adjs[d][0]) k4.ad.emplace_back(d, cycle_adjs[d][0]); // the only minor difference here is that we go all the way around the big cycle up to the root of the DFS tree
						k4.ad.emplace_back(root, tree_bottom);									   
						for (int d = tree_bottom; d != k4.a; d = cycle_adjs[d][0]) k4.ad.emplace_back(d, cycle_adjs[d][0]);

						k4.ac.emplace_back(k4.a, k4.c); // both interlacing ears are trivial
						k4.bd.emplace_back(k4.b, k4.d);

						if (k4.authenticate(g)) { // auth our K4 to know for sure our auth of the K23 failed due to the main code for outerplanar detection being wrong
							L_LOG("====== AUTH OF K23 FAILED: ear (" << w << ", " << vertex_stacks[w].top().first << ") interlaces with ear (" << seq_source << ", " << k4.d << ") ======\n\n")
							return false;
						} else { // if we can't even generate a K4, then we don't even know if it was the main code or this authentication code which is wrong
							L_LOG("====== AUTH OF K23 FAILED ======\n\n")
							return false;
						}
					}
				} // end of extending seq

				if (w != root) { // don't loop over children of root, that would not be O(|E|) time but O(|V||E|) (since if there are multiple bicomps rooted at the same vertex, we'd go through the adjacency list of the vertex once per rooted bicomp, leading to very bad performance on, say, K_(1, 100000))
								 // fortunately, we can get away with not doing so; every other vertex in this bicomp will send their back-edges to the root, and when we process the parent bicomp we will go through the adjacencies of the root of this bicomp exactly once for a nice O(|E|) time
					for (int u : g.adjLists[w]) { // loop over children in G
						if (u == cycle_adjs[w][0]) { // if it's on the cycle, don't process it and instead mark the relevant edge in the cycle as seen
							seen[dfs_path.size() - 2] = true;
						} else if (u == cycle_adjs[w][1]) {
							seen[dfs_path.size() - 1] = true;
						} else if (dfs_no[u] < dfs_no[w]) { // if u is a back edge (not the parent, and not the only tree child)
																			  						   // note that this misses the back edge of the big cycle, but that's fine since that back edge will never be part of interlacing ears
							if (u == root || component[u] == cnum) { // if the child is in this biconnected component (either the root or some other vertex in the component), the child edge is internal
								V_LOG("BACK EDGE (" << w << ", " << u << ")\n")
								vertex_stacks[u].emplace(w, -1); // instantly push this entry ending at w onto the relevant vertex stack (we already know it's the losing ear, since the big cycle is always the winning ear)
								
								if (dfs_no[u] < dfs_no[earliest_outgoing]) earliest_outgoing = u; // there is a new earliest outgoing
							} else if (comp_parent[component[u]] != cnum || w != root_vertices[component[u]]) {		   // otherwise, u isn't in this biconnected component
																													   // if the parent of the bicomp it's in isn't this bicomp, then this edge does not belong to any bicomp (since even if w is the root vertex of a child bicomp the child bicomp can't contain u)
																													   // similarly, if w is not the root vertex of the component it's in, the edge doesn't belong (since w would then only be in this bicomp and no other bicomp, and u is not in this bicomp)
								L_LOG("====== AUTH FAILED: edge (" << w << ", " << u << ") does not belong to any bicomp ======\n\n")
								return false;
							}
						}
					}
				} // end of children loop
				
				if (earliest_outgoing != g.n) { // move the seq to the .tail of the relevant stack entry if there's an outgoing back edge
					vertex_stacks[earliest_outgoing].top().second = seq_source;
					seq_source = w;
				}

				dfs_path.pop();
			} // end of interlacing check; at this point, the bicomp has every edge not on its exterior boundary embeddable in the interior of the face formed by it, and has no incident edges (except for possibly at the root) that don't belong to any bicomp

			for (int i = 0; i < comp_verts; i++) { // check every vertex in the cycle is in G
				if (!seen[i]) {
					L_LOG("====== AUTH FAILED: edges in bicomp " << cnum << " cycle do not all belong to graph ======\n\n")
					return false;
				}
			}

			N_LOG("bicomp " << cnum << " boundary good\n")

			dfs_no[root_vertices[cnum]] = 0; // reset the root vertex auxiliary data so another bicomp can use it
			cycle_adjs[root_vertices[cnum]][0] = -1;
			cycle_adjs[root_vertices[cnum]][1] = -1;
		} // end of bicomp verification; at this point, every vertex's adj list has been examined exactly once, except for the root vertex of the final bicomp with no parent


		N_LOG((boundaries.size() > 1 ? "all bicomp boundaries" : "bicomp boundary") << " good... final root left\n")

		for (int u : g.adjLists[root_vertices[boundaries.size() - 1]]) { // check the root of the final bicomp for errant edges
			if ((size_t)(component[u]) != boundaries.size() - 1 && 
			   ((size_t)(comp_parent[component[u]]) != boundaries.size() - 1 || root_vertices[boundaries.size() - 1] != root_vertices[component[u]])) { // errant edge
				L_LOG("====== AUTH FAILED: edge (" << u << ", " << root_vertices[boundaries.size() - 1] << ") does not belong to any bicomp ======\n\n")
				return false;
			}
		}

		N_LOG("final root good...\n")

		vert_count -= (int)(boundaries.size() - 1); // for every bicomp, minus the final one, there is an overlap of a vertex

		if (vert_count != g.n) { // must include every vertex in graph in the bicomps
			L_LOG("====== AUTH FAILED: exterior boundaries span " << vert_count << " vertices, but graph has " << g.n << " vertices ======\n\n")
			return false;	
		}

		N_LOG((boundaries.size() > 1 ? "bicomps span" : "bicomp spans") << " whole graph...\n")
		L_LOG("====== AUTH SUCCESS ======\n\n")

		verified = true;
		return true;
	}
};

struct gsp_sp_op_result {
	bool is_gsp;
	bool is_sp;
	bool is_op;
	std::shared_ptr<certificate> gsp_reason; // it's a shared_ptr because the reasons may point to the same object (e.g. a K4 will be the GSP, SP, and the OP non-reason, and all of these shared_ptrs will point to the same K4)
	std::shared_ptr<certificate> sp_reason;
	std::shared_ptr<certificate> op_reason;

	bool authenticate(graph const& g) {
		//int useless = 1;
		L_LOG("================== AUTHENTICATING GSP-SP-OP RESULT ==================\n") 
		V_LOG(g)
		V_LOG("=============================================================\n")
		L_LOG("\n")

		if (!gsp_reason) {
			L_LOG("ERROR: gsp_reason not given")
			return false;
		}
		if (!gsp_reason->authenticate(g)) return false;

		if (!sp_reason) {
			L_LOG("ERROR: sp_reason not given")
			return false;
		}
		if (!sp_reason->authenticate(g)) return false;

		if (!op_reason) {
			L_LOG("ERROR: op_reason not given")
			return false;
		}
		if (!op_reason->authenticate(g)) return false;

		L_LOG("this graph is " << (is_gsp ? "" : "NOT ") << "GSP, " << (is_sp ? "" : "NOT ") << "SP, and " << (is_op ? "" : "NOT ") << "OP\n")
		return true;
	}
};

#endif