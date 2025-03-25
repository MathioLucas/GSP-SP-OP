// this file contains the meat of the GSP-SP-OP implementation

#ifndef __GSP_SP_OP_HXX__
#define __GSP_SP_OP_HXX__

#include "graph.hxx"
#include "sp-tree.hxx"
#include "logging.hxx"
#include "gsp-sp-op-certificates.hxx"
#include <vector>
#include <stack>
#include <ostream>
#include <memory>
#include <algorithm>

// ====================================================== GSP-SP-OP ===========================================================

std::vector<std::pair<int, int>> get_bicomps(graph const&, std::vector<int>&, gsp_sp_op_result&, int = 0);
void report_K4_non_stack_pop_case(gsp_sp_op_result&, std::vector<int> const&, std::vector<std::stack<sp_chain_stack_entry>>&, int, int, int, int, int, int);
void K23_test(std::shared_ptr<certificate>&, std::vector<int>&, std::vector<int> const&, edge_t, edge_t, int);
int path_contains_edge(std::vector<edge_t> const&, edge_t);

gsp_sp_op_result GSP_SP_OP(graph const& g) {
	gsp_sp_op_result retval{};
	std::shared_ptr<positive_cert_op> op{new positive_cert_op{}}; // the edges making up the exterior boundary of an outerplanar embedding of G

	std::vector<int> cut_verts(g.n, -1);										// cut_verts[i] is -1 if vertex i is not a cut vertex, and a unique number from 0 to (the number of biconnected components of G) - 1 otherwise
	std::vector<edge_t> bicomps = get_bicomps(g, cut_verts, retval);			// get the bicomps of G
	int n_bicomps = (int)(bicomps.size());
	std::vector<sp_tree> cut_vertex_attached_tree((size_t)(n_bicomps));			// cut_vertex_attached_tree[i] is the SP decomposition tree attached to cut vertex i, which we use to stitch together the decomposition trees of bicomps
	std::vector<int> comp(g.n, -1);												// comp[i] is the biconnected component vertex i belongs to (if a vertex belongs to two or more bicomps, then comp[i] is the unique bicomp whose root vertex is not at vertex i)

	std::vector<std::stack<sp_chain_stack_entry>> vertex_stacks((size_t)(g.n)); // vertex_stacks[i] is the per-vertex stack we store for vertex i
	std::vector<int> dfs_no((size_t)(g.n + 1), 0);								// dfs_no[i] is the DFS number of vertex i
	std::vector<int> parent((size_t)(g.n), 0);									// parent[i] is the parent of vertex i

	std::vector<edge_t> ear((size_t)(g.n), edge_t{g.n, g.n});					// ear[i] is the back edge associated with the ear containing the edge between vertex i and its parent in the DFS tree
																				// the lexicographic infinity is defined as (g.n, g.n), and the DFS number of either of its endpoints is also g.n
																				// the first entry in the pair is the source of the back edge of the ear and the second is the sink of the back edge of the ear
	std::vector<sp_tree> seq((size_t)(g.n));									// seq[i] is the last SP subgraph in the i-SP chain (after the DFS backs from vertex i to its parent, it will be finished and ready to use)
	std::vector<int> earliest_outgoing((size_t)(g.n), g.n);						// earliest_outgoing[i] is the source of the lexicographically earliest outgoing ear whose sink is at vertex i, or s_w in the paper
																				// vertex g.n is the ancestral infinity, and all other vertices are proper ancestors of it

	std::vector<char> num_children((size_t)(g.n), 0);							// num_children[i] is the number of children of vertex i in the DFS tree
																				// only needs to be 0, 1, or 2, so it's a 1-byte char to save memory
	std::vector<int> alert((size_t)(g.n), -1);									// alert[i] -1 if there isn't any, or the source of the back edge of that ear if there is a non-trivial ear whose sink is i and whose source is parent[i] (this back edge is 'b' in the paper)
																				// if there are two of such ears with sink on a vertex, there is a 3.5(b) violation and we report a K23 for non-outerplanar-ness (but may still be SP)
																				// NOTE: it might be possible to recover b, alert, and the number of children just by examining the vertex stack of v (seeing if the top entry on that stack has an end at w, which is O(1) time)
																				// this would save a fair bit of memory (you wouldn't need either of the above two arrays), but I haven't proven the correctness of it so I'll just stick with the approach given in the paper

																				// NOTE: technically many of the above could all go in the DFS stack frames but my implementation keeps an array of them for every vertex to keep the stack frames simple
																				// the space usage is still O(|V| + |E|), but maybe if I were to optimize the memory use in the future I'd get rid of most of these and put them in the DFS stack frames 
																				// this would probably around halve the memory use for sparser graphs where the DFS frequently backtracks, and massively cut it down for early exits due to K4

	std::stack<std::pair<int, int>> dfs;  // I implement the DFS iteratively rather than recursively to avoid stack overflows on exceptionally large input graphs
										  // the first entry in each pair of this stack represents the w of one recursive call, and the second entry the current index in the adjacency list we're looking at for that vertex
										  // my predecessor has done something similar in their "dfs2-using-2d-vector.cpp", but they scan over the entire adjacency list every time to find a not-yet-visited descendant instead of keeping track of the adj address (which is O(|E^2|) time instead of O(|E|) like my approach is)

	dfs_no[g.n] = g.n; // dfs_no[g.n] is g.n, a special value of the DFS number for the lexicographic and ancestral infinity

	bool do_k23_edge_replacement = true; // flag to avoid uselessly redoing the work of checking if we've already successfully replaced a fake edge on a K23 to ensure we only make that check once

	for (int bicomp = 0; bicomp < n_bicomps; bicomp++) {
		N_LOG("BICOMP " << bicomp << "\n")

		if (!retval.op_reason) {
			op->boundaries.emplace_back(); // add a new bicomp to the exterior boundary if we haven't found K23 yet
		}
		std::vector<edge_t>& ext_boundary = op->boundaries.back();

		int root = bicomps[bicomp].first;
		int next;
		if (!retval.sp_reason && bicomp > 0 && bicomp < n_bicomps - 1) { // case where 1 < i < h and sp != false in the paper, we generate a potentially fake edge between the two cut vertices attached to this bicomp
			next = bicomps[bicomp - 1].first; // go to the cut vertex of the previous bicomp
		} else {
			next = bicomps[bicomp].second; // proceed as normal
		}

		dfs.emplace(root, -1); // we will never actually iterate over the root; we return from the bicomp immediately before backing into the root, since we're done then
		dfs.emplace(next, 0); // force the first edge of the DFS
		bool fake_edge = true; // keep track of whether the edge we added was fake or not; we'll check this at the end of the processing of this bicomp

		dfs_no[root] = 1;
		parent[root] = -1;
		dfs_no[next] = 2;
		parent[next] = root;
		comp[next] = bicomp;
		int curr_dfs = 3;

		// ==================== SP-OP begins here ====================

		while (!dfs.empty()) { // the DFS is implemented iteratively rather than recursively to avoid stack overflows on exceptionally large input graphs, but is equivalent to what's in the paper
			std::pair<int, int> p = dfs.top();
			int v = parent[p.first]; // v is the parent of w in the DFS tree, u is the current vertex being examined in w's adjacency list (which I will call the "child" but this is not quite accurate since the edge between u and w could be a back edge)
			int w = p.first;
			int u = g.adjLists[p.first][p.second];

			if (comp[u] == -1 || comp[u] == bicomp) { // skip over child if it isn't part of this bicomp
				V_LOG("v: " << v << " w: " << w << " u: " << u << "\n")
				V_LOG("seq_w: " << seq[w] << ", seq_u: " << seq[u] << "\n")
				if (dfs_no[u] == 0) { // the first time we see this child node, if it's unvisited and it's in this bicomp, make a recursive call
					dfs.push(std::pair{u, 0});
					parent[u] = w;
					dfs_no[u] = curr_dfs++;
					comp[u] = bicomp;
					num_children[w]++; // will never exceed 3 if the graph is outerplanar (we are guaranteed to find a K23 violation if it becomes 3)
									   // if it's not outerplanar, the char might suffer from integer overflow, but this is ok since we don't care about the contents this array once G is not outerplanar

					continue;
				}
				// otherwise, we've seen this node before
				bool child_back_edge = (dfs_no[u] < dfs_no[w] && u != v); // keep track of whether the edge between w and u is an incoming back-edge or tree-edge 
																		  // note that my implementation does not explicitly store the associated ear of a back-edge, since ear(u ↶ w) is trivially u ↶ w
				#ifdef __LOGGING__
					if (child_back_edge) N_LOG("BACK EDGE (" << w << ", " << u << ")\n")
				#endif	
																  
				if (parent[u] == w) { // w->u is a tree edge, meaning last time we saw this node we made a recursive call and are now returning from it after processing the DFS subtree rooted at u
					N_LOG("tree edge (" << w << ", " << u << ")\n")
					// --- update-seq in the paper begins here ---
					for (; !vertex_stacks[w].empty(); vertex_stacks[w].pop()) { // combine all of the links in the u-SP chain whose SPs have a source at w into one (extending its seq)
						if (seq[u].source() != vertex_stacks[w].top().end) { // if the next link in the chain isn't the one on the top of the stack, then the source of seq won't be equal to the end of the top entry
																			 // this results in a theorem 3.4(b) violation and a K4 is reported here
							N_LOG("OOPS, 3.4b due to POPPING STACK child seq " << seq[u] << " parent seq " << seq[w] << "\n")
							std::shared_ptr<negative_cert_K4> k4{new negative_cert_K4{}};

							k4->b = seq[u].source(); 					 // the sinks of the two interlacing ears are the source of our seq (the first link in the w-SP chain) and the .end of the stack entry we're popping
							k4->a = vertex_stacks[w].top().end;			 // these are part of the K4 subdivision, so we include them in our certificate
							k4->c = w;									 // the source of the interlacing ear we're popping off the stack is, of course, w
							edge_t holding_ear = ear[u]; 				 // this is the ear the two interlacing ears are s-attached to, forming a 3.4(b) violation
																		 // the paper suggests using ear[vertex_stacks[w].top().end] instead, but these two are the same (if they weren't we'd have a 3.4(a) violation rather than a 3.4(b) one, which would've been detected when we backed off of vertex_stacks[w].top().end)
																		 // I chose ear[u] because if I ever optimize away the per-vertex ear array into the DFS stack frames we already have ear[u] on hand, whereas we'd need to compute ear[vertex_stacks[w].top().end] with another DFS

							for (int a = k4->a; a != k4->b; a = parent[a]) k4->ab.emplace_back(a, parent[a]); // trace the ear up; I've ordered the terminating vertices of the K4 subdivision a, b, c, d from ancestrally latest to earliest
							for (int b = k4->b; b != k4->c; b = parent[b]) k4->bc.emplace_back(b, parent[b]);

							k4->d = -1;
							int c = k4->c; // cd is a bit more complicated; we don't know where d, the source of the ear whose sink is the source of seq[u], is at
										   // we do know, though, that the ear is stored on some vertex stack (since we've already processed it when we backed up from its sink), and this stack's vertex must be d
										   // we examine every stack as we walk up the tree until we find an entry whose ear ends at k4->b, and the vertex of that stack must be a valid d
										   // we need to check the entire stack instead of just the top since there may be another violating ear we haven't seen yet that messes up the stack ordering
							while (k4->d == -1) {
								k4->cd.emplace_back(c, parent[c]);
								c = parent[c];

								for (; !vertex_stacks[c].empty(); vertex_stacks[c].pop()) {
									if (vertex_stacks[c].top().end == k4->b) {
										k4->d = c;
										break; 	// break to preserve the stack entry so we can generate bd with its SP
									}
								}
							}

							for (int d = k4->d; d != holding_ear.second; d = parent[d]) k4->ad.emplace_back(d, parent[d]); // keep walking up the tree until we reach the source of the ear that holds the two interlacing ears
							k4->ad.emplace_back(holding_ear.second, holding_ear.first);									   // add the back edge of that ear (in reverse direction since the sink of the back edge is the source of the ear), forming a cycle
							for (int d = holding_ear.first; d != k4->a; d = parent[d]) k4->ad.emplace_back(d, parent[d]);  // walk up the remainder of the holding ear until we loop back around to a, completing the cycle

							int ear1 = vertex_stacks[k4->d].top().SP.underlying_tree_path_source(); // add the interlacing ear between b and d
							k4->bd.emplace_back(k4->d, ear1);
							for (; ear1 != k4->b; ear1 = parent[ear1]) k4->bd.emplace_back(ear1, parent[ear1]);
							int ear2 = vertex_stacks[k4->c].top().SP.underlying_tree_path_source(); // add the interlacing ear between a and c and we are done
							k4->ac.emplace_back(k4->c, ear2);
							for (; ear2 != k4->a; ear2 = parent[ear2]) k4->ac.emplace_back(ear2, parent[ear2]);

							retval.gsp_reason = k4;
							break;
						} // end of K4 reporting

						seq[u].compose(std::move(vertex_stacks[w].top().SP), c_type::antiparallel); // the antiparallel composition is used here to 'mark' the SP_(x, y) being popped off the stack as flipped, but is otherwise identical to a parallel composition; this will come in handy for the authentication algorithm
						seq[u].l_compose(std::move(vertex_stacks[w].top().tail), c_type::series);
					}
					// ---- update-seq in the paper ends here ----

					if (retval.gsp_reason) break; // stop if we found a K4
				}

				if (parent[u] == w || child_back_edge) { // bit of an ugly check here but it's better than copying update-ear-of-parent twice
					// ---- update-ear-of-parent in the paper begins here ----
					edge_t ear_f = (child_back_edge ? edge_t{w, u} : ear[u]); // if the edge to the child is a back edge then its ear is itself; otherwise if it's a tree edge it's the previously-computed value of ear[u]
					sp_tree seq_u = (child_back_edge ? sp_tree{u, w} : std::move(seq[u]));

					if (dfs_no[ear_f.second] < dfs_no[ear[w].second]) { // case (b). note that by the nature of DFS and the definition of ear, it's impossible for neither t(ear(v->w)) and t(ear(w->u)) to be ancestors of the other, so it suffices to compare the DFS numbers to determine ancestorship
						if (ear[w].first != g.n) { // check if ear[w] isn't the lexicographic infinity (i.e we've already processed at least one child or outgoing back edge of w)
							if (!retval.op_reason && ear[w].first != w) K23_test(retval.op_reason, alert, parent, ear[w], ear_f, w); // if we haven't found a K23 yet and the ear we completed is non-trivial (i.e its back edge is outgoing from w), we need to test for a 3.5(a) or (b) violation
							if (seq[w].source() != ear[w].second) { // in the above case, the child we're looking at "wins" and its ear continues onto the parent, cutting off the ear of the previous winner
																   // this means that the ear of the previous winner is done and must be completely converted to an SP subgraph (i.e. its seq must be the whole ear with no other links in its SP-chain)
																   // if this isn't true, then there is a theorem 3.4(a) violation and we need to report a K4, which is done in this branch
								N_LOG("OOPS, 3.4a due to CASE B prev winner " << seq[w] << " prev winner ear (" << ear[w].first << ", " << ear[w].second << ")\n")
								report_K4_non_stack_pop_case(retval, parent, vertex_stacks, seq[w].source(), w, ear[w].second, ear[w].first, ear_f.second, ear_f.first); // generate k4
								break;
							} // otherwise the completed ear is good, we can store it on the relevant stack

							N_LOG("CASE B (ear exists): placed " << seq[w] << " onto stk " << ear[w].second << "\n")
							vertex_stacks[ear[w].second].emplace(std::move(seq[w]), w, sp_tree{});
							earliest_outgoing[w] = ear[w].second; // the ear we just finished is guaranteed to be lexicographically earlier than all other previous finished ears, so we can update s_w right away
						}
						ear[w] = ear_f;
						seq[w] = std::move(seq_u); // the partial sequence of the winning ear becomes the sequence of the child ear
						N_LOG("CASE B (replace seq): current winning seq " << seq[w] << "\n")
					} else { // case (a), or case (c)
						if (seq_u.source() != ear_f.second) { // same as before, if the ear that got cut off isn't completely processed (in case (a)) or if the sources and sinks of the seqs don't match up (in case (c)) there's a violation and we report a K4
															  // this time, it may be a 3.4(a) violation (if we're in case (a), or if we're in case (c)) or it may be a 3.4(b) violation (if we're in case (c))
							N_LOG("OOPS, 3.4a/b due to CASE A/C child seq " << seq_u << " child ear (" << ear_f.first << ", " << ear_f.second << ")\n")
							report_K4_non_stack_pop_case(retval, parent, vertex_stacks, seq_u.source(), w, ear_f.second, ear_f.first, ear[w].second, ear[w].first); // generate and return k4
							break;
						}

						if (dfs_no[ear_f.second] == dfs_no[ear[w].second]) { // case (c). in this case, the sources of the child ear and the currently winning ear are the same, so which wins depends entirely on parts (ii) and (iii) of the lexicographic ordering
																			 // regardless of which ear wins, neither of the SP chains of the ears can have any links in them; if the sequence of the losing ear is not the entire SP-chain for its sink vertex, there's a 3.4(a) violation, and if the sequence of the winning ear isn't, there's a 3.4 (b) violation
							if (!retval.op_reason && !child_back_edge && ear[w].first != w) K23_test(retval.op_reason, alert, parent, ear_f, ear[w], w); // if both ears are non-trivial, test for a K23

							if (seq[w].source() != ear[w].second) { // we already checked the ear of the child, now we just need to check the currently winning ear
																	// if it's not the whole chain, there's a K4 to be reported
								N_LOG("OOPS, 3.4a/b due to CASE C parent seq " << seq[w] << " parent ear (" << ear[w].first << ", " << ear[w].second << ")\n")
								report_K4_non_stack_pop_case(retval, parent, vertex_stacks, seq[w].source(), w, ear[w].second, ear[w].first, ear_f.second, ear_f.first); // generate and return k4
								break;
							}
							seq[w].compose(std::move(seq_u), c_type::parallel); // both sequences share a source and sink if they're both complete, so whichever ear wins, that ear's sequence must be merged with the one of the losing ear
							N_LOG("CASE C: current winning seq after merge " << seq[w] << "\n")

							if ((ear[w].first == w || dfs_no[ear_f.first] < dfs_no[ear[w].first]) && ear_f.first != w) { // if ear[w].first = w, then s(ear(v->w)) is a proper ancestor of s(ear(w->u)) and thus is lexicographically later (part (ii) of lex-ordering); if ear_f.first = w, then it's the same but the other way around; if neither are the case, then neither is an ancestor of the other and we compare DFS numbers of the sources (part (iii))
								ear[w] = ear_f;
							}
						} else { // case (a). similar to case (b), but the child we're looking at does not cut off the ear of the current winner and instead it is the one that gets cut off
								 // we already checked that the child ear is completely processed, so no need to do that
							if (!retval.op_reason && !child_back_edge) K23_test(retval.op_reason, alert, parent, ear_f, ear[w], w); // if the conpleted ear (child ear) is non-trivial, test for K23

							if (!vertex_stacks[ear_f.second].empty() && vertex_stacks[ear_f.second].top().end == w) { // if there's already an entry in the relevant stack whose SP ends here, we merge the completed ear with a parallel composition instead of adding a new entry
																													  // we didn't need to do this in case (b) because in that case the vertex whose stack we were adding to was a proper ancestor of any other vertex whose stack we've added an entry ending at w to before
								N_LOG("CASE A (merge onto existing stack entry for stk " << ear_f.second << "): current child seq before merge " << seq_u << "\n")
								vertex_stacks[ear_f.second].top().SP.compose(std::move(seq_u), c_type::parallel);
							} else {
								N_LOG("CASE A (new stack entry): placed " << seq_u << " onto stk " << ear_f.second << " (earliest outgoing " << earliest_outgoing[w] << ")\n")
								vertex_stacks[ear_f.second].emplace(std::move(seq_u), w, sp_tree{});
								if (dfs_no[ear_f.second] < dfs_no[earliest_outgoing[w]]) { // unlike in case (b), the ear we just finished is not guaranteed to be lexicographically earlier than the previous finished ears, so we need to lexi-compare them
																							// as before, the source of one ear is guaranteed to be an ancestor of the source of the other ear, so it suffices to compare DFS numbers
									earliest_outgoing[w] = ear_f.second;
								}
							}
						}
					}
					// ----- update-ear-of-parent in the paper ends here -----
				}
			}

			if ((size_t)(++dfs.top().second) >= g.adjLists[p.first].size()) { // if we're done processing this vertex's adj list, we're ready to return to the parent
				if (w != root) {
					if (earliest_outgoing[w] != g.n) { // if there's an ear whose sink is at w, we've created a new link in the SP chain for all proper ancestors of this vertex until we reach the source of that ear
													   // we need to move our current sequence to the tail of that new link to reflect this
						N_LOG("EARLIEST OUTGOING " << earliest_outgoing[w] << ": moved current winning seq " << seq[w] << " to vertex stack entry tail with SP " << vertex_stacks[earliest_outgoing[w]].top().SP << "\n")
						vertex_stacks[earliest_outgoing[w]].top().tail = std::move(seq[w]);
					}

					if (!retval.op_reason) { // if we haven't found a K23, we extend the exterior boundary when we back up to parent
						if (num_children[w] == 0) { // case 1; if the vertex has no children, we add the parent edge and the lexicographically earliest outgoing back-edge (the ear of the parent edge) to the boundary
							ext_boundary.emplace_back(v, w); // add parent edge
							ext_boundary.emplace_back(ear[w].first, ear[w].second); // add back edge
						} else if (num_children[w] == 1) { // case 2; the vertex has one child
							if (ear[w].first == w) { // case 2(a); the lexicographically earliest outgoing back-edge of w is lexicographically earlier than the back-edge of the ear of the child of w (this is true if the back-edge is the winning ear and trivial)
													 // in this case, the parent edge is internal, since we can embed it inside the cycle formed by this back edge, ear[u], and the tree paths between them, where u is the sole child of w
													 // we add the back edge to the boundary, not the parent edge
								ext_boundary.emplace_back(ear[w].first, ear[w].second);
							} else { // case 2(b); the lexicographically earliest outgoing back edge of w either does not exist or is not lexicographically earlier than the back-edge corresponding to the ear of the child
									 // in this case, the back edge is internal, since we can embed it inside the cycle formed by the ear of the child and the tree-path from the sink of that ear up to its source
									 // we add the parent edge to the boundary since it is part of this cycle, not the back edge
								ext_boundary.emplace_back(v, w);
							}
						} // if the vertex has two children, we don't add any edges to the exterior boundary, and if it has more than two there's a K23 which we'd've already found
					}

					for (int u1 : g.adjLists[next]) { // iterate over the adj list of the tree child of the root an additional time to check if the first edge in the DFS tree doesn't exist (as it might when the bicomps form a chain)
													  // over the whole algorithm this takes O(|E|) time, we are guaranteed to run this loop at most once per vertex (since we always run it on a non-root-vertex of a bicomp, and when bicomps overlap exactly one of those bicomps do not have a root vertex on the overlap)
						if (u1 == root) {
							fake_edge = false;
							break;
						}
					}

					if (v == root) { // there will never be more than one tree edge going out of the root, otherwise the root would be a cut vertex (we assume G is biconnected in SP-OP)
									 // thus this happens exactly once
						seq[w].compose((fake_edge ? sp_tree{} : sp_tree{v, w}), c_type::parallel); // do the last composition  (if the first tree edge is fake don't include it)					
																								   // it's different from all the others since the root ear is the only cyclic ear, so we need a parallel composition (source to source) rather than a series one (sink to source) to reduce it to an edge

						if (cut_verts[w] != -1) { // add the tree hanging off this cut vertex if it exists, it needs to be composed in series with seq[w] since this is the final edge
							seq[w].compose(std::move(cut_vertex_attached_tree[cut_verts[w]]), c_type::series);
						}
						break; // at this point we have assembled the entire exterior boundary and decomposition tree assuming no K23 or K4, so there's no point doing anything more (we will never detect a K23 or K4 at the root)

					} else { // otherwise, it's just a regular old vertex
						if (cut_verts[w] != -1) { // add the tree hanging off this cut vertex if it exists, for all edges except for the final edge this tree will be dangling along with the parent edge
							cut_vertex_attached_tree[cut_verts[w]].l_compose(sp_tree{w, v}, c_type::dangling);
							seq[w].compose(std::move(cut_vertex_attached_tree[cut_verts[w]]), c_type::series);
						} else {
							seq[w].compose(sp_tree{w, v}, c_type::series); // extend the sequence to include the parent edge, completing the w-SP chain
						}		
					}
				}

				dfs.pop(); // end the recursive call
			}
		}

		// ==================== SP-OP ends here ====================

		dfs_no[root] = 0; // clear the dfs_no of the root so the DFS for future bicomps can revisit it

		if (!retval.op_reason) {
			N_LOG("no K23 found\n")

			if (ext_boundary.size() == 2 && ext_boundary.back() == edge_t{g.n, g.n}) { // if the bicomp is trivial (2 vertices and 1 edge), we need to handle it with special care; in SP-OP the paper assumes our graph is biconnected and thus has at least two edges, and in particular in the case of a trivial bicomp the single edge will not have a defined ear
																					   // this messes up the exterior boundary, since the single edge has 0 children and thus we add both the single edge and its nonexistent ear to the boundary
																					   // fortunately, it's an easy fix; we detect the exterior boundary being messed up in this particular way and just eliminate the nonexistent ear edge from it to correct it
				ext_boundary.pop_back(); // delete nonexistent edge from exterior boundary if the bicomp is trivial
			}

			ext_boundary.shrink_to_fit(); // if we haven't found a K23, shrink the ext boundary of this bicomp to save memory
		}

		if (fake_edge) { // handle fake edge
			edge_t fake = edge_t{root, next};

			if (retval.gsp_reason) { // if there is a fake edge and a reported K4 contains it, it's not the end of the world; we remove the fake edge and put the resulting T4 as an sp reason, and we can try again to see if the graph is still GSP
				std::shared_ptr<negative_cert_K4> k4 = std::dynamic_pointer_cast<negative_cert_K4>(retval.gsp_reason); // this is a horrible mess of spaghetti; I really should have put the paths in an array earlier and probably designed the certificates so that I don't need to do these annoying casts but oh well
				std::vector<edge_t> * k4_paths[6] = {&k4->ab, &k4->ac, &k4->ad, &k4->bc, &k4->bd, &k4->cd};
				int k4_verts[4] = {k4->a, k4->b, k4->c, k4->d};
				static const int k4_t4_translation[6][5] = {{1, 3, 2, 4, 5}, {0, 3, 2, 5, 4}, {0, 4, 1, 5, 3}, {0, 1, 4, 5, 2}, {0, 2, 3, 5, 1}, {1, 2, 3, 4, 0}}; // c1a, c2a, c1b, c2b, ab respectively
				static const int k4_t4_endpoint_translation[6][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}, {1, 2, 0, 3}, {1, 3, 0, 2}, {2, 3, 0, 1}}; // c1, c2, a, b respectively

				int pnum = 0;
				for (; pnum < 6; pnum++) {
					if (path_contains_edge(*(k4_paths[pnum]), fake) != -1) break; // if this path contains the fake edge, eliminate it and put the other paths in a T4
																				  // note that by a similar argument to the one about replacing edges in K23's in the paper, either a) the fake edge is between two terminating vertices of the K4 or b) there is a K4 in the bicomp without the fake edge
				}

				if (pnum != 6) {
					N_LOG("FAKE EDGE IN K4 (pnum " << pnum << "), GENERATE T4\n")
					std::shared_ptr<negative_cert_T4> t4{new negative_cert_T4{}}; // generate t4

					t4->c1a = std::move(*(k4_paths[k4_t4_translation[pnum][0]]));
					t4->c2a = std::move(*(k4_paths[k4_t4_translation[pnum][1]]));
					t4->c1b = std::move(*(k4_paths[k4_t4_translation[pnum][2]]));
					t4->c2b = std::move(*(k4_paths[k4_t4_translation[pnum][3]]));
					t4->ab = std::move(*(k4_paths[k4_t4_translation[pnum][4]]));
					t4->c1 = k4_verts[k4_t4_endpoint_translation[pnum][0]];
					t4->c2 = k4_verts[k4_t4_endpoint_translation[pnum][1]];
					t4->a = k4_verts[k4_t4_endpoint_translation[pnum][2]];
					t4->b = k4_verts[k4_t4_endpoint_translation[pnum][3]];

					retval.sp_reason = t4;
					retval.gsp_reason.reset(); // remove k4

					for (int i = 0; i < g.n; i++) { // reset bicomp data (this requires checking over the whole graph because of how I represented the bicomps, but fortunately it is guaranteed to happen at most once so this step is O(|V|) time)
						if (comp[i] == bicomp) {
							dfs_no[i] = 0;
							parent[i] = 0;
							ear[i] = edge_t{g.n, g.n};
							num_children[i] = 0;
							alert[i] = -1;
							earliest_outgoing[i] = g.n;
							seq[i] = sp_tree{};
							vertex_stacks[i] = std::stack<sp_chain_stack_entry>{};
						}
					}
					op->boundaries.pop_back(); // reset ext boundary

					bicomp--; // reprocess this bicomp
				}
			} // end of fake edge K4 handling

			if (retval.op_reason && do_k23_edge_replacement) { // if there is a fake edge and a reported K23 contains it, we need to correct the K23 by replacing that edge with an ear corresponding to a tree child of this vertex not in the K23
				std::shared_ptr<negative_cert_K23> k23 = std::dynamic_pointer_cast<negative_cert_K23>(retval.op_reason);
				std::vector<edge_t> * k23_paths[3] = {&k23->one, &k23->two, &k23->three}; 

				int pnum = 0;
				int path_ind;
				for (; pnum < 3; pnum++) {
					path_ind = path_contains_edge(*(k23_paths[pnum]), fake);
					if (path_ind != -1) break;
				}

				if (pnum != 3) {
					std::vector<edge_t>& violating_path = *(k23_paths[pnum]);
					N_LOG("FAKE EDGE IN K23 (" << violating_path[path_ind].first << ", " << violating_path[path_ind].second << "), REPLACE WITH PATH\n")

					std::vector<edge_t> splice_path;
					bool in_k23[g.n];
					for (int i = 0; i < g.n; i++) in_k23[i] = false;

					for (std::vector<edge_t> * path : k23_paths) { // mark every vertex in the K23
						for (edge_t e : *path) {
							in_k23[e.first] = true;
							in_k23[e.second] = true;
							V_LOG("(" << e.first << ", " << e.second << ") in K23\n")
						}
					}

					for (int u2 : g.adjLists[next]) { // find tree child not in K23
						if (comp[u2] == bicomp && parent[u2] == next && !in_k23[u2]) { // found it, generate tree path
							V_LOG("FOUND TREE CHILD OF NEXT " << next << " NOT IN K23: " << u2 << ", ear (" << ear[u2].first << ", " << ear[u2].second << ")\n")
							splice_path.emplace_back(ear[u2].first, root);
							for (int i = ear[u2].first; i != next; i = parent[i]) splice_path.emplace_back(parent[i], i);
							break;
						}
					}

					std::reverse(splice_path.begin(), splice_path.end()); // reverse generated path; the violating edge goes from next to root, but the path goes from root to next, so we need to reverse it to bring the endpoints in line
																		  // I can't just construct the path in reverse; that would not be O(|V|) time since inserting at the front of a vector is O(n) time, while inserting at the end is O(1)
					violating_path.erase(violating_path.begin() + path_ind); // remove fake edge
					violating_path.insert(violating_path.begin() + path_ind, splice_path.begin(), splice_path.end()); // splice in new path (this is O(|V|) time according to cppreference.com)
				}
			} // end of fake edge K23 handling
		} // end of fake edge handling

		if (retval.op_reason) do_k23_edge_replacement = false; // if we found a K23 this bicomp, there's no need to check for K23 fake edge replacement anymore (either it didn't have a fake edge or it did have a fake edge and that edge got replaced above)

		if (retval.gsp_reason) { // if K4, we are done
			retval.sp_reason = retval.gsp_reason;
			retval.op_reason = retval.gsp_reason;
			break;
		}

		if (cut_verts[root] != -1) { // if there's a cut vertex here, attach the relevant tree to our seq (if there is no tree attached to this cut vertex yet, then nothing will happen)
			#ifdef __VERBOSE_LOGGING__
			if (cut_vertex_attached_tree[cut_verts[root]].root) {
				V_LOG("combine tree " << cut_vertex_attached_tree[cut_verts[root]] << " with " << seq[next] << " (bicomp " << bicomp << ")\n");
			}
			#endif

			seq[next].compose(std::move(cut_vertex_attached_tree[cut_verts[root]]), c_type::dangling);
		}

		if (bicomp < n_bicomps - 1) { // if this isn't the root bicomp, attach our finished tree to the relevant cut vertex
			V_LOG("ATTACH " << seq[next] << " to cut vertex " << root << " (bicomp " << bicomp << ")\n");
			cut_vertex_attached_tree[cut_verts[root]] = std::move(seq[next]);
		} else { // otherwise, this is the last bicomp
			if (!retval.gsp_reason) { // if there are no K4s preventing GSP, output the completed certificate
				std::shared_ptr<positive_cert_gsp> gsp{new positive_cert_gsp{}}; // move the finished decomp tree into place if we have not found k4

				gsp->decomposition = std::move(seq[next]);
				//gsp->decomposition.deantiparallelize(); // optional step for deantiparallelization (it's still O(|V| + |E|) time to do this, but the paper doesn't suggest it)
														  // UPDATE: deantiparallelization no longer works on non-biconnected graphs for reasons explained in gsp-sp-op-certificates.hxx (in short, the dangling composition is attached to the wrong edge if there's a swap)
				retval.gsp_reason = gsp;
				retval.is_gsp = true;
				N_LOG("graph is GSP\n")

				if (!retval.sp_reason) { // if we haven't generated a negative cert for SP yet, also move the positive cert in as the reason for SP
					retval.sp_reason = gsp;
					gsp->is_sp = true;
					retval.is_sp = true;
					N_LOG("graph is SP\n")
				}
			}
		}
	} // processing of all bicomps ends here

	if (!retval.op_reason) { // if we haven't found a K23, move in the exterior boundary
		retval.op_reason = op;
		retval.is_op = true;
		N_LOG("graph is OP\n")
	}

	#ifdef __VERBOSE_LOGGING__
		for (int i = 0; i < g.n; i++) {
			V_LOG("vertex " << i << " ear: (" << ear[i].first << ", " << ear[i].second << ")\n")
			V_LOG("vertex " << i << " parent: " << parent[i] << "\n")
			V_LOG("vertex " << i << " dfs_no: " << dfs_no[i] << "\n")
		}
	#endif

	return retval;
}

std::vector<edge_t> get_bicomps(graph const& g, std::vector<int>& cut_verts, gsp_sp_op_result& cert_out, int root) {	 // determine the biconnected components of G (rather, their root vertices (the first element of the pairs) and an outgoing edge from those root vertices (whose other vertex is the second element of those pairs))
																														 // the paper suggests using "Tarjan's algorithm" to accomplish this
																														 // I do not know what that is (and cannot figure out how to access the reference in the paper), so instead I have implemented the algorithm in the COMP-4540 textbook to do the same thing, which seems to have the required property that the bicomps are produced in the order of a rooted tree													   							   				 	 			 // this bit of the algorithm also generates negative certificates for a bicomp containing three or more cut vertices or a cut vertex contained in three or more bicomps																													 // the paper does not elaborate on how to generate the tri-cut-comp and tri-comp-cut certificates, so I've filled in the details there a bit; I haven't written up a correctness or time complexity proof for my algorithm but if you want I can
	std::vector<int> dfs_no((size_t)(g.n), 0); // dfs_no[i] is the DFS number of vertex i
	std::vector<int> parent((size_t)(g.n), 0); // parent[i] is the parent of vertex i in its adjacency list
	std::vector<int> low((size_t)(g.n), 0);    // low[i] is the DFS number of the sink of the back edge whose sink has a lowest DFS number among all back edges outgoing from a descendant of vertex i, or just dfs_no[i] if no such back edge exists

	std::vector<edge_t> retval; // this list stores the outgoing edges from the root vertices of the bicomps
								// there's no need to construct whole subgraphs for every bicomp when we'll just end up DFSing from the root vertices of the bicomps anyway in GSP_SP_OP
	std::stack<std::pair<int, int>> dfs;

	dfs.emplace(root, 0);
	dfs_no[root] = 1;
	low[root] = 1;
	parent[root] = -1;
	int curr_dfs = 2;
	bool root_cut = false; // root_cut is true if we have seen two or more root vertices of bicomps at the root of the DFS; for all other vertices, we can have at most two bicomps rooted at a vertex before we need to generate a tri-comp-cut negative certificate, but for the root of the DFS we can have at most three (since the root bicomp has no parent)

	while (!dfs.empty()) {
		std::pair<int, int> p = dfs.top();
		int w = p.first;
		int u = g.adjLists[p.first][p.second];
		if (dfs_no[u] == 0) { // recurse if we haven't seen u
			dfs.push(std::pair{u, 0});
			parent[u] = w;
			dfs_no[u] = curr_dfs++;
			low[u] = dfs_no[u]; // initialize low[u]; this will never actually cause a bicomp to be detected (since dfs_no[w] < dfs_no[u] by the nature of DFS and dfs_no[u] = low[u] by this statement implies low[u] >/= dfs_no[w] if there are no outgoing back edges of a descendant of u) 
								// we do it anyway, though, so that low[u] is always defined
			continue;
		} // otherwise, we have visited u before

		if (parent[u] == w) { // w->u is a tree edge
			if (low[u] >= dfs_no[w]) { // detect a biconnected component; if the lowest back edge coming out of a descendant of the child isn't a proper ancestor of the parent then no back edge coming out of a descendant of the child is and w must be a cut vertex (u is part of a biconnected component incident on w)
				if (cut_verts[w] != -1) { // we already have a cut vertex here
					if (w != root || root_cut) { // if it isn't at the root (or we have found three bicomps rooted at the root) generate a tri-comp-cut
						if (!cert_out.sp_reason) {
							N_LOG("NON-SP, three component cut vertex at " << w << "\n")
							std::shared_ptr<negative_cert_tri_comp_cut> cut{new negative_cert_tri_comp_cut{}};
							cut->v = w;
							cert_out.sp_reason = cut;
						}
					} else { // if it is at the root, flag root-cut so a third bicomp will generate the negative certificate
						root_cut = true;
					}
				} else { // this is the first cut vertex here, give it a unique number
					cut_verts[w] = retval.size();
				}

				retval.emplace_back(w, u); // store the root vertex of this bicomp (w) and an adjacent vertex in the same component (u), so we can start off the DFS using this information later in GSP-SP-OP
										   // there is no need to generate the whole bicomp
			}

			if (low[u] < low[w]) low[w] = low[u]; // update low, if this child has a lower value then it follows that w also has at most that value, since the descendants of a child are also descendants of its parent
												  // after processing all the children of w, we will have correctly computed low[w] (by induction on the height of w in the DFS tree)
		} else if (dfs_no[u] < dfs_no[w] && u != parent[w]) { // u ↶ w is a back edge
			if (dfs_no[u] < low[w]) low[w] = dfs_no[u]; // update low, if the outgoing back edge has a lower DFS number at its sink then it follows that w has at most that value as its low, since w is trivially a descendant of w
		}

		if ((size_t)(++dfs.top().second) >= g.adjLists[p.first].size()) { // we are done processing the adjacency list of w and are backing up to the parent
			dfs.pop(); // end the recursive call
		}
	}

	int n_bicomps = (int)(retval.size()); // cast it to an int to avoid annoying compiler warnings
	N_LOG(n_bicomps << " bicomp" << (n_bicomps == 1 ? "" : "s") << " found\n")
	for (int i = 0; i < n_bicomps; i++) {
		V_LOG("bicomp " << i << ": root " << retval[i].first << ", edge " << retval[i].second << "\n")
	}

	if (!root_cut) cut_verts[root] = -1; // if the root is not a cut vertex, don't mark it as one (the bicomp detection algorithm marks the root of every bicomp found as a cut-vertex, even the root bicomp)

	retval.shrink_to_fit();
	if (cert_out.sp_reason) return retval; // if we found tri-comp-cut, return early (no need to check for tri-cut-comp or reorder the bicomps)

	N_LOG("no tri-comp-cut found\n")

	// now that we've found the bicomps, we will check if they form a chain (generating a tri-cut-comp if they don't)
	// we do this by walking up the DFS tree from the root vertex of each bicomp in the list until we reach root, marking the cut vertices we visit along the way
	// if we encounter an already marked root vertex that isn't right at the start of our walk, then we report a tri-cut-comp

	N_LOG("scanning for bicomp with three cut vertices:\n")

	std::vector<int> prev_cut((size_t)(n_bicomps), -1); // prev_cut[i] is the previous cut vertex we saw when we hit the cut vertex with cut_verts value i, or -1 if there is none
	int root_one = -1; 									// two additional child bicomps for the root; the root bicomp is the only bicomp which can have two children in the bicomp tree if there's no tri-comp-cut, so we need some extra variables to keep track of these similar to root_cut
	int root_two = -1;

	for (int i = 0; i < n_bicomps - 1; i++) { // walk up from the root vertex of every bicomp to the next cut vertex along the DFS tree (except for the root bicomp, which has no parent)
		int w = retval[i].first;
		int u = -1;
		int start = w; // used to report a tri-cut-comp and update prev_cut

		while (w != root) {
			u = w;
			w = parent[w];
			V_LOG("walking up tree for bicomp " << i << ", w: " << w << ", u: " << u <<"\n")

			if (cut_verts[w] != -1 && u == retval[cut_verts[w]].second) { // if we backed into the root vertex of a bicomp from that bicomp
				V_LOG("found child bicomp: vertex " << start << " (bicomp " << i << ") child of vertex " << w << " (bicomp " << cut_verts[w] << ")\n")
				if (prev_cut[cut_verts[w]] == -1) { // if it doesn't already have a child, add this bicomp as its child
					prev_cut[cut_verts[w]] = start;
				} else { // otherwise, this bicomp has two children; report a tri-cut-comp
					std::shared_ptr<negative_cert_tri_cut_comp> cut{new negative_cert_tri_cut_comp{}};
					cut->c1 = w;
					cut->c2 = start;
					cut->c3 = prev_cut[cut_verts[w]];
					N_LOG("NON-SP, bicomp (not at root) with three cut vertices: " << cut->c1 << ", " << cut->c2 << ", " << cut->c3 << "\n")

					cert_out.sp_reason = cut;
					return retval;
				}
				break;
			}
		}

		if (w == root && (u == retval.back().second || u == -1)) { // special handling for the root bicomp, which is allowed to have two children before we report a tri-cut-comp
			V_LOG("found child bicomp of root: vertex " << start << " (bicomp " << i << ") child of vertex " << w << " (bicomp " << n_bicomps - 1 << ")\n")
			if (root_one == -1) {
				root_one = start;
			} else if (root_two == -1) {
				root_two = start;
			} else { // if the root bicomp has three children bicomps, report a tri-cut-comp
				std::shared_ptr<negative_cert_tri_cut_comp> cut{new negative_cert_tri_cut_comp{}};
				cut->c1 = root_one;
				cut->c2 = root_two;
				cut->c3 = start;
				N_LOG("NON-SP, bicomp (at root) with three cut vertices: " << cut->c1 << ", " << cut->c2 << ", " << cut->c3 << "\n")

				cert_out.sp_reason = cut;
				return retval;
			}
		}
	} // end of finding tri-cut-comp, at this point the bicomps are guaranteed to form a chain

	N_LOG("no tri-cut-comp found\n")

	if (n_bicomps > 1) { // only reorder if there are at least two bicomps to reorder
		N_LOG("ordering bicomps as chain: ")
		int second_endpoint = n_bicomps - 1;

		for (int i = 1; i < n_bicomps - 1; i++) { // locate the second bicomp with no child (i.e. no previous cut vertex)
			if (prev_cut[i] == -1) {
				second_endpoint = i;
				break;
			}
		}

		N_LOG("bicomp " << second_endpoint << " is the other bicomp with no child\n")

		std::reverse(retval.begin() + second_endpoint, retval.end() - 1); // swap the order of the bicomps after the one with no child, but don't include the root bicomp
		if (second_endpoint != n_bicomps - 1) { // if the root isn't already at the end
			retval.back().second = retval[n_bicomps - 2].first; // move the root bicomp over so that it points into the root vertex of its only child bicomp (fulfilling the case i = h in GSP-SP-OP in the paper)
			retval.back().first = retval[n_bicomps - 2].second;
		} else { // handle edge case where root is at end
			if (retval.back().first == retval[n_bicomps - 2].first) { // if the root bicomp's root vertex overlaps with its child, then just flip it around so it points into its child
				retval.back().first = retval.back().second;
			} else { // otherwise, the child of the root bicomp has a defined parent, which is in the root bicomp (and hence we can use the edge between the parent of the root vertex of the child bicomp and the root vertex of the child bicomp as the edge)
				retval.back().first = parent[retval[n_bicomps - 2].first];
			}
			retval.back().second = retval[n_bicomps - 2].first;
		}


		for (int i = second_endpoint; i < n_bicomps - 1; i++) { // invert the pointing directions of all cut vertices after the swap point, realigning the bicomps so they point in a chain
			retval[i].second = parent[retval[i].first];
		}

		#ifdef __VERBOSE_LOGGING__
			for (int i = 0; i < n_bicomps; i++) {
				V_LOG("bicomp " << i << " after reordering: root " << retval[i].first << ", edge " << retval[i].second << "\n")
			}
		#endif
	}

	return retval;
}


void report_K4_non_stack_pop_case(gsp_sp_op_result& cert_out,
								  std::vector<int> const& parent, 
								  std::vector<std::stack<sp_chain_stack_entry>>& vertex_stacks, 
								  int a, 
								  int b,
								  int d,
								  int elose,
								  int ewin_src,
								  int ewin_sink) { // reports a K4, for the more general case where a seq does not match its entire ear (rather than a stack popping violation)
												   // as in the other case, I order the 4 vertices of the K4 subdivision a, b, c, d from ancestrally latest to earliest
												   // ewin_src and ewin_sink are the source and sink of the back edges corresponding to the non-cut-off ear, so its tree path can be traced
												   // elose is the source of the back edge corresponding to the ear that got cut off, so its tree path can be traced
												   // I've already explained in more detail a lot of the processes I'm doing here in the other K4 subdivision reporting branch

	std::shared_ptr<negative_cert_K4> k4{new negative_cert_K4{}};
	k4->a = a;  // in this case, a is the source of the sequence of the vertex whose parent edge's ear got cut off (representing a link in the SP-chain for that ear); it is also the sink of the ear not s-attached to any other ear
	k4->b = b;  // b is the sink of the cut-off ear
				// c is the to-be-determined source of the ear not s-attached to any other ear, causing a 3.4a violation (in case C it could be a 3.4b violation, but everything works out the same regardless and we can arbitrarily consider the ear containing the sink of the violating ear to be the cut off one)
	k4->d = d;  // d is the source of the cut-off ear

	sp_tree earliest_violating_ear; // ear not s-attached to any other ear
	for (int b = parent[k4->b]; b != k4->d; b = parent[b]) { // walk up the tree-path from the sink of the cut-off ear to the source (not counting the source since then c would equal d), examining the vertex stacks as we go up to find an ear whose sink is at a, like we do in the other case
															 // technically any ear whose source is strictly ancestrally earlier than b and strictly ancestrally later than d should suffice here, but the paper says to take the lexicographically earliest one so I will do that
		for (; !vertex_stacks[b].empty(); vertex_stacks[b].pop()) {
			if (vertex_stacks[b].top().end == k4->a) { // if this ear sinks at a, its source s-belongs to the winning ear but its sink s-belongs to the losing ear, which is a violation
				earliest_violating_ear = std::move(vertex_stacks[b].top().SP); // note that we are moving up the DFS tree, so every time we go to a new candidate c the potential violating ear gets lexicographically earlier, so there's no need to perform any check to see if the ear is earlier and we can just immediately grab it
				k4->c = b;
			}
		}
	}


	for (int a = k4->a; a != k4->b; a = parent[a]) k4->ab.emplace_back(a, parent[a]); // now it's a straightforward path trace up the tree from ab to bc to cd
	for (int b = k4->b; b != k4->c; b = parent[b]) k4->bc.emplace_back(b, parent[b]);
	for (int c = k4->c; c != k4->d; c = parent[c]) k4->cd.emplace_back(c, parent[c]);

	k4->ad.emplace_back(k4->d, elose); // add the edge of the cut-off ear, which loops back around to a when we walk up the tree
	for (int d = elose; d != k4->a; d = parent[d]) k4->ad.emplace_back(d, parent[d]);
	for (int e = k4->d; e != ewin_src; e = parent[e]) {
		k4->bd.emplace_back(e, parent[e]); // follow the tree path up from the source of the cut-off ear until we meet the source of the non-cut-off ear (which may be the same), then loop down and around to b
	}
	k4->bd.emplace_back(ewin_src, ewin_sink);
	for (int e = ewin_sink; e != k4->b; e = parent[e]) k4->bd.emplace_back(e, parent[e]);
	int ear_path = earliest_violating_ear.underlying_tree_path_source(); // add the ear not s-attached to any ear, completing the sixth and final path
	k4->ac.emplace_back(k4->c, ear_path);
	for (; ear_path != k4->a; ear_path = parent[ear_path]) k4->ac.emplace_back(ear_path, parent[ear_path]);

	cert_out.gsp_reason = k4;
}

void K23_test(std::shared_ptr<certificate>& cert_ptr, std::vector<int>& alert, std::vector<int> const& parent, edge_t ear_found, edge_t ear_winning, int w) { // tests for K23, and puts the produced subdivision in cert_ptr if there is a K23
																																	// ear_found is the back-edge corresponding to the non-trivial ear we are testing for violation (note in the paper they just pass in an edge and index the ear array at that edge in this procedure, but I pass in the back-edge of that ear directly and index the array when calling)
																																	// ear_winning is the ear that cut off that ear (or, in case (c), chosen arbitrarily to be the other ear)
																																	// w is the sink of that ear
																																	// parent is the parent of w (so we can check if it's a 3.5(a) violation)
	V_LOG("testing K23: found ear (" << ear_found.first << ", " << ear_found.second << "), winning ear (" << ear_winning.first << ", " << ear_winning.second << ")\n")
	if (ear_found.second != parent[w]) { // the sink of the ear is at w, but its source isn't the parent of w, so there is a 3.5(a) violation
								// we report a K23 here
		N_LOG("OOPS, 3.5(a) violation, nonouterplanar\n")
		std::shared_ptr<negative_cert_K23> k23{new negative_cert_K23{}};
		k23->a = w;			       // the ends of the non-trivial ear we made are the two vertices on the "2" of the K23
		k23->b = ear_found.second;

		k23->one.emplace_back(k23->b, ear_found.first); // one path between these two vertices of length at least 2 is the non-trivial ear, found by tracing up the tree path from the source of its back edge to w
		for (int i = ear_found.first; i != k23->a; i = parent[i]) k23->one.emplace_back(i, parent[i]);

		for (int i = k23->a; i != k23->b; i = parent[i]) k23->two.emplace_back(i, parent[i]); // another path of length at least 2 is the tree path from the sink of the ear to its source (which is length at least 2 because of the ear_found.second != parent[w] check)

		for (int i = k23->b; i != ear_winning.second; i = parent[i]) k23->three.emplace_back(i, parent[i]); // the third and final path of length at least 2 forming the K23 is obtained by walking up the tree from the sink of the violating ear until we get to the source of the ear that cut it off and looping back around to the sink of the violating ear 
		k23->three.emplace_back(ear_winning.second, ear_winning.first);
		for (int i = ear_winning.first; i != k23->a; i = parent[i]) k23->three.emplace_back(i, parent[i]); 

		cert_ptr = k23;
		return;
	} // end of reporting K23

	if (alert[w] != -1) { // we have already had a non-trivial ear whose sink is at w, and both of these ears have a source equal to the parent of w, so there is a 3.5(b) violation
						  // we report a K23 here
		N_LOG("OOPS, 3.5(b) violation, nonouterplanar\n")
		std::shared_ptr<negative_cert_K23> k23{new negative_cert_K23{}};
		k23->a = w;			       // the ends of the non-trivial ear we made are the two vertices on the "2" of the K23
		k23->b = ear_found.second;

		k23->one.emplace_back(k23->b, ear_found.first); // path one is the same, it's the non-trivial ear we made
		for (int i = ear_found.first; i != k23->a; i = parent[i]) k23->one.emplace_back(i, parent[i]);

		k23->two.emplace_back(k23->b, alert[w]); // path two is the second non-trivial ear we found with sink at w and source at parent[w]
											// we already stored the source of the back edge of that ear in alert[w], so we can recover the relevant tree path
		for (int i = alert[w]; i != k23->a; i = parent[i]) k23->two.emplace_back(i, parent[i]);

		for (int i = k23->b; i != ear_winning.second; i = parent[i]) k23->three.emplace_back(i, parent[i]); // path three is the same, we walk up the tree from the source of the violating ear and loop back around to the sink using the winning ear
		k23->three.emplace_back(ear_winning.second, ear_winning.first);
		for (int i = ear_winning.first; i != k23->a; i = parent[i]) k23->three.emplace_back(i, parent[i]); 

		cert_ptr = k23;
		return;
	} else { // end of reporting K23
		// this is the first non-trivial ear sinked at w
		alert[w] = ear_found.first; // mark the back edge corresponding to that ear so we can report a K23 if we find another one
									// we only need to mark the source of the back edge since the sink of it is the parent of w
	}
}

int path_contains_edge(std::vector<edge_t> const& path, edge_t test) { // utility for finding if a path contains an edge, used for replacing K4's with T4's and fake edges in K23's
	for (size_t i = 0; i < path.size(); i++) {
		edge_t e = path[i];
		if (e == test || (e.first == test.second && e.second == test.first)) return (int)(i);
	}

	return -1;
}

#endif