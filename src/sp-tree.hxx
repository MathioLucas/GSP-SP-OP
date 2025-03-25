// this file contains the definition of an SP tree, a SP chain stack entry, and associated operations on SP trees

#ifndef __SP_TREE_HXX__
#define __SP_TREE_HXX__

#include "logging.hxx"
#include <ostream>
#include <stack>

enum class c_type { // composition type
	edge, series, parallel, antiparallel, dangling // an antiparallel composition is the exact same as a parallel composition, but the right node is flipped around
												   // the source of the left subgraph is equal to the sink of the right one, the sink of the left one is the source of the right one, and the overall source and sink are equal to the source and sink of the left one 
												   // I use this to "mark" nodes, as the paper suggests, in the authentication algorithm for SP subgraphs, but other than that antiparallel nodes function identically to parallel ones
}; 

char c_type_char(c_type comp) { // get the character corresponding to a composition type, for display purposes
	switch (comp) {
		case c_type::edge:
			return 'e';
		case c_type::series:
			return 'S';
		case c_type::parallel:
			return 'P';
		case c_type::antiparallel:
			return 'Q';
		case c_type::dangling:
			return 'D';
	}
}

struct sp_tree_node { // a single node in the sp tree
	int source;
	int sink;
	sp_tree_node * l;			// pointers to left and right subtrees
	sp_tree_node * r;
	c_type comp;			// the type of composition this node represents, or c_type::edge if it's not a composition but a single edge

	sp_tree_node(int source_, int sink_) : source{source_}, sink{sink_}, comp{c_type::edge} {} // construct edge node

	sp_tree_node(sp_tree_node * l_, sp_tree_node * r_, c_type comp_) : l{l_}, r{r_}, comp{comp_} {  // construct composed node
		switch (comp) {
			case c_type::series:
				source = l->source;
				sink = r->sink;
				break;
			case c_type::dangling:
			case c_type::parallel:
			case c_type::antiparallel:
				source = l->source;
				sink = l->sink;
				break;
			case c_type::edge:
				break;
		}
	}
};

struct sp_tree { // an SP decomposition tree
	sp_tree_node * root;

	void compose(sp_tree&& other, c_type comp) { // compose two trees, with the other tree going on the right of this one
		if (!root) {
			root = other.root;
			other.root = nullptr;
			return;
		} else if (!other.root) {
			return;
		}

		root = new sp_tree_node{root, other.root, comp};
		other.root = nullptr;
	}

	void l_compose(sp_tree&& other, c_type comp) { // compose two trees, with the other tree going on the left of this one
		if (!root) {
			root = other.root;
			other.root = nullptr;
			return;
		} else if (!other.root) {
			return;
		}

		root = new sp_tree_node{other.root, root, comp};
		other.root = nullptr;
	}

	void deantiparallelize() { // removes all the antiparallel compositions in an SP tree by swapping around the sources and sinks of their right children, converting the antiparallel compositions to parallel ones
							   // runs in O(|V| + |E|) time
							   // the paper suggests keeping a "swap" switch and doing this while verifying the decomposition tree, but the approach here works equally well (and might be nicer since antiparallel compositions don't actually exist in SP graphs)
		std::stack<std::pair<sp_tree_node *, int>> hist;
		bool swap = false;

		if (!root) return;

		hist.emplace(root, 0);

		while (!hist.empty()) {
			sp_tree_node * curr = hist.top().first;

			if (hist.top().second == 0) {
				hist.top().second++;
				if (curr->r) hist.emplace(curr->r, 0);

				if (curr->comp == c_type::antiparallel) swap = !swap;
			} else {
				hist.pop();
				if (curr->l) hist.emplace(curr->l, 0);

				if (curr->comp == c_type::antiparallel) {
					swap = !swap;
					curr->comp = c_type::parallel;
				}

				if (swap) {
					sp_tree_node * temp = curr->l;
					curr->l = curr->r;
					curr->r = temp;
					int temp_src = curr->source;
					curr->source = curr->sink;
					curr->sink = temp_src;
				}
			}
		}
	}

	int source() {return root->source;}
	int sink() {return root->sink;}
	int underlying_tree_path_source() { // get the source of a tree path which, when combined with an outgoing back edge, forms a path between the source and sink of this SP tree
									    // this is only used to generate subdivisions, and only meaningful for the .SPs of stack entries, since the .tails aren't complete ears but only sections of ears (and their source edges are not guaranteed to be back edges)
									    // I do this by traversing the SP tree which is O(|E|) time rather than O(|V|) time as the paper suggests generating a K4 or K2,3 should be
									    // I could store the back edge corresponding to the underlying ear as a data member of the tree to be able to find this in O(1) time, but that would be a waste 
									    // it'd take O(|E|) time to keep that data member up to date anyway as we went through the algorithm, since we'd need to update it with every SP tree we create and we create one per edge
		sp_tree_node * leftmost = root;
		for (; leftmost->comp != c_type::edge; leftmost = leftmost->l);
		return leftmost->sink;
	}

	sp_tree() {
		root = nullptr;
	}
	sp_tree(int source_, int sink_) : root{new sp_tree_node{source_, sink_}} {} // construct edge from source and sink

	~sp_tree();

	sp_tree(sp_tree const& other) = delete; // should never copy an sp tree, that is not O(1) time; we should only move them
											// this way the compiler will complain if I ever attempt to copy an sp tree, preventing me from accidentally making the implementation non-linear time

	sp_tree& operator=(sp_tree const& other) = delete;

	sp_tree(sp_tree&& other) {
		root = other.root;
		other.root = nullptr;
	}


	sp_tree& operator=(sp_tree&& other) {
		if (this != &other) {
			delete root;
			root = other.root;
			other.root = nullptr;	
		}
		return *this;
	}
};

std::ostream& operator<<(std::ostream& os, sp_tree_node const& t) {
	#ifdef __VERBOSE_LOGGING__
	os << "{";
	if (t.l) os << *(t.l);
	os << t.source << c_type_char(t.comp) << t.sink;
	if (t.r) os << *(t.r);
	os << "}";
	#else
	os << "{" << t.source << c_type_char(t.comp) << t.sink << "}";
	#endif

	return os;
}

std::ostream& operator<<(std::ostream& os, sp_tree const& t) { // output an SP tree (for debugging purposes)
	if (t.root) {
		os << *(t.root);
	} else {
		os << "(null tree)";
	}
	return os;
}

sp_tree::~sp_tree() {
	if (!root) return;

	std::stack<std::pair<sp_tree_node *, int>> hist;
	hist.emplace(root, 0);

	while (!hist.empty()) {
		sp_tree_node * curr = hist.top().first;

		if (hist.top().second == 0) {
			hist.top().second = 1;
			if (curr->r) hist.emplace(curr->r, 0);
			if (curr->l) hist.emplace(curr->l, 0);
		} else {
			delete curr;
			hist.pop();
		}
	}
}

struct sp_chain_stack_entry {
	sp_tree SP;   // an ear with source y and sink x and all ears s*-attached to that ear, represented as an SP decomposition tree SP_(x, y)
				    // note that though the ear has source y and sink x, the produced SP decomposition has source x and sink y
				    // y is always equal to the vertex whose stack this entry is stored on and x is equal to this entry's end

	int end;        // the end of the above ears

	sp_tree tail; // a path with source z and sink x making up a section of an ear and all ears s*-attached to that section, represented as an SP decomposition tree SP_(z -> x)
				    // x is equal to this entry's end, and z is some other vertex

	sp_chain_stack_entry(sp_tree SP_, int end_, sp_tree tail_) : SP{std::move(SP_)}, end{end_}, tail{std::move(tail_)} {} // construct stack entry
																														  	  // this invalidates its arguments afterwards
	sp_chain_stack_entry() = default;
};

#endif