# GSP-SP-OP

This is Nathan Levy's COMP-4990 final project at the University of Windsor, supervised by Dr. Tsin.
It is a C++ implementation of a certifying algorithm given in <a href="https://www.sciencedirect.com/science/article/abs/pii/S0166218X22003900?via%3Dihub">this paper</a> for determining if a given input graph is generalized series-parallel, and if it is, whether or not it is series-parallel and/or outerplanar. It also contains the relevant authentication algorithms for verifying the certificates. The algorithm and all the authentication algorithms run in O(|V| + |E|) time and space, and if I've done everything right this implementation should as well.

## Usage

The `src` directory contains a header file `gsp-sp-op.hxx`. To use the implementation, simply include this file. It defines a `graph` struct, which represents a graph the implementation can be run on.
To create a `graph`, read in from a `std::basic_istream` (e.g. file, standard input) using the `>>` operator. The input must have the following format:
* On the first line of the input, there are two numbers. The first is the number of vertices in the graph (n) and the second is the number of edges (e).
* On all of the next e lines of the input, there are two numbers from 0 to (n - 1) inclusive, which represent the two endpoints of an edge in the graph.

To execute the implementation on a `graph`, use `gsp_sp_op_result GSP_SP_OP(graph const& g)`. This returns a `gsp_sp_op_result` struct, which has three members:
* `is_gsp` is true if and only if the graph is generalized series-parallel
* `is_sp` is true if and only if the graph is series-parallel
* `is_op` is true if and only if the graph is outerplanar

Additionally, a `gsp_sp_op_result` contains three pointers to `certificate`s: `gsp_reason`, `sp_reason`, and `op_reason`. Seven types inherit from `certificate`:
* `positive_cert_gsp` represents a GSP decomposition tree, which shows that a graph is GSP. It also has a member `is_sp` if the decomposition tree is the special case of an SP decomposition tree (making the graph also SP).
* `positive_cert_op` represents the exterior boundary of an outerplanar embedding, which shows a graph is OP.
* `negative_cert_K4` represents a subdivision of K4, which shows a graph is not GSP, not SP, and not OP.
* `negative_cert_K23` represents a subdivision of K23, which shows a graph is not OP, but may still be GSP and/or SP.
* `negative_cert_T4` represents a subdivision of T4 (T4 is K4, but with an edge missing) with two cut vertices at the endpoints of the removed path, which shows a graph is not GSP, but may still be SP.
* `negative_cert_tri_comp_cut` represents a cut vertex contained in three distinct biconnected components, which shows a graph is not GSP, but may still be SP.
* `negative_cert_tri_cut_comp` represents a biconnected component with three cut vertices contaned in it, which shows a graph is not GSP, but may still be SP.

These certificates may be authenticated using `bool certificate::authenticate(graph const& g)`, to ensure they are well-formed and verify the result produced by the implementation. This will return `true` if and only if the authentication was successful. The graph passed into this function must be the same graph that generated the certificate.

You can use `bool gsp_sp_op_result::authenticate(graph const& g)` to authenticate all three of a result's `gsp_reason`, `sp_reason`, and `op_reason` for a given graph. Note that the pointers to `certificates` may point to the same certificate (e.g. if there is a K4 subdivision in the graph, then all three of `gsp_reason`, `sp_reason`, and `op_reason` will point to the same `negative_cert_K4`).

## Demo compilation and execution

`clang++ -std=c++20 -Wall -Wextra src/(filename) -o tester`

There are three files which I've provided to be compiled in this way:
### tester.cxx
tester.cxx takes a directory as its single command line argument. Every .txt file which is an immediate child of that directory will be converted into a graph, have the implementation executed on it, and then have the result authenticated. There are three directories of manually-constructed test cases given as examples:
* test cases/biconnected contains several small biconnected graphs. These have been manually constructed to hit several corner cases.
* test cases/non_biconnected contains several small non-biconnected graphs. These have also been manually constructed with corner cases in mind.
* test cases/massive contains two very large graphs (~1000000 vertices), to demonstrate that the implementation can handle them.

### random_tester.cxx
random_tester.cxx will generate 100000 random moderately-sized graphs (~400 vertices or so on average), run the implementation on all of them, and authenticate all the results. It might take around 30 seconds to run. If any random test fails to authenticate, it will stop and output the parameters used for the random graph generator, which can be given as command line arguments to recreate_random_failed_test.cxx to replicated the failed test case.

### recreate_random_failed_test.cxx
recreate_random_failed_test.cxx takes as command line arguments the parameters output by random_tester.cxx, and uses them to recreate the test case that failed. It re-runs the implementation on the given failed test case, and re-authenticates.

## Graph generator
In addition to the three demo files, I've also provided a file `GraphGenerator.hxx`, which may be included to generate random `graph`s for testing the implementation on. It defines a single function `graph generate_graph(long nC, long lC, long nK, long lK, long three_edges, long seed)`. The algorithm this function uses to generate a graph is equivalent to the following:
* First, generate `nC` cycle subgraphs on `lC` vertices and `nK` complete subgraphs on `lK` vertices.
* Label the vertices of these subgraphs from `0` to `nC * lC + nK * lK - 1`. Shuffle the labels.
* Randomly order all of these subgraps and number them from `1` to `nC + nK`.
* For `i` from `2` to `nC + nK`:
  * Choose a random subgraph from the `1`st to the `i-1`th.
  * Select two random distinct vertices in the `i`th subgraph (or three, if `three_edges` is true).
  * Connect each of these vertices to a randomly selected vertex in the chosen subgraph with a single edge, ensuring the two (or three) vertices such selected are distinct as well.
* Shuffle the order of every edge in the adjacency lists.

The parameters are as follows:
* `nC` is the number of cycle subgraphs in the generated output
* `lC` is the length of all of these cycles (must be at least 3)
* `nK` is the number of complete subgraphs in the generated output
* `lK` is the size of these complete subgraphs (must be at least 3)
* `three_egdes` is whether to connect each generated subgraph to the rest with three rather than two edges. A "no" is represented by 0, and anything else is a "yes".
* `seed` is the seed to be passed to the (C-style) random number generator. If not specified, it will use a "random" seed based off the current time.

It directly returns a `graph`, which can then be passed to `GSP_SP_OP`. Unlike the previous version, it does not accept command line arguments, and all arguments are instead passed to the above function. (This was done to make it easier to generate multiple graphs in one run of a test program).

## Logging
Some command-line options can be given to the compiler when the code is compiled to affect the implementation's logging behaviour, for debugging purposes:
* `-D__LIGHT_LOGGING__` will cause a minimal amount of extra information to be printed to standard output, which takes constant time to print.
* `-D__LOGGING__` will cause some extra information to be printed, including some details about what the implementation is doing (which may take time linear in the size of the input to print).
* `-D__VERBOSE_LOGGING__` will cause a very large amount of extra information to be printed and show what the implementation is doing in excruciating detail for debugging purposes. It is not recommended for large graphs, as printing all this information takes time quadratic in the size of the input.
* `-D__DEBUG_LOGGING__` is unused and does nothing. If you ever want to debug the code, though, a statement printed out with debug logging will flush standard output afterwards (so you can still see what went wrong if there's a segfault).

## Caveats
* The implementation assumes the graph is simple (contains no multiple edges or self-loops) and connected, and will fail if this is not the case. It would be trivial to modify it to handle disconnected/multigraphs, so I have been told not to bother doing this
* The random graph generator only generates biconnected graphs (since it connects every subgraph it generates to the rest of the graph with two edges). This means the code is pretty poorly tested on non-biconnected graphs, and some bugs may still exist in the part of the code which handles non-biconnected graphs. The implementation works for every non-biconnected corner case I could come up with, though, and I've gone ahead and manually modified a few big randomly generated graphs to be non-biconnected, so it should hopefully be fine.
