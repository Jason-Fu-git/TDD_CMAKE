//
// Created by Jason Fu on 24-7-14.
//

#ifndef TDD_C_CONTRACTIONOPTIMIZER_HPP
#define TDD_C_CONTRACTIONOPTIMIZER_HPP

#include "ContractionTree.hpp"
#include "Tdd.hpp"
#include "Graph.hpp"
#include "Package.hpp"
#include "Export.hpp"


#include <string>
#include <map>
#include <regex>
#include <cassert>
#include <random>
#include <ctime>
#include <vector>
#include <libkahypar.h>

using GateSet = std::map<int, gate>;
using IndexSet = std::map<int, std::vector<dd::Index>>;
using Node = ContractionTree::Node;
using Graph = dd::Graph;
using GraphNode = dd::Graph::GraphNode;
using GraphEdge = dd::Graph::GraphEdge;

class ContractionOptimizer {

public:

    /**
     * Construct an optimizer with the given gate set and index set;
     * @param num_qubits The number of qubits in the circuit;
     * @param gates The pointer to the gate set;
     * @param indexes The pointer to the index set;
     * @param _index_width The width of the index, default is 2;
     * @author Jason Fu
     */
    explicit ContractionOptimizer(int num_qubits, GateSet *gates, IndexSet *indexes, int _index_width = 2)
            : qubits_num(num_qubits), gate_set(gates), index_set(indexes), index_width(_index_width) {}

    virtual ~ContractionOptimizer() = default;

    /**
     * Optimize the contraction order of the given circuit using the given method;
     * @author Jason Fu
     */
    virtual ContractionTree *optimize() = 0;

    /**
     * Perform contraction on the given contraction tree;
     * @return a size-2 array. The first element is the maximum number of nodes in the TDD during the simulation,
     * and the second element is the final number of nodes in the TDD.
     * @author Jason Fu
     */
    int *contract(ContractionTree *tr, std::unique_ptr<dd::Package<>> &dd, bool release, bool print = false) {
        auto gateSet = *gate_set;
        auto indexSet = *index_set;
        int max_nodes = 0, final_nodes = 0;
        std::vector<Node *> stack = {tr->getRoot()};
        auto tdd = contractNode(stack, max_nodes, final_nodes, dd, release, print);
//        auto tdd1 = contractNode(tr->getRoot(), max_nodes, final_nodes, dd, release);
        // do something with the tdd
//        dd::serialize(tdd.e, std::cout, false);
        std::cout << tdd.e.w << std::endl;
//        dd::TDD tdd2 = dd->cont(tdd, tdd1);
//        std::cout << dd->size(tdd2.e) << std::endl;
//        std::cout << tdd2.e.w.r->value << " " << tdd2.e.w.i->value << std::endl;
        // delete it
        if (release)
            dd->decRef(tdd.e);
        dd->garbageCollect();
        // return the statistics
        return new int[2]{max_nodes, final_nodes};
    }

protected:
    // parse the circuit into a graph
    Graph *constructGraph();

    // contract using stack
    dd::TDD
    contractNode(std::vector<Node *> &stack, int &max_nodes, int &final_nodes, std::unique_ptr<dd::Package<>> &dd,
                 bool release, bool print = false);

    // construct a TDD from a gate
    // the parser is copied from Cir_import.h
    static dd::TDD
    constructGate(std::string nam, const std::vector<dd::Index> &indexSet, std::unique_ptr<dd::Package<>> &dd);

    // string utils, copied from Cir_import.h

    static float match_a_string(std::string s);

    static std::vector<std::string> split(const std::string &s, const std::string &seperator);

    inline static int unionFind(const std::vector<int> &unionFindSet, int s);

    int qubits_num;
    int index_width;
    GateSet *gate_set;
    IndexSet *index_set;
};

/**
 * perform exhaustive search on the optimal contraction tree
 * using dynamic programming and only considering connected sub-graphs
 *
 * Recommended for small circuits (e.g. gate set size < 18, which can be calculated in 10 minutes)
 * @see https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.033315
 */
class ExhaustiveSearchOptimizer : public ContractionOptimizer {
public:

    // Bitset used in bipartite process
    struct Bitset {
        unsigned long long *values;
        int length;

        // default constructor
        explicit Bitset(int length) : length(length) {
            values = new unsigned long long[length / 64 + 1];
            memset(values, 0, (length / 64 + 1) * sizeof(unsigned long long));
        }

        // destructor
        ~Bitset() { delete[] values; }

        // check if the bit at position i is set
        bool test(int i) const { return values[i / 64] & (1ULL << (i % 64)); }

        // set the bit at position i
        void set(int i) { values[i / 64] |= (1ULL << (i % 64)); }

        // clear the bit at position i
        void clear(int i) { values[i / 64] &= ~(1ULL << (i % 64)); }

        // add 1
        void add1() {
            int i = 0;
            while (values[i] == ULLONG_MAX && i < length / 64 + 1) {
                values[i++] = 0;
            }
            if (i == length / 64 + 1) {
                throw std::runtime_error("Bitset overflow");
            }
            values[i]++;
        }

        // print
        void print() const { for (int i = 0; i < length; i++) std::cout << (short) test(i); }
    };

    explicit ExhaustiveSearchOptimizer(int num_qubits, GateSet *gates, IndexSet *indexes, int _index_width = 2)
            : ContractionOptimizer(
            num_qubits, gates, indexes, _index_width) {
        // construct a leaf node for all the gates in the gate set
        for (int i = 0; i < gate_set->size(); i++) {
            auto idxVec = index_set->at(i);
            auto *idxSet = new std::set<Index>();
            for (auto &idx: idxVec) {
                idxSet->insert(idx);
            }
            auto node = new Node(idxSet, i);
            assert(idxSet->size() == idxVec.size());
            delete idxSet;
            node_set.push_back(node);
        }
    }

    ContractionTree *optimize() override;

    ~ExhaustiveSearchOptimizer() override {
        for (auto node: node_set)
            delete node;
        node_set.clear();
    }

private:
    // find the optimal partition of the given node, using dynamic programming
    ContractionTree *findOptimalPartition(std::vector<int> &nodeIdxs, unsigned long long minCost);

    // find the temped optimal partition from the dp_table
    ContractionTree *searchDPTable(std::vector<int> &nodeIdxs);

    // dump the current result to the dp_table
    void dumpDPTable(std::vector<int> &nodeIdxs, ContractionTree *result);

    unsigned long long count = 0;

    // the DP table
    std::map<std::string, ContractionTree *> dp_table;

    // the node set
    std::vector<Node *> node_set;
};

/**
 * Partition scheme 1 : Horizontally split the circuit. Introduce a vertical cut until the horizontal cut
 * of CNOTs has reached cx_cut_max.
 * @note This partition scheme does not cut the CNOTs. It puts the entire CNOT into one block instead.
 * @see https://arxiv.org/abs/2009.02618
 * @author Jason Fu
 */
class PartitionScheme1Optimizer : public ContractionOptimizer {
public:
    explicit PartitionScheme1Optimizer(int num_qubits, GateSet *gates, IndexSet *indexes, int cx_cut_max,
                                       int _index_width = 2)
            : ContractionOptimizer(num_qubits, gates, indexes, _index_width), max_cx_cut(cx_cut_max) {}

    ~PartitionScheme1Optimizer() override = default;

    ContractionTree *optimize() override;

private:
    int max_cx_cut;
};

/**
 * Partition scheme 2 : Horizontally split the circuit. When the number of CNOTs across the cut reaches cx_cut_max,
 * introduce a small block to contain the following CNOTs. When the width of the block reaches c_part_width,
 * introduce a vertical cut.
 * @note This partition scheme does not cut the CNOTs. It puts the entire CNOT into one block instead.
 * @see https://arxiv.org/abs/2009.02618
 * @author Jason Fu
 */
class PartitionScheme2Optimizer : public ContractionOptimizer {
public:
    explicit PartitionScheme2Optimizer(int num_qubits, GateSet *gates, IndexSet *indexes,
                                       int cx_cut_max, int c_part_width, int _index_width = 2)
            : ContractionOptimizer(num_qubits, gates, indexes, _index_width),
              cx_cut_max(cx_cut_max), c_part_width(c_part_width) {}

    ~PartitionScheme2Optimizer() override = default;

    ContractionTree *optimize() override;

private:
    int cx_cut_max;
    int c_part_width;
};

/**
 * Perform community contraction first. The community detection algorithm is by Girvan and Newman.
 * @see https://www.pnas.org/doi/full/10.1073/pnas.122653799
 * @author Jason Fu
 */
class GNCommunityOptimizer : public ContractionOptimizer {
public:
    explicit GNCommunityOptimizer(int num_qubits, GateSet *gates, IndexSet *indexes, int _index_width = 2)
            : ContractionOptimizer(num_qubits, gates, indexes, _index_width) {
    }

    ContractionTree *optimize() override;

private:
    // build the contraction tree from the graph (GN algorithm)
    ContractionTree *buildTreeFromGraph(const std::vector<GraphEdge> &edgeOrder);

    // build the contraction tree from the graph (fast GN algorithm)
    ContractionTree *buildTreeFromGraphFast(Graph *graph);
};

/**
 * Optimizer based on a bottom-up greedy scheme.
 * Considering that outer product is rarely beneficial, we only consider nodes that have indices in common.
 * @see https://quantum-journal.org/papers/q-2021-03-15-410/
 * @author Jason Fu
 */
class GreedyOptimizer : public ContractionOptimizer {
public:
    struct Matrix {
        int rows;
        int cols;
        double *data;

        Matrix(int rows, int cols) : rows(rows), cols(cols) {
            data = new double[rows * cols];
            memset(data, 0, sizeof(double) * rows * cols);
        }

        ~Matrix() {
            delete[] data;
        }

        inline double at(int i, int j) const {
            return data[i * cols + j];
        }

        inline void set(int i, int j, double value) {
            data[i * cols + j] = value;
        }

        inline void setNull(int i, int j) {
            set(i, j, 1e20);
        }
    };

    struct Edge {
        int u;
        int v;
        double w;

        bool operator<(const Edge &e) const {
            return v < e.v;
        }

        bool operator==(const Edge &e) const {
            return u == e.u && v == e.v;
        }
    };

    explicit GreedyOptimizer(int num_qubits, GateSet *gates, IndexSet *indexes, double alpha, double tau,
                             int _index_width = 2)
            : ContractionOptimizer(num_qubits, gates, indexes, _index_width), alpha(alpha), tau(tau) {
    }

    ContractionTree *optimize() override;

private:
    /**
     * cost(lc, rc) = SIZE(lc * rc) - alpha * (SIZE(lc) + SIZE(rc))
     * @return the Boltzmann weight
     */
    inline double calculateCost(Node *lc, Node *rc);

    static inline std::pair<int, int> getBestPair(const std::vector<std::set<Edge> > &m);;

    static inline bool isAdjacent(Node *lc, Node *rc);

    // hyper parameters
    double alpha;
    double tau;
    // random generator
    std::default_random_engine randGen;
};

/**
 * Contraction Optimizer based on Kahypar
 * @see https://kahypar.org/
 * @author Jason Fu
 *
 */
class KahyparOptimizer : public ContractionOptimizer {
public:
    explicit KahyparOptimizer(int num_qubits, GateSet *gates, IndexSet *indexes,
                              double imbalance, kahypar_partition_id_t k,
                              int _index_width = 2)
            : ContractionOptimizer(num_qubits, gates, indexes, _index_width),
              imbalance(imbalance), k(k) {
    }

    ContractionTree *optimize() override;

private:
    ContractionTree *partition(const std::vector<int> &nodeIndexes,
                               kahypar_context_t *context);

    double imbalance;
    kahypar_partition_id_t k;
};


#endif //TDD_C_CONTRACTIONOPTIMIZER_HPP
