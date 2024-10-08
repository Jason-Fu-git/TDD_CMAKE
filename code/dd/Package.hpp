#pragma once

#include "Complex.hpp"
#include "ComplexCache.hpp"
#include "ComplexNumbers.hpp"
#include "ComplexTable.hpp"
#include "ComplexValue.hpp"
#include "ComputeTable.hpp"
#include "Control.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"
#include "GateMatrixDefinitions.hpp"
#include "Package_fwd.hpp"

#include "UniqueTable.hpp"

#include "Tdd.hpp"
#include "Tensor.hpp"


#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <numeric>

#include <xtensor/xarray.hpp>
#include <xtensor/xshape.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xslice.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>

namespace dd {

    /**
     * Packs basic TDD operations, like generation, normalization, addition and contraction.
     * Memory use will also be handled here.
     *
     */
    template<class Config>
    class Package {
        static_assert(std::is_base_of_v<DDPackageConfig, Config>, "Config must be derived from DDPackageConfig");

        ///
        /// Complex number handling
        ///
    public:
        ComplexNumbers cn{};

        ///
        /// Construction, destruction, information and reset
        ///

        static constexpr std::size_t MAX_POSSIBLE_QUBITS =
                static_cast<std::make_unsigned_t<Qubit>>(std::numeric_limits<Qubit>::max()) + 1U;
        static constexpr std::size_t DEFAULT_QUBITS = 300;


        bool to_test = false;

        std::map<std::string, int> varOrder;


        explicit Package(std::size_t nq = DEFAULT_QUBITS) : nqubits(nq) {
            resize(nq);
        };

        ~Package() = default;

        Package(const Package &package) = delete;

        Package &operator=(const Package &package) = delete;

        // resize the package instance
        void resize(std::size_t nq) {
            if (nq > MAX_POSSIBLE_QUBITS) {
                throw std::invalid_argument("Requested too many qubits from package. "
                                            "Qubit datatype only allows up to " +
                                            std::to_string(MAX_POSSIBLE_QUBITS) +
                                            " qubits, while " + std::to_string(nq) +
                                            " were requested. Please recompile the "
                                            "package with a wider Qubit type!");
            }
            nqubits = nq;
            nodeUniqueTable.resize(nqubits);
        }

        // reset package state
        void reset() {
            clearUniqueTables();
            clearComputeTables();
            cn.clear();
        }

        // getter for qubits
        [[nodiscard]] auto qubits() const { return nqubits; }

    private:
        std::size_t nqubits;

        ///
        /// Vector nodes, edges and quantum states
        ///
    public:

        Edge<mNode> xarray_2_edge(
                const xt::xarray<ComplexValue> &array,
                std::vector<std::size_t> *order = nullptr
        ) {
            std::size_t sum_of_dim = std::accumulate(array.shape().begin(),
                                                     array.shape().end(),
                                                     0);

            if (sum_of_dim == array.dimension()) {
                std::vector<Edge<mNode>> edges;
                for (auto num: array) {
                    if (num == complex_zero) {
                        edges.push_back(Edge<mNode>::zero);
                    } else if (num == complex_one) {
                        edges.push_back(Edge<mNode>::one);
                    } else {
                        edges.push_back(Edge<mNode>::terminal(cn.lookup(num)));
                    }
                }
                return makeDDNode(0, edges, false);

            }

            if (order == nullptr) {
                // list(range(dim))
                order = new std::vector<std::size_t>(array.dimension());
                std::iota(order->begin(), order->end(), 0);
            } else if (order->empty()) {
                order->reserve(array.dimension());
                std::iota(order->begin(), order->end(), 0);
            }

            auto split_pos = std::max_element(order->begin(), order->end()) - order->begin();
            Qubit x = (*order)[split_pos];
            (*order)[split_pos] = -1;
            std::vector<xt::xarray<ComplexValue>> split_U = xt::split(array, array.shape(split_pos), split_pos);


            std::vector<Edge<mNode>> edges;
            for (auto u: split_U) {
                edges.push_back(xarray_2_edge(u, order));
            }

            return makeDDNode((Qubit) x + 1, edges, false);

        }

        /**
         * Extract string keys from index set `var`
         * @param var
         * @return
         */
        std::vector<std::string> generate_key(std::vector<Index> var) {
            std::vector<std::string> res(var.size());
            for (auto &index: var) {
                res.push_back(index.key);
            }
            return res;
        }

        TDD Tensor_2_TDD(const Tensor tn) {
            if (tn.data.dimension() != tn.index_set.size()) {
                throw "action non definies";
            }

            TDD res;
            res.e = xarray_2_edge(tn.data);
            res.index_set = tn.index_set;
            res.key_2_index = generate_key(tn.index_set);

            return res;
        }


        TDD Matrix2TDD(const GateMatrix mat, std::vector<Index> var_out) {

            TDD low, high, res;
            Edge<mNode> e_temp[4];

            std::vector<Edge<mNode>> e_low(2), e_high(2), e(2);

            int Radix = 2;

            // 00, 01, 10, 11
            // convert complex values to edges
            for (int i = 0; i < Radix; i++) {
                for (int j = 0; j < Radix; j++) {
                    if (mat[2 * i + j] == complex_zero) {
                        e_temp[i * Radix + j] = Edge<mNode>::zero;
                    } else if (mat[2 * i + j] == complex_one) {
                        e_temp[i * Radix + j] = Edge<mNode>::one;
                    } else {
                        e_temp[i * Radix + j] = Edge<mNode>::terminal(cn.lookup(mat[2 * i + j]));
                    }
                }
            }


            std::vector<std::string> key_2_index;
            if (varOrder[var_out[0].key] < varOrder[var_out[1].key]) {
                // first split index var_out[0]
                e_low[0] = e_temp[0];
                e_low[1] = e_temp[1];
                e_high[0] = e_temp[2];
                e_high[1] = e_temp[3];
                key_2_index = {var_out[0].key, var_out[1].key};
            } else {
                // first split index var_out[1]
                e_low[0] = e_temp[0];
                e_low[1] = e_temp[2];
                e_high[0] = e_temp[1];
                e_high[1] = e_temp[3];
                key_2_index = {var_out[1].key, var_out[0].key};
            }


            if (e_low[0].p == e_low[1].p and e_low[0].w == e_low[1].w) {
                low.e = e_low[0];
            } else {
                low.e = makeDDNode(0, e_low, false);
            }
            if (e_high[0].p == e_high[1].p and e_high[0].w == e_high[1].w) {
                high.e = e_high[0];
            } else {
                high.e = makeDDNode(0, e_high, false);
            }
            if (low.e.p == high.e.p and low.e.w == high.e.w) {
                res.e = low.e;
            } else {
                e[0] = low.e;
                e[1] = high.e;
                res.e = makeDDNode(1, e, false);;
            }
            res.index_set = var_out;
            res.key_2_index = key_2_index;
            return res;
        }

        TDD diag_matrix_2_TDD(const GateMatrix mat, std::vector<Index> var_out) {

            TDD res;
            int Radix = 2;

            std::vector<Edge<mNode>> e_temp(2);
            for (int i = 0; i < Radix; i++) {

                if (mat[2 * i + i] == complex_zero) {
                    e_temp[i] = Edge<mNode>::zero;
                } else {
                    if (mat[2 * i + i] == complex_one) {
                        e_temp[i] = Edge<mNode>::one;
                    } else {
                        e_temp[i] = Edge<mNode>::terminal(cn.lookup(mat[2 * i + i]));
                    }
                }
            }

            res.e = makeDDNode(0, e_temp, false);
            res.index_set = var_out;
            res.key_2_index = {var_out[0].key};
            return res;
        }

        /**
         * Convert CNOT gate to a TDD, given the index set var.
         * @param var
         * @param ca
         * @return
         */
        TDD cnot_2_TDD(std::vector<Index> var, int ca = 1) {


            TDD low, high, res;
            std::vector<Edge<mNode>> e(2);
            if (ca == 1) {
                if (varOrder[var[0].key] > varOrder[var[3].key] && varOrder[var[0].key] > varOrder[var[4].key]) {
                    low = Matrix2TDD(Imat, {var[3], var[4]});
                    high = Matrix2TDD(Xmat, {var[3], var[4]});
                    e[0] = low.e;
                    e[1] = high.e;
                    res.e = makeDDNode(2, e, false);
                    res.index_set = {var[0], var[2], var[3], var[4]};
                    low.key_2_index.push_back(var[0].key);
                    res.key_2_index = low.key_2_index;
                } else if (varOrder[var[3].key] > varOrder[var[0].key] && varOrder[var[3].key] > varOrder[var[4].key]) {
                    low = Matrix2TDD(Imat, {var[0], var[4]});
                    high = Matrix2TDD(Xmat, {var[0], var[4]});
                    e[0] = low.e;
                    e[1] = high.e;
                    res.e = makeDDNode(2, e, false);
                    res.index_set = {var[0], var[2], var[3], var[4]};
                    low.key_2_index.push_back(var[3].key);
                    res.key_2_index = low.key_2_index;
                } else {
                    low = Matrix2TDD(Imat, {var[0], var[3]});
                    high = Matrix2TDD(Xmat, {var[0], var[3]});
                    e[0] = low.e;
                    e[1] = high.e;
                    res.e = makeDDNode(2, e, false);
                    res.index_set = {var[0], var[2], var[3], var[4]};
                    low.key_2_index.push_back(var[4].key);
                    res.key_2_index = low.key_2_index;
                }
                return res;
            }
            return res;

        }

        ///
        /// Matrix nodes, edges and quantum gates
        ///
        template<class Node>
        Edge<Node> normalize(const Edge<Node> &e, bool cached) {

            auto argmax = -1;

            int R = e.p->e.size();
            std::vector<bool> zero;
            for (int k = 0; k < R; k++) {
                zero.push_back(e.p->e[k].w.approximatelyZero());
            }

            // make sure to release cached numbers approximately zero, but not exactly zero
            if (cached) {
                for (auto i = 0U; i < R; i++) {
                    if (zero[i] && e.p->e[i].w != Complex::zero) {
                        cn.returnToCache(e.p->e[i].w);
                        e.p->e[i] = Edge<Node>::zero;
                    }
                }
            }

            fp max = 0;
            auto maxc = Complex::one;
            // determine max amplitude
            for (auto i = 0U; i < R; ++i) {
                if (zero[i]) {
                    continue;
                }
                if (argmax == -1) {
                    // argmax has not been initialized
                    argmax = static_cast<decltype(argmax)>(i);
                    max = ComplexNumbers::mag2(e.p->e[i].w);
                    maxc = e.p->e[i].w;
                } else {
                    // compare with existing max
                    auto mag = ComplexNumbers::mag2(e.p->e[i].w);
                    if (mag - max > ComplexTable<>::tolerance()) {
                        argmax = static_cast<decltype(argmax)>(i);
                        max = mag;
                        maxc = e.p->e[i].w;
                    }
                }
            }

            // all equal to zero
            if (argmax == -1) {
                if (!cached && !e.isTerminal()) {
                    // If it is not a cached computation, the node has to be put back into the chain
                    getUniqueTable<Node>().returnNode(e.p);
                }
                return Edge<Node>::zero;
            }

            auto r = e;
            // multiply the weight of the edge by maxc
            // divide each entry by maxc
            for (auto i = 0U; i < R; ++i) {
                if (static_cast<decltype(argmax)>(i) == argmax) {
                    if (cached) {
                        if (r.w.exactlyOne()) {
                            r.w = maxc;
                        } else {
                            ComplexNumbers::mul(r.w, r.w, maxc);
                        }
                    } else {
                        if (r.w.exactlyOne()) {
                            r.w = maxc;
                        } else {
                            auto c = cn.getTemporary();
                            ComplexNumbers::mul(c, r.w, maxc);
                            r.w = cn.lookup(c);
                        }
                    }
                    r.p->e[i].w = Complex::one;
                } else {
                    if (zero[i]) {
                        if (cached && r.p->e[i].w != Complex::zero) {
                            cn.returnToCache(r.p->e[i].w);
                        }
                        r.p->e[i] = Edge<Node>::zero;
                        continue;
                    }
                    if (cached && !zero[i] && !r.p->e[i].w.exactlyOne()) {
                        cn.returnToCache(r.p->e[i].w);
                    }
                    if (r.p->e[i].w.approximatelyOne()) {
                        r.p->e[i].w = Complex::one;
                    }
                    auto c = cn.getTemporary();
                    ComplexNumbers::div(c, r.p->e[i].w, maxc);
                    r.p->e[i].w = cn.lookup(c);
                }
            }
            return r;

        }


    private:

        ///
        /// Unique tables, Reference counting and garbage collection
        ///
    public:
        // unique tables
        template<class Node>
        [[nodiscard]] auto &getUniqueTable() {
            return nodeUniqueTable;
        }

        template<class Node>
        void incRef(const Edge<Node> &e) {
            getUniqueTable<Node>().incRef(e);
        }

        template<class Node>
        void decRef(const Edge<Node> &e) {
            getUniqueTable<Node>().decRef(e);
        }

        UniqueTable<mNode, Config::UT_MAT_NBUCKET, Config::UT_MAT_INITIAL_ALLOCATION_SIZE> nodeUniqueTable{nqubits};

        bool garbageCollect(bool force = false) {
            // return immediately if no table needs collection
            if (!force &&
                !nodeUniqueTable.possiblyNeedsCollection() &&
                !cn.complexTable.possiblyNeedsCollection()) {
                return false;
            }

            auto cCollect = cn.garbageCollect(force);
            if (cCollect > 0) {
                // Collecting garbage in the complex numbers table requires collecting the
                // node tables as well
                force = true;
            }

            auto mCollect = nodeUniqueTable.garbageCollect(force);

            // invalidate all compute tables where any component of the entry contains
            // numbers from the complex table if any complex numbers were collected
            if (mCollect > 0) {

                addTable.clear();
                contTable.clear();
            }
            return mCollect > 0;
        }

        void clearUniqueTables() {
            nodeUniqueTable.clear();

        }

        /**
         *\brief create a normalized, reduced DD node and return an edge pointing to it. The node will
        * not be recreated if it already exists.
        */
        template<class Node>
        Edge<Node> makeDDNode(
                Qubit var,
                const std::vector<Edge<Node>> &edges,
                bool cached = false) {

            auto &uniqueTable = getUniqueTable<Node>();
            Edge<Node> e{uniqueTable.getNode(), Complex::one};
            e.p->v = var;
            e.p->e = edges;

            assert(e.p->ref == 0);

            bool all_equal = true;
            for (int k = 0; k < edges.size(); k++) {
                if (edges[k].p != edges[0].p || !e.p->e[k].w.approximatelyEquals(e.p->e[0].w)) {
                    all_equal = false;
                    break;
                }

            }

            if (all_equal) {
                if (cached) {
                    for (int k = 1; k < edges.size(); k++) {
                        if (e.p->e[k].w != Complex::zero && e.p->e[k].w != Complex::one) {
                            cn.returnToCache(e.p->e[k].w);
                            return edges[0];
                        }
                    }
                    return edges[0];
                }
            }

            e = normalize(e, cached);

            assert(e.p->v == var || e.isTerminal());

            // look it up in the unique tables
            auto l = uniqueTable.lookup(e, false);

            assert(l.p->v == var || l.isTerminal());

            return l;
        }


        ///
        /// Compute table definitions
        ///
    public:
        void clearComputeTables() {

        }


    public:

        ComputeTable<mCachedEdge, mCachedEdge, mCachedEdge, Config::CT_VEC_ADD_NBUCKET> addTable{};
        ComputeTable2<mEdge, mEdge, mCachedEdge, Config::CT_MAT_MAT_MULT_NBUCKET> contTable{};

        key_2_new_key_node key_2_new_key_tree_header_element = {-1, -1, {}, nullptr};
        key_2_new_key_node *key_2_new_key_tree_header = &key_2_new_key_tree_header_element;

        key_2_new_key_node *append_new_key(key_2_new_key_node *self, float new_key) {

            auto it = self->next.find(new_key);

            if (it != self->next.end()) {
                return self->next[new_key];
            } else {
                self->next[new_key] = new key_2_new_key_node{static_cast<short>(self->level + 1), new_key, {}, self};
                return self->next[new_key];
            }
        }

        // todo : determine the use of this function
        template<class Edge>
        Edge T_add(const Edge &x, const Edge &y) {

            return y;
        }


        /**
         * Contract two TDDs, whose common indices will be contracted.
         *
         */
        TDD cont(TDD tdd1, TDD tdd2) {

            TDD res;

            std::vector<Index> var_out;
            std::vector<std::string> var_cont_temp;
            std::vector<std::string> var_cont;
            std::vector<std::string> var_out_key;

            // search for common indices
            int k;
            int k1;
            for (k = 0; k < tdd1.index_set.size(); ++k) {
                bool flag = true;

                for (k1 = 0; k1 < tdd2.index_set.size(); ++k1) {
                    if (tdd2.index_set[k1].idx == tdd1.index_set[k].idx &&
                        tdd2.index_set[k1].key == tdd1.index_set[k].key) {
                        var_cont_temp.push_back(tdd1.index_set[k].key);
                        flag = false;
                        break;
                    }
                }

                if (flag) {
                    var_out.push_back(tdd1.index_set[k]);
                    var_out_key.push_back(tdd1.index_set[k].key);
                }
            }

            for (k = 0; k < tdd2.index_set.size(); ++k) {
                bool flag = true;
                for (k1 = 0; k1 < tdd1.index_set.size(); ++k1) {
                    if (tdd1.index_set[k1].idx == tdd2.index_set[k].idx &&
                        tdd1.index_set[k1].key == tdd2.index_set[k].key) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    var_out.push_back(tdd2.index_set[k]);
                    var_out_key.push_back(tdd2.index_set[k].key);
                }
            }

            // avoid appending indices repeatedly
            for (k = 0; k < var_cont_temp.size(); ++k) {
                if (find(var_out_key.begin(), var_out_key.end(), var_cont_temp[k]) == var_out_key.end()) {
                    if (find(var_cont.begin(), var_cont.end(), var_cont_temp[k]) == var_cont.end()) {
                        var_cont.push_back(var_cont_temp[k]);
                    }
                }
            }


            // print information
            if (to_test) {
                std::cout << "TDD1: ";
                for (const auto &element: tdd1.key_2_index) {
                    std::cout << element << " ";
                }
                std::cout << std::endl;
                std::cout << "TDD2: ";
                for (const auto &element: tdd2.key_2_index) {
                    std::cout << element << " ";
                }
                std::cout << std::endl;

            }


            key_2_new_key_node *key_2_new_key1 = key_2_new_key_tree_header; //
            key_2_new_key_node *key_2_new_key2 = key_2_new_key_tree_header; //

            std::vector<std::string> new_key_2_index; // map : int(key) -> string(index)

            k1 = 0;
            int k2 = 0;
            int new_key = 0;
            int m1 = tdd1.key_2_index.size();
            int m2 = tdd2.key_2_index.size();
            int repeat_time = 1;
            float last_cont_idx = -2;

            // iterate through all the indices
            while (k1 < m1 || k2 < m2) {

                // When one of the TDDs is exhausted, append the remaining indices from the other TDD.
                if (k1 == m1) {
                    for (k2; k2 < m2; ++k2) {
                        key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
                        new_key_2_index.push_back(tdd2.key_2_index[k2]);
                        new_key++;
                    }
                    break;
                }
                if (k2 == m2) {
                    for (k1; k1 < m1; ++k1) {
                        key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
                        new_key_2_index.push_back(tdd1.key_2_index[k1]);
                        new_key++;
                    }
                    break;
                }

                // Otherwise, compare the order of the indices and append the smaller one.
                if (varOrder[tdd1.key_2_index[k1]] < varOrder[tdd2.key_2_index[k2]]) {
                    key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
                    new_key_2_index.push_back(tdd1.key_2_index[k1]);
                    new_key++;
                    k1++;
                } else if (varOrder[tdd1.key_2_index[k1]] > varOrder[tdd2.key_2_index[k2]]) {
                    key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
                    new_key_2_index.push_back(tdd2.key_2_index[k2]);
                    new_key++;
                    k2++;
                } else if (find(var_out_key.begin(), var_out_key.end(), tdd1.key_2_index[k1]) == var_out_key.end()) {
                    // The index is going to be contracted
                    if (new_key - last_cont_idx <= 0.5) {
                        last_cont_idx = last_cont_idx + 1 / (3 * nqubits) * repeat_time;
                        repeat_time += 1;
                        key_2_new_key1 = append_new_key(key_2_new_key1, last_cont_idx);
                        key_2_new_key2 = append_new_key(key_2_new_key2, last_cont_idx);
                        k1++;
                        k2++;
                    } else {
                        key_2_new_key1 = append_new_key(key_2_new_key1, new_key - 0.5);
                        key_2_new_key2 = append_new_key(key_2_new_key2, new_key - 0.5);
                        last_cont_idx = new_key - 0.5;
                        repeat_time = 1;
                        k1++;
                        k2++;
                    }

                } else {
                    // index with the same order but not going to be contracted
                    key_2_new_key1 = append_new_key(key_2_new_key1, new_key);
                    key_2_new_key2 = append_new_key(key_2_new_key2, new_key);
                    new_key_2_index.push_back(tdd1.key_2_index[k1]);
                    new_key++;
                    k1++;
                    k2++;
                }
            }

            res.index_set = var_out;
            res.key_2_index = new_key_2_index;

            if (to_test) {
                std::cout << "TDD1: ";
                for (const auto &element: tdd1.key_2_index) {
                    std::cout << element << " ";
                }
                std::cout << std::endl;

                std::cout << "TDD2: ";
                for (const auto &element: tdd2.key_2_index) {
                    std::cout << element << " ";
                }
                std::cout << std::endl;
            }


            [[maybe_unused]] const auto before = cn.cacheCount();

            res.e = cont2(tdd1.e, tdd2.e, key_2_new_key1, key_2_new_key2, var_cont.size());

            if (to_test) {
                std::cout << "TDD: ";
                for (const auto &element: res.key_2_index) {
                    std::cout << element << " ";
                }
                std::cout << std::endl;

            }

            var_out.clear();
            var_cont_temp.clear();
            var_out_key.clear();
            var_cont.clear();

            if (!res.e.w.exactlyZero() && !res.e.w.exactlyOne()) {
                cn.returnToCache(res.e.w);
                res.e.w = cn.lookup(res.e.w);
            }

            [[maybe_unused]] const auto after = cn.cacheCount();
            assert(before == after);

            return res;
        }


    private:

        /**
         * Add two TDDs.
         * @param x root edge for TDD F
         * @param y root edge for TDD G
         * @return
         */
        template<class Node>
        Edge<Node> T_add2(const Edge<Node> &x, const Edge<Node> &y) {

            // guarantee the order
            if (x.p > y.p) {
                return T_add2(y, x);
            }

            // trivial cases
            if (x.p == nullptr) {
                return y;
            }
            if (y.p == nullptr) {
                return x;
            }

            if (x.w.exactlyZero()) {
                if (y.w.exactlyZero()) {
                    return Edge<Node>::zero;
                }
                auto r = y;
                r.w = cn.getCached(CTEntry::val(y.w.r), CTEntry::val(y.w.i));
                return r;
            }
            if (y.w.exactlyZero()) {
                auto r = x;
                r.w = cn.getCached(CTEntry::val(x.w.r), CTEntry::val(x.w.i));
                return r;
            }

            // rF = rG, return G with new weight F.w + G.w
            if (x.p == y.p) {
                auto r = y;
                r.w = cn.addCached(x.w, y.w);
                if (r.w.approximatelyZero()) {
                    cn.returnToCache(r.w);
                    return Edge<Node>::zero;
                }
                return r;
            }

            auto xCopy = x;
            auto yCopy = y;
            if (x.w != Complex::one) {
                // guarantee the uniqueness of look up
                xCopy.w = Complex::one;
                yCopy.w = cn.divCached(y.w, x.w);
            }


            auto r = addTable.lookup({xCopy.p, xCopy.w}, {yCopy.p, yCopy.w});

            if (r.p != nullptr) {
                if (r.w.approximatelyZero()) {
                    return Edge<Node>::zero;
                }
                auto c = cn.getCached(r.w);
                if (x.w != Complex::one) {
                    // multiply back
                    cn.mul(c, c, x.w);
                    cn.returnToCache(yCopy.w);
                }
                return {r.p, c};
            }

            // w is the smaller index between rF and rG
            const Qubit w = (x.isTerminal() || (!y.isTerminal() && y.p->v > x.p->v))
                            ? y.p->v
                            : x.p->v;

            int n = (x.p->v != w) ? y.p->e.size() : x.p->e.size();

            std::vector<Edge<Node>> edge(n);
            for (std::size_t i = 0U; i < n; i++) {
                // e1 <- rF_{x=i} if x.p->v == w else F
                // Note : We use xCopy instead of x, in order to guarantee the uniqueness in lookup table
                Edge<Node> e1{};
                if (!x.isTerminal() && x.p->v == w) {
                    e1 = x.p->e[i];
                    if (e1.w != Complex::zero) {
                        e1.w = cn.mulCached(e1.w, xCopy.w);
                    }
                } else {
                    e1 = xCopy;
                    if (y.p->e[i].p == nullptr) {
                        e1 = {nullptr, Complex::zero};
                    }
                }
                // e2 <- rG_{y=i} if y.p->v == w else G
                Edge<Node> e2{};
                if (!y.isTerminal() && y.p->v == w) {
                    e2 = y.p->e[i];

                    if (e2.w != Complex::zero) {
                        e2.w = cn.mulCached(e2.w, yCopy.w);
                    }
                } else {
                    e2 = yCopy;
                    if (x.p->e[i].p == nullptr) {
                        e2 = {nullptr, Complex::zero};
                    }
                }

                // recurse here
                edge[i] = T_add2(e1, e2);


                if (!x.isTerminal() && x.p->v == w && e1.w != Complex::zero) {
                    cn.returnToCache(e1.w);
                }

                if (!y.isTerminal() && y.p->v == w && e2.w != Complex::zero) {
                    cn.returnToCache(e2.w);
                }
            }

            auto e = makeDDNode(w, edge, true);

            addTable.insert({xCopy.p, xCopy.w}, {yCopy.p, yCopy.w}, {e.p, e.w});
            if (x.w != Complex::one) {
                // multiply back
                cn.mul(e.w, e.w, x.w);
                cn.returnToCache(yCopy.w);
            }
            return e;
        }

        /**
         * Contract two edges.
         *
         * @note If one of the edges are empty, the function will return an empty edge.
         * @param x root edge for TDD F
         * @param y root edge for TDD G
         * @param key_2_new_key1 new index set for F
         * @param key_2_new_key2 new index set for G
         * @param var_num The size of the index set to be contracted
         * @return
         */
        Edge<mNode> cont2(const Edge<mNode> &x, const Edge<mNode> &y, key_2_new_key_node *key_2_new_key1,
                          key_2_new_key_node *key_2_new_key2, const int var_num) {


            using ResultEdge = Edge<mNode>;

            // Empty edge
            if (x.p == nullptr) {
                return {nullptr, Complex::zero};
            }
            if (y.p == nullptr) {
                return y;
            }

            // Contraction between two zero tensors always results in zero tensor
            if (x.w.exactlyZero() || y.w.exactlyZero()) {
                return ResultEdge::zero;
            }

            // Both edges are trivial, return a terminal node, whose weight is wx * wy * 2^len(var)
            if (x.p->v == -1 && y.p->v == -1) {
                auto c = cn.mulCached(x.w, y.w);

                if (var_num > 0) {
                    ComplexNumbers::mul(c, c, cn.getTemporary(pow(2, var_num), 0));
                }

                return ResultEdge::terminal(c);
            }

            key_2_new_key_node *temp_key_2_new_key2 = key_2_new_key2;
            while (temp_key_2_new_key2->level > y.p->v) {
                temp_key_2_new_key2 = temp_key_2_new_key2->father;
            }

            // x is trivial, y is not
            if (x.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key2->new_key - y.p->v) < 1e-10) {
                return ResultEdge{y.p, cn.mulCached(x.w, y.w)};
            }

            key_2_new_key_node *temp_key_2_new_key1 = key_2_new_key1;
            while (temp_key_2_new_key1->level > x.p->v) {
                temp_key_2_new_key1 = temp_key_2_new_key1->father;
            }

            // y is trivial, x is not
            if (y.p->v == -1 && var_num == 0 && std::abs(temp_key_2_new_key1->new_key - x.p->v) < 1e-10) {
                return ResultEdge{x.p, cn.mulCached(x.w, y.w)};
            }


            auto xCopy = x;
            xCopy.w = Complex::one;
            auto yCopy = y;
            yCopy.w = Complex::one;

            auto res = contTable.lookup(xCopy, yCopy, temp_key_2_new_key1, temp_key_2_new_key2);
            if (res.e.p != nullptr) {
                if (res.e.w.approximatelyZero()) {
                    return ResultEdge::zero;
                }
                auto e = ResultEdge{res.e.p, cn.getCached(res.e.w)};
                ComplexNumbers::mul(e.w, e.w, x.w);
                ComplexNumbers::mul(e.w, e.w, y.w);
                if (e.w.approximatelyZero()) {
                    cn.returnToCache(e.w);
                    return ResultEdge::zero;
                }
                if (res.cont_num != var_num) {
                    ComplexNumbers::mul(e.w, e.w,
                                        cn.getTemporary(pow(2, var_num - res.cont_num), 0));//对于一般形状的tensor,以2为底数可能有问题
                }
                return e;
            }


            float newk1 = temp_key_2_new_key1->new_key;
            float newk2 = temp_key_2_new_key2->new_key;

            ResultEdge e1{}, e2{}, r{};

            // recursively compute contraction
            if (newk1 > newk2) {
                if (std::abs(int(newk1 * 2)) % 2 == 1) {
                    // newk1 is going to be contracted
                    r = ResultEdge::zero;
                    ResultEdge etemp{};
                    // r <- \sum_k cont(F_{x=k}, G, var/{x})
                    for (int k = 0; k < x.p->e.size(); ++k) {
                        e1 = x.p->e[k];
                        e2 = yCopy;
                        etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
                        if (!etemp.w.exactlyZero()) {
                            if (r != ResultEdge::zero) {
                                auto temp = r.w;
                                r = T_add2(r, etemp);
                                cn.returnToCache(temp);
                                cn.returnToCache(etemp.w);
                            } else {
                                r = etemp;
                            }
                        }
                    }
                } else {
                    // newk1 is not going to be contracted, compute outer product (create a new node)
                    // qubit case: v.low <- L, v.high <- H
                    std::vector<ResultEdge> e;
                    for (int k = 0; k < x.p->e.size(); ++k) {
                        e1 = x.p->e[k];
                        e2 = yCopy;
                        e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
                    }
                    r = makeDDNode(Qubit(newk1), e, true);
                }
            } else if (newk1 < newk2) {
                if (std::abs(int(newk2 * 2)) % 2 == 1) {
                    r = ResultEdge::zero;
                    ResultEdge etemp{};
                    for (int k = 0; k < y.p->e.size(); ++k) {
                        e1 = xCopy;
                        e2 = y.p->e[k];
                        etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
                        if (!etemp.w.exactlyZero()) {
                            if (r != ResultEdge::zero) {
                                auto temp = r.w;
                                r = T_add2(r, etemp);
                                cn.returnToCache(temp);
                                cn.returnToCache(etemp.w);
                            } else {
                                r = etemp;
                            }
                        }
                    }
                } else {
                    std::vector<ResultEdge> e;
                    for (int k = 0; k < y.p->e.size(); ++k) {
                        e1 = xCopy;
                        e2 = y.p->e[k];
                        e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
                    }
                    r = makeDDNode(Qubit(newk2), e, true);
                }

            } else { // newk1 == newk2
                if (std::abs(int(newk2 * 2)) % 2 == 1) {
                    r = ResultEdge::zero;
                    ResultEdge etemp{};
                    for (int k = 0; k < x.p->e.size(); ++k) {
                        e1 = x.p->e[k];
                        e2 = y.p->e[k];
                        etemp = cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num - 1);
                        if (!etemp.w.exactlyZero()) {
                            if (r != ResultEdge::zero) {
                                auto temp = r.w;
                                r = T_add2(r, etemp);
                                cn.returnToCache(temp);
                                cn.returnToCache(etemp.w);
                            } else {
                                r = etemp;
                            }
                        }
                    }
                } else {
                    std::vector<ResultEdge> e;
                    for (int k = 0; k < x.p->e.size(); ++k) {
                        e1 = x.p->e[k];
                        e2 = y.p->e[k];
                        e.push_back(cont2(e1, e2, temp_key_2_new_key1, temp_key_2_new_key2, var_num));
                    }
                    r = makeDDNode(Qubit(newk1), e, true);
                }
            }

            contTable.insert(xCopy, yCopy, {r.p, r.w}, temp_key_2_new_key1, temp_key_2_new_key2, var_num);

            // r.w <- r.w * x.w * y.w
            if (!r.w.exactlyZero() && (x.w.exactlyOne() || !y.w.exactlyZero())) {
                if (r.w.exactlyOne()) {
                    r.w = cn.mulCached(x.w, y.w);
                } else {
                    ComplexNumbers::mul(r.w, r.w, x.w);
                    ComplexNumbers::mul(r.w, r.w, y.w);
                }
                if (r.w.approximatelyZero()) {
                    cn.returnToCache(r.w);
                    return ResultEdge::zero;
                }
            }

            return r;

        }


    public:
        ///
        /// Decision diagram size
        ///
        template<class Edge>
        unsigned int size(const Edge &e) {
            static constexpr unsigned int NODECOUNT_BUCKETS = 200000;
            static std::unordered_set<decltype(e.p)> visited{NODECOUNT_BUCKETS}; // 2e6
            visited.max_load_factor(10);
            visited.clear();
            return nodeCount(e, visited);
        }

    private:
        template<class Edge>
        unsigned int nodeCount(const Edge &e,
                               std::unordered_set<decltype(e.p)> &v) const {
            v.insert(e.p);
            unsigned int sum = 1;
            if (!e.isTerminal()) {
                for (const auto &edge: e.p->e) {
                    if (edge.p != nullptr && !v.count(edge.p)) {
                        sum += nodeCount(edge, v);
                    }
                }
            }
            return sum;
        }


        ///
        /// Printing and Statistics
        ///
    public:
        // print information on package and its members
        static void printInformation() {
            std::cout << "\n  compiled: " << __DATE__ << " " << __TIME__
                      << "\n  Complex size: " << sizeof(Complex) << " bytes (aligned "
                      << alignof(Complex) << " bytes)"
                      << "\n  ComplexValue size: " << sizeof(ComplexValue)
                      << " bytes (aligned " << alignof(ComplexValue) << " bytes)"
                      << "\n  ComplexNumbers size: " << sizeof(ComplexNumbers)
                      << " bytes (aligned " << alignof(ComplexNumbers) << " bytes)"
                      << "\n  mEdge size: " << sizeof(mEdge) << " bytes (aligned "
                      << alignof(mEdge) << " bytes)"
                      << "\n  mNode size: " << sizeof(mNode) << " bytes (aligned "
                      << alignof(mNode) << " bytes)"
                      << "\n  Package size: " << sizeof(Package) << " bytes (aligned "
                      << alignof(Package) << " bytes)"
                      << "\n"
                      << std::flush;
        }

        // print unique and compute table statistics

        void statistics() {
            std::cout << "DD statistics:\n";
            std::cout << "[UniqueTable] ";
            nodeUniqueTable.printStatistics();
            std::cout << "[Add] ";
            addTable.printStatistics();
            std::cout << "[Cont] ";
            contTable.printStatistics();
            std::cout << "[ComplexTable] ";
            cn.complexTable.printStatistics();
        }

    };

} // namespace dd
