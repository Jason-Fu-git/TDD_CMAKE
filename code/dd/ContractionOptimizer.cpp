//
// Created by Jason Fu on 24-7-14.
//

#include "ContractionOptimizer.hpp"
#include "Definitions.hpp"

// =============
// Base Class
// =============

dd::TDD
ContractionOptimizer::contractNode(Node *node, int &max_nodes, int &final_nodes, std::unique_ptr<dd::Package<>> &dd,
                                   bool release) {
    assert(node != nullptr);

    dd::TDD tdd;
    if (node->isLeaf()) {
        tdd = constructGate(gate_set->at(node->gate_idx).name, index_set->at(node->gate_idx), dd);
        if (release)
            dd->incRef(tdd.e);
    } else {
        if (node->lc == nullptr)
            return contractNode(node->rc, max_nodes, final_nodes, dd, release);

        if (node->rc == nullptr)
            return contractNode(node->lc, max_nodes, final_nodes, dd, release);

        auto lc = contractNode(node->lc, max_nodes, final_nodes, dd, release);
        auto rc = contractNode(node->rc, max_nodes, final_nodes, dd, release);
        tdd = dd->cont(lc, rc);
        if (release) {
            dd->decRef(lc.e);
            dd->decRef(rc.e);
            dd->incRef(tdd.e);
        }
    }
    // calculate node num
    final_nodes = dd->size(tdd.e);
    max_nodes = std::max(max_nodes, final_nodes);
    // collect garbage
    dd->garbageCollect();
    return tdd;
}

dd::TDD ContractionOptimizer::constructGate(std::string nam, const std::vector<dd::Index> &indexSet,
                                            std::unique_ptr<dd::Package<>> &dd) {
    std::map<std::string, int> gate_type;
    gate_type["x"] = 1;
    gate_type["y"] = 2;
    gate_type["z"] = 3;
    gate_type["h"] = 4;
    gate_type["s"] = 5;
    gate_type["sdg"] = 6;
    gate_type["t"] = 7;
    gate_type["tdg"] = 8;

    dd::TDD temp_tdd;

    if (nam == "cx") {
        temp_tdd = dd->cnot_2_TDD(indexSet, 1);
    } else {
        switch (gate_type[nam]) {
            case 1:
                temp_tdd = dd->Matrix2TDD(dd::Xmat, indexSet);
                break;
            case 2:
                temp_tdd = dd->Matrix2TDD(dd::Ymat, indexSet);
                break;
            case 3:
                temp_tdd = dd->diag_matrix_2_TDD(dd::Zmat, indexSet);
                break;
            case 4:
                temp_tdd = dd->Matrix2TDD(dd::Hmat, indexSet);
                break;
            case 5:
                temp_tdd = dd->diag_matrix_2_TDD(dd::Smat, indexSet);
                break;
            case 6:
                temp_tdd = dd->diag_matrix_2_TDD(dd::Sdagmat, indexSet);
                break;
            case 7:
                temp_tdd = dd->diag_matrix_2_TDD(dd::Tmat, indexSet);
                break;
            case 8:
                temp_tdd = dd->diag_matrix_2_TDD(dd::Tdagmat, indexSet);
                break;
            default:
                if (nam[0] == 'r' and nam[1] == 'z') {
                    std::regex pattern(R"(rz\((-?\d.\d+)\))");
                    std::smatch result;
                    regex_match(nam, result, pattern);
                    float theta = stof(result[1]);
                    temp_tdd = dd->diag_matrix_2_TDD(dd::Phasemat(theta), indexSet);
                    break;
                }
                if (nam[0] == 'u' and nam[1] == '1') {
                    std::regex para(".*?\\((.*?)\\)");
                    std::smatch result;
                    regex_match(nam, result, para);
                    float theta = match_a_string(result[1]);
                    temp_tdd = dd->diag_matrix_2_TDD(dd::Phasemat(theta), indexSet);
                    break;
                }
                if (nam[0] == 'u' and nam[1] == '3') {
                    std::regex para(".*?\\((.*?)\\)");
                    std::smatch result;
                    regex_match(nam, result, para);
                    std::vector<std::string> para2 = split(result[1], ",");
                    float theta = match_a_string(para2[0]);
                    float phi = match_a_string(para2[1]);
                    float lambda = match_a_string(para2[2]);
                    //dd::GateMatrix  U3mat = { { cos(theta / 2), 0 }, { -cos(lambda) * sin(theta / 2),-sin(lambda) * sin(theta / 2)} , { cos(phi) * sin(theta / 2),sin(phi) * sin(theta / 2) }, { cos(lambda + phi) * cos(theta / 2),sin(lambda + phi) * cos(theta / 2) }  };
                    temp_tdd = dd->Matrix2TDD(dd::U3mat(lambda, phi, theta), indexSet);
                    break;
                }
        }
    }
    return temp_tdd;
}

float ContractionOptimizer::match_a_string(std::string s) {
    std::smatch result;
    std::regex pattern("(-?\\d+.\\d+)");
    std::regex pattern2(R"((-?\d+.\d+)\*?pi/(\d+))");
    std::regex pattern3(R"((-?\d+.\d+)\*?pi)");
    std::regex pattern4("pi/(\\d+)");
    std::regex pattern5("(\\d+)");
    std::regex pattern6("-pi/(\\d+)");
    if (regex_match(s, result, pattern)) {
        return stof(result[1]);
    } else if (regex_match(s, result, pattern2)) {
        return stof(result[1]) * dd::PI / stof(result[2]);
    } else if (regex_match(s, result, pattern3)) {
        return stof(result[1]) * dd::PI;
    } else if (regex_match(s, result, pattern4)) {
        return dd::PI / stof(result[1]);
    } else if (regex_match(s, result, pattern5)) {
        return stof(result[1]);
    } else if (regex_match(s, result, pattern6)) {
        return -dd::PI / stof(result[1]);
    }
    std::cout << s << std::endl;
    std::cout << "Not Match" << std::endl;
    return 0.0;
}


std::vector<std::string> ContractionOptimizer::split(const std::string &s, const std::string &seperator) {
    std::vector<std::string> result;
    typedef std::string::size_type string_size;
    string_size i = 0;

    while (i != s.size()) {
        //找到字符串中首个不等于分隔符的字母；
        int flag = 0;
        while (i != s.size() && flag == 0) {
            flag = 1;
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[i] == seperator[x]) {
                    ++i;
                    flag = 0;
                    break;
                }
        }

        //找到又一个分隔符，将两个分隔符之间的字符串取出；
        flag = 0;
        string_size j = i;
        while (j != s.size() && flag == 0) {
            for (string_size x = 0; x < seperator.size(); ++x)
                if (s[j] == seperator[x]) {
                    flag = 1;
                    break;
                }
            if (flag == 0)
                ++j;
        }
        if (i != j) {
            result.push_back(s.substr(i, j - i));
            i = j;
        }
    }
    return result;
}

// =========
// ExhaustiveSearch
// =========

// O(2^{n^2})
ContractionTree *ExhaustiveSearchOptimizer::optimize() {
    std::vector<int> leaves;
    for (int i = 0; i < node_set.size(); i++) {
        leaves.push_back(i);
    }
    auto tr = findOptimalPartition(leaves, std::numeric_limits<unsigned long long>::max());
    printf("Exhaustive search performed %lld times\n", count);
    return tr;
}

ContractionTree *ExhaustiveSearchOptimizer::findOptimalPartition(std::vector<int> &nodeIdxs,
                                                                 unsigned long long minCost) {
    // trivial condition
    if (nodeIdxs.empty()) {
        return nullptr;
    }
    // first lookup in the dp table
    auto tr = searchDPTable(nodeIdxs);
    if (tr) {
        return tr;
    }
    // if not found, do exhaustive search
    tr = nullptr;
    if (nodeIdxs.size() == 1) {
        tr = new ContractionTree(index_width);
        auto node = tr->constructNode(*node_set[nodeIdxs[0]]);
        tr->getRoot() = node;
    } else if (nodeIdxs.size() == 2) {
        tr = new ContractionTree(index_width);
        auto lc = tr->constructNode(*node_set[nodeIdxs[0]]);
        auto rc = tr->constructNode(*node_set[nodeIdxs[1]]);
        auto p = tr->constructParent(lc, rc);
        tr->getRoot() = p;
    } else {
        // perform partition
        unsigned long long best_cost = minCost;
        ContractionTree *mLc = nullptr, *mRc = nullptr;
        // todo : not necessarily even partition
        Bitset mask((int) nodeIdxs.size() + 1);
        mask.set(0);
        while (true) {
            ++count;
            // generate partition
            std::vector<int> left, right;
            for (int i = 0; i < nodeIdxs.size(); i++) {
                if (mask.test(i)) {
                    left.push_back(nodeIdxs[i]);
                } else {
                    right.push_back(nodeIdxs[i]);
                }
            }
            if (right.empty())
                break;
            // compute cost
            auto lc = findOptimalPartition(left, best_cost);
            auto rc = findOptimalPartition(right, best_cost);
            if (lc && rc) {
                auto cost = ContractionTree::contractionCost(lc, rc, index_width);
                if (cost < best_cost) {
                    best_cost = cost;
                    mLc = lc;
                    mRc = rc;
                }
//                if (nodeIdxs.size() == gate_set->size()) {
//                    std::cout << "cost: " << cost << std::endl;
//                }
            }
            // next
            mask.add1();
        }
        if (mLc && mRc) {
            tr = ContractionTree::concat(mLc, mRc, index_width);
        }
    }
    // update dp table
    dumpDPTable(nodeIdxs, tr);

    return tr;
}

ContractionTree *ExhaustiveSearchOptimizer::searchDPTable(std::vector<int> &nodeIdxs) {
    // convert to string
    std::string key;
    for (auto i: nodeIdxs) {
        key += std::to_string(i) + ",";
    }

    // lookup
    if (dp_table.find(key) != dp_table.end()) {
        return dp_table[key];
    } else {
        return nullptr;
    }
}

void ExhaustiveSearchOptimizer::dumpDPTable(std::vector<int> &nodeIdxs, ContractionTree *result) {
    // convert to string
    std::string key;
    for (auto i: nodeIdxs) {
        key += std::to_string(i) + ",";
    }

    // insert
    dp_table[key] = result;
}

// ==============
// Partition Scheme 1
// ==============

// copied from Cir_import.h
ContractionTree *PartitionScheme1Optimizer::optimize() {
    auto gateSet = *gate_set;
    // perform partition
    // @var par : block index -> block gate map f, f(0) contains up part gates, f(1) contains down part gates.
    // Note all the gates are indices in the `gate_set`.
    std::map<int, std::map<int, std::vector<int>>> par;
    int cx_cut = 0;
    int block = 0;
    int gates_num = (int) gateSet.size();
    for (int k = 0; k < gates_num; k++) {
        std::string nam = gateSet[k].name;
        if (nam != "cx") {
            if (gateSet[k].qubits[0] <= qubits_num / 2) {
                // up part
                par[block][0].push_back(k);
            } else {
                // down part
                par[block][1].push_back(k);
            }
        } else {
            if (gateSet[k].qubits[0] <= qubits_num / 2 && gateSet[k].qubits[1] <= qubits_num / 2) {
                par[block][0].push_back(k);
            } else if (gateSet[k].qubits[0] > qubits_num / 2 && gateSet[k].qubits[1] > qubits_num / 2) {
                par[block][1].push_back(k);
            } else {
                // CNOT across the cut
                if (cx_cut <= max_cx_cut) {
                    if (gateSet[k].qubits[1] > qubits_num / 2) {
                        par[block][1].push_back(k);
                    } else {
                        par[block][0].push_back(k);
                    }
                    cx_cut += 1;
                } else {
                    block += 1;
                    cx_cut = 1;
                    if (gateSet[k].qubits[1] > qubits_num / 2) {
                        par[block][1].push_back(k);
                    } else {
                        par[block][0].push_back(k);
                    }
                }
            }
        }
    }
    // construct tree
    auto *tr = new ContractionTree(index_width);
    // iterate through each vertical cut
    Node *root = nullptr;
    for (auto &i: par) {
        auto block_gates = i.second;
        Node *vBlock_node = nullptr;
        // iterate through each horizontal cut
        for (auto &j: block_gates) {
            auto gates = j.second;
            // contract these gates
            Node *hBlock_node = nullptr;
            for (auto &g: gates) {
                auto gate = gateSet[g];
                auto indexSet = (*index_set)[g];
                auto node = tr->constructNode(indexSet, g);
                if (hBlock_node == nullptr) {
                    hBlock_node = node;
                } else {
                    hBlock_node = tr->constructParent(node, hBlock_node);
                }
            }
            // contract horizontal cuts
            if (vBlock_node == nullptr) {
                vBlock_node = hBlock_node;
            } else {
                vBlock_node = tr->constructParent(hBlock_node, vBlock_node);
            }
        }
        // contract vertical cuts
        if (root == nullptr) {
            root = vBlock_node;
        } else {
            root = tr->constructParent(vBlock_node, root);
        }
    }
    tr->getRoot() = root;
    return tr;
}

ContractionTree *PartitionScheme2Optimizer::optimize() {
    auto gateSet = *gate_set;
    // perform partitioning
    // @var par : block index -> block gate map f, f(0) contains up part gates, f(1) contains down part gates, and f(2)
    // contains the small block.
    // Note all the gates are indices in the `gate_set`.
    std::map<int, std::map<int, std::vector<int>>> par;
    int cx_cut = 0;
    int block = 0;
    int c_part_min = qubits_num / 2;
    int c_part_max = qubits_num / 2;
    for (int k = 0; k < gateSet.size(); k++) {
        std::string nam = gateSet[k].name;

        if (cx_cut <= cx_cut_max) {

            if (nam != "cx") {
                if (gateSet[k].qubits[0] <= qubits_num / 2) {
                    par[block][0].push_back(k);
                } else {
                    par[block][1].push_back(k);
                }
            } else {
                if (gateSet[k].qubits[0] <= qubits_num / 2 && gateSet[k].qubits[1] <= qubits_num / 2) {
                    par[block][0].push_back(k);
                } else if (gateSet[k].qubits[0] > qubits_num / 2 && gateSet[k].qubits[1] > qubits_num / 2) {
                    par[block][1].push_back(k);
                } else {
                    if (gateSet[k].qubits[1] > qubits_num / 2) {
                        par[block][1].push_back(k);
                    } else {
                        par[block][0].push_back(k);
                    }
                    cx_cut += 1;
                }
            }
        } else {
            if (nam != "cx") {
                if (gateSet[k].qubits[0] < c_part_min) {
                    par[block][0].push_back(k);
                } else if (gateSet[k].qubits[0] > c_part_max) {
                    par[block][1].push_back(k);
                } else {
                    par[block][2].push_back(k);
                }
            } else if (gateSet[k].qubits[0] >= c_part_min && gateSet[k].qubits[0] <= c_part_max &&
                       gateSet[k].qubits[1] >= c_part_min && gateSet[k].qubits[1] <= c_part_max) {
                par[block][2].push_back(k);
            } else if (gateSet[k].qubits[0] < c_part_min && gateSet[k].qubits[1] < c_part_min) {
                par[block][0].push_back(k);
            } else if (gateSet[k].qubits[0] > c_part_max && gateSet[k].qubits[1] > c_part_max) {
                par[block][1].push_back(k);
            } else {
                // CNOT gate across c_part_min or c_part_max
                int temp_c_min = std::min(c_part_min, (int) std::min(gateSet[k].qubits[0], gateSet[k].qubits[1]));
                int temp_c_max = std::max(c_part_max, (int) std::max(gateSet[k].qubits[0], gateSet[k].qubits[1]));
                if ((temp_c_max - temp_c_min) > c_part_width) {
                    block += 1;
                    cx_cut = 0;
                    c_part_min = qubits_num / 2;
                    c_part_max = qubits_num / 2;
                    if (gateSet[k].qubits[0] <= qubits_num / 2 && gateSet[k].qubits[1] <= qubits_num / 2) {
                        par[block][0].push_back(k);
                    } else if (gateSet[k].qubits[0] > qubits_num / 2 && gateSet[k].qubits[1] > qubits_num / 2) {
                        par[block][1].push_back(k);
                    } else {
                        if (gateSet[k].qubits[1] > qubits_num / 2) {
                            par[block][1].push_back(k);
                        } else {
                            par[block][0].push_back(k);
                        }
                        cx_cut += 1;
                    }
                } else {
                    par[block][2].push_back(k);
                    c_part_min = temp_c_min;
                    c_part_max = temp_c_max;
                }
            }
        }
    }
    // construct tree
    auto *tr = new ContractionTree(index_width);
    // iterate through each vertical cut
    Node *root = nullptr;
    for (auto &i: par) {
        auto block_gates = i.second;
        Node *vBlock_node = nullptr;
        // iterate through each horizontal cut
        for (auto &j: block_gates) {
            auto gates = j.second;
            // contract these gates
            Node *hBlock_node = nullptr;
            for (auto &g: gates) {
                auto gate = gateSet[g];
                auto indexSet = (*index_set)[g];
                auto node = tr->constructNode(indexSet, g);
                if (hBlock_node == nullptr) {
                    hBlock_node = node;
                } else {
                    hBlock_node = tr->constructParent(node, hBlock_node);
                }
            }
            // contract horizontal cuts
            if (vBlock_node == nullptr) {
                vBlock_node = hBlock_node;
            } else {
                vBlock_node = tr->constructParent(hBlock_node, vBlock_node);
            }
        }
        // contract vertical cuts
        if (root == nullptr) {
            root = vBlock_node;
        } else {
            root = tr->constructParent(vBlock_node, root);
        }
    }
    tr->getRoot() = root;
    return tr;
}

// ==================================================
// GNCommunityOptimizer
// ==================================================

ContractionTree *GNCommunityOptimizer::optimize() {
    // convert the circuit to a graph
    Graph *graph = constructGraph();
    // perform GN algorithm to find communities
    auto edgeOrder = graph->girvanNewman();
    // build a contraction tree based on the order
    auto *tr = buildTreeFromGraph(edgeOrder);
    delete graph;
    return tr;
}

// O(n^2), where n is the number of gates
Graph *GNCommunityOptimizer::constructGraph() {
    std::vector<int> data(gate_set->size());
    for (int i = 0; i < gate_set->size(); i++) {
        data.push_back(i);
    }
    auto *graph = new Graph(data);
    // iterate through each gate
    for (int i = 0; i < gate_set->size(); i++) {
        for (int j = i + 1; j < gate_set->size(); j++) {
            // find common indexes
            for (auto &index: (*index_set)[i]) {
                if (std::find((*index_set)[j].begin(), (*index_set)[j].end(), index) != (*index_set)[j].end()) {
                    graph->addEdge(i, j);
                    break;
                }
            }
        }
    }
    return graph;
}

ContractionTree *GNCommunityOptimizer::buildTreeFromGraph(const std::vector<GraphEdge> &edgeOrder) {
    auto *tr = new ContractionTree(index_width);
    // first construct gate nodes
    std::vector<Node *> nodes;
    std::vector<int> unionFindSet; // unionFind set
    for (int i = 0; i < gate_set->size(); i++) {
        auto *node = tr->constructNode(index_set->at(i), i);
        nodes.push_back(node);
        unionFindSet.push_back(i);
    }
    // then construct contraction tree
    for (int i = edgeOrder.size() - 1; i >= 0; i--) {
        auto &edge = edgeOrder[i];
        auto u = nodes[unionFind(unionFindSet, edge.u)];
        auto v = nodes[unionFind(unionFindSet, edge.v)];
        auto p = tr->constructParent(u, v);
        unionFindSet[edge.u] = edge.v;
        if (edge.u != edge.v) {
            nodes[edge.u] = nullptr;
        }
        nodes[unionFind(unionFindSet, edge.v)] = p;
        // the last node is the root
        if (i == 0) {
            // in the end, connect disjoint nodes to the root
            auto *prev = p;
            for (auto &node: nodes) {
                if (node != nullptr && node != prev) {
                    prev = tr->constructParent(prev, node);
                }
            }
            tr->getRoot() = prev;
        }
    }
    return tr;
}

int GNCommunityOptimizer::unionFind(const std::vector<int> &unionFindSet, int s) {
    while (unionFindSet[s] != s)
        s = unionFindSet[s];
    return s;
}
