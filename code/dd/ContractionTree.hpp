//
// Created by Jason Fu on 24-7-14.
//

#ifndef TDD_C_CONTRACTIONTREE_HPP
#define TDD_C_CONTRACTIONTREE_HPP

#include "string"
#include "iostream"
#include "map"
#include "set"
#include "vector"
#include "dd/Tdd.hpp"


struct gate {
    std::string name;
    short int qubits[2];

    // default constructor
    gate() : qubits{0, 0} {}

    // copy constructor
    gate(const gate &g) : name(g.name), qubits{g.qubits[0], g.qubits[1]} {}

    // copy assignment
    gate &operator=(const gate &g) {
        name = g.name;
        qubits[0] = g.qubits[0];
        qubits[1] = g.qubits[1];
        return *this;
    }
};

using Index = dd::Index;

/**
 * Contraction tree used to optimise the contraction order of a circuit
 * @see https://doi.org/10.22331/q-2021-03-15-410
 * @author Jason Fu
 *
 */
class ContractionTree {
public:
    struct Node {
        Node *parent;
        Node *lc;
        Node *rc;
        std::set<Index> indexes; // indexes the node contains
        double cost; // contraction cost for the subtree represented by the node (log2)
        int gate_idx; // non-negative value represents a leaf, -1 represents a non-leaf
        int tdd_idx; // index of the TDD node in the TDD vector

        // constructors
        explicit Node(std::set<Index> *idxSet, int _gate_idx = -1) :
                parent(nullptr), lc(nullptr), rc(nullptr), cost(0), gate_idx(_gate_idx), tdd_idx(-1) {
            if (idxSet != nullptr) {
                for (const auto &idx: *idxSet) {
                    indexes.insert(idx);
                }
            }
        }

        explicit Node(const std::vector<Index> &idxSet, int _gate_idx = -1) :
                parent(nullptr), lc(nullptr), rc(nullptr), cost(0), gate_idx(_gate_idx), tdd_idx(-1),
                indexes(idxSet.begin(), idxSet.end()) {}

        Node(const Node &n) : parent(n.parent), lc(n.lc), rc(n.rc), indexes(n.indexes), cost(n.cost),
                              gate_idx(n.gate_idx) {}


        Node &operator=(const Node &n) {
            assert(this != &n);
            parent = n.parent;
            lc = n.lc;
            rc = n.rc;
            indexes = n.indexes;
            cost = n.cost;
            gate_idx = n.gate_idx;
            tdd_idx = n.tdd_idx;
            return *this;
        }

        // destructor, in a recursive manner
        ~Node() {
            delete lc;
            delete rc;
        }

        /**
         * Whether the node is a leaf node
         *
         */
        bool isLeaf() {
            return gate_idx >= 0;
        }

        /**
         * get the size of the tensor the node represents
         * @returns log_{index_width}(size)
         */
        int getTensorSize() {
            return indexes.size();
        }

        /**
         * construct a parent node based on the index set of its left and right children
         * the cost is calculated at the same time
         */
        static Node *constructParent(Node *lc, Node *rc, int index_width) {
            if (lc == nullptr && rc == nullptr)
                printf("ERROR: both children are null\n");
            assert(lc || rc);
            if (!lc || rc == lc) {
                return rc;
            } else if (!rc) {
                return lc;
            } else {
                // construct the parent node
                Node *parent = new Node(nullptr);
                // construct the index set for the parent node
                int commonIdxNum = 0;
                for (const auto &index: lc->indexes) {
                    if (rc->indexes.find(index) == rc->indexes.end()) {
                        parent->indexes.insert(index);
                    } else {
                        // the index will get contracted
                        ++commonIdxNum;
                    }
                }
                for (const auto &index: rc->indexes) {
                    if (lc->indexes.find(index) == lc->indexes.end()) {
                        parent->indexes.insert(index);
                    }
                }
                // perform the cost calculation
                double cost = logSum(lc->cost, rc->cost, index_width);
                cost = logSum(cost, (double) (parent->indexes.size() + commonIdxNum), index_width);
                parent->cost = cost;
                // link the parent node to its children
                parent->lc = lc;
                parent->rc = rc;
                lc->parent = parent;
                rc->parent = parent;
                return parent;
            }
        }

    };

    explicit ContractionTree(int _index_width = 2) {
        root = nullptr;
        size = 0;
        index_width = _index_width;
    }

    // destructor, in a recursive manner
    ~ContractionTree() {
        delete root;
    }

    // getters

    Node *&getRoot() {
        return root;
    }

    int getSize() const {
        return size;
    }

    /**
     *  Get the total contraction cost
     *  @note The cost is calculated based on typical matrix representation, not on exact tdd operation count!
     */
    double getCost() {
        return root->cost;
    }

    // node construction and deletion

    /**
     * Construct a node using Node's default constructor
     */
    Node *constructNode(std::set<dd::Index> *idxSet, int _gate_idx = -1) {
        Node *node = new Node(idxSet, _gate_idx);
        ++size;
        return node;
    }

    Node *constructNode(const std::vector<Index> &idxSet, int _gate_idx = -1) {
        Node *node = new Node(idxSet, _gate_idx);
        ++size;
        return node;
    }

    /**
     * Construct a node using Node's copy constructor
     */
    Node *constructNode(const Node &n) {
        Node *node = new Node(n);
        ++size;
        return node;
    }

    /**
     * construct a parent node based on the index set of its left and right children
     * the cost is calculated at the same time
     */
    Node *constructParent(Node *lc, Node *rc) {
        auto p = Node::constructParent(lc, rc, index_width);
        if (p != lc && p != rc) ++size;
        return p;
    }

    /**
     * Delete a single node, with the assumption that it has no children
     * @param node
     */
    void deleteNode(Node *node) {
        assert(node);
        assert(node->lc == nullptr && node->rc == nullptr);
        delete node;
        --size;
    }

    /**
     * concat two contraction trees
     */
    static ContractionTree *concat(ContractionTree *t1, ContractionTree *t2, int index_width) {
        auto *tr = new ContractionTree(index_width);
        tr->size += (t1->size + t2->size);
        tr->root = tr->constructParent(t1->root, t2->root);
        return tr;
    }

    /**
     * calculate the contraction cost of two contraction trees
     */
    static double contractionCost(ContractionTree *t1, ContractionTree *t2, int index_width) {
        std::vector<Index> idxSet;
        int commonIdxNum = 0;
        // iterate through t1's indexes
        for (const auto &index: t1->getRoot()->indexes) {
            if (t2->getRoot()->indexes.find(index)
                == t2->getRoot()->indexes.end()) {
                idxSet.push_back(index);
            } else {
                // the index will get contracted
                ++commonIdxNum;
            }
        }
        // iterate through t2's indexes
        for (const auto &index: t2->getRoot()->indexes) {
            if (t1->getRoot()->indexes.find(index)
                == t1->getRoot()->indexes.end()) {
                idxSet.push_back(index);
            }
        }
        // perform the cost calculation
        double cost = logSum(t1->getCost(), t2->getCost(), index_width);
        cost = logSum(cost, (double) (idxSet.size() + commonIdxNum), index_width);
        return cost;
    }

    /**
     * calculate the contraction cost of two nodes
     */
    static double contractionCost(Node *lc, Node *rc, int index_width) {
        std::vector<Index> idxSet;
        int commonIdxNum = 0;
        // iterate through t1's indexes
        for (const auto &index: rc->indexes) {
            if (lc->indexes.find(index)
                == lc->indexes.end()) {
                idxSet.push_back(index);
            } else {
                // the index will get contracted
                ++commonIdxNum;
            }
        }
        // iterate through t2's indexes
        for (const auto &index: lc->indexes) {
            if (rc->indexes.find(index)
                == rc->indexes.end()) {
                idxSet.push_back(index);
            }
        }
        // perform the cost calculation
        double cost = logSum(rc->cost, lc->cost, index_width);
        cost = logSum(cost, (double) (idxSet.size() + commonIdxNum, index_width), index_width);
        return cost;
    }

    /**
     * solve function w^a + w^b = w^x
     * @return x
     */
    inline static double logSum(double a, double b, int w) {
        // ensure a < b
        if (a > b) {
            return logSum(b, a, w);
        }
        if (b > a + 10)
            return b;
        // w^x = w^b (1 + w^{a-b})
        return b + log(1 + pow(w, a - b)) * log(w);
    }

    /**
     * solve function w^a - w^b = w^x
     * @return x
     */
    inline static double logSub(double a, double b, int w) {
        assert(a >= b);
        if (a - b > 10)
            return a;
        return a + log(1 - pow(w, b - a)) * log(w);
    }

    /**
     * Print the quantum circuit in a human-readable format
     * @param gateSet
     * @param qubitNum
     * @param printIndex
     * @author Jason Fu
     *
     */
    static void printQuantumCircuit(const std::map<int, gate> &gateSet, int qubitNum, bool printIndex) {
        gate placeholder;
        int maxNameWidth = 0;
        int gateIdx = 0;
        int maxVectorLength = 0;
        // reorganize the gate set by qubit
        std::vector<std::vector<gate>> qubitGateSet;
        qubitGateSet.reserve(qubitNum);
        for (int i = 0; i < qubitNum; ++i) {
            qubitGateSet.emplace_back();
        }

        for (auto &g: gateSet) {
            auto gate = g.second;

            if (gate.name.size() > maxNameWidth) {
                maxNameWidth = gate.name.size();
            }

            // only support cx gate as 2-qubit gate
            if (gate.name == "cx") {
                int control = gate.qubits[0];
                int target = gate.qubits[1];
                // pad the wires with placeholders
                for (int i = 0; i < qubitNum; ++i) {
                    while (qubitGateSet[i].size() < maxVectorLength) {
                        qubitGateSet[i].push_back(placeholder);
                    }
                }
                // push the cx gate
                qubitGateSet[control].push_back(gate);
                qubitGateSet[target].push_back(gate);
                qubitGateSet[control].back().name = "cx";
                qubitGateSet[target].back().name = printIndex ? "x " + std::to_string(gateIdx++) : "x";

                // pad the wires in between
                int lesser = control < target ? control : target;
                int greater = control > target ? control : target;
                for (int i = lesser + 1; i < greater; i++) {
                    qubitGateSet[i].push_back(placeholder);
                    qubitGateSet[i].back().name = "cr"; // cross
                }

                // update the max vector length
                if (qubitGateSet[control].size() > maxVectorLength) {
                    maxVectorLength = qubitGateSet[control].size();
                }
            } else {
                qubitGateSet[gate.qubits[0]].push_back(gate);
                if (printIndex) {
                    qubitGateSet[gate.qubits[0]].back().name = gate.name + " " + std::to_string(gateIdx++);
                }

                // update the max vector length
                if (qubitGateSet[gate.qubits[0]].size() > maxVectorLength) {
                    maxVectorLength = qubitGateSet[gate.qubits[0]].size();
                }
            }
        }

        // pad all the vectors
        for (int i = 0; i < qubitNum; ++i) {
            while (qubitGateSet[i].size() < maxVectorLength) {
                qubitGateSet[i].push_back(placeholder);
            }
        }

        if (printIndex) {
            maxNameWidth += 2 + (int) log10(gateIdx);
        }
        // print the circuit

        // gate figure
        std::string upperLeft = "┌";
        std::string upper = "─";
        std::string upperRight = "┐";
        std::string left = "┤";
        std::string right = "├";
        std::string lowerLeft = "└";
        std::string lower = "─";
        std::string lowerRight = "┘";

        // cx figure
        std::string control = "■";
        std::string cross = "┼";

        // wire
        std::string hwire = "─";
        std::string vwire = "│";

        // space
        std::string space = " ";

        // draw the circuit
        int width = maxNameWidth % 2 == 0 ? maxNameWidth + 3 : maxNameWidth + 2;
        int center = width / 2;
        int maxDisplay = 150 / width;
        int blockNum = maxVectorLength / maxDisplay + 1;
        int qubitLabelWidth = 2 + (int) log10(qubitNum) + 1;

        for (int blockIdx = 0; blockIdx < blockNum; blockIdx++) {
            for (int qubit = 0; qubit < qubitNum; qubit++) {
                // each wire occupies three lines
                for (int line = 0; line < 3; line++) {
                    // print arrow
                    if (blockIdx > 0)
                        std::cout << "«";
                    else
                        std::cout << space;

                    // print qubit label
                    int spaceCount = qubitLabelWidth + 1;
                    if (line == 1) {
                        std::cout << "q_" << std::to_string(qubit);
                        if (qubit > 0)
                            spaceCount -= 2 + (int) log10(qubit) + 1;
                        else
                            spaceCount -= 3;
                    }
                    for (int i = 0; i < spaceCount; i++) {
                        std::cout << space;
                    }

                    // print gate
                    for (gateIdx = blockIdx * maxDisplay;
                         gateIdx < blockIdx * maxDisplay + maxDisplay; gateIdx++) {
                        if (gateIdx >= maxVectorLength)
                            break;
                        auto &gate = qubitGateSet[qubit][gateIdx];
                        // placeholder
                        if (gate.name.empty()) {
                            if (line == 1) {
                                for (int i = 0; i < width; i++) {
                                    std::cout << hwire;
                                }
                            } else {
                                for (int i = 0; i < width; i++) {
                                    std::cout << space;
                                }
                            }
                        }
                            // cx (control side)
                        else if (gate.name == "cx") {
                            if (line == 0) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center && gate.qubits[0] > gate.qubits[1]) {
                                        std::cout << vwire;
                                    } else {
                                        std::cout << space;
                                    }
                                }
                            }

                            if (line == 1) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center) {
                                        std::cout << control;
                                    } else {
                                        std::cout << hwire;
                                    }
                                }
                            }

                            if (line == 2) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center && gate.qubits[0] < gate.qubits[1]) {
                                        std::cout << vwire;
                                    } else {
                                        std::cout << space;
                                    }
                                }
                            }
                        }
                            // cr (cross)
                        else if (gate.name == "cr") {
                            if (line == 0) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center) {
                                        std::cout << vwire;
                                    } else {
                                        std::cout << space;
                                    }
                                }
                            }

                            if (line == 1) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center) {
                                        std::cout << cross;
                                    } else {
                                        std::cout << hwire;
                                    }
                                }
                            }

                            if (line == 2) {
                                for (int i = 0; i < width; i++) {
                                    if (i == center) {
                                        std::cout << vwire;
                                    } else {
                                        std::cout << space;
                                    }
                                }
                            }
                        }

                            // single qubit gates
                        else {
                            if (line == 0) {
                                std::cout << upperLeft;
                                for (int i = 0; i < width - 2; i++) {
                                    std::cout << upper;
                                }
                                std::cout << upperRight;
                            }

                            if (line == 1) {
                                std::cout << left;
                                int spacePad = (width - 2 - gate.name.size()) / 2;
                                for (int i = 0; i < spacePad; i++) {
                                    std::cout << space;
                                }
                                std::cout << gate.name;
                                for (int i = 0; i < width - 2 - spacePad - gate.name.size(); i++) {
                                    std::cout << space;
                                }
                                std::cout << right;
                            }

                            if (line == 2) {
                                std::cout << lowerLeft;
                                for (int i = 0; i < width - 2; i++) {
                                    std::cout << lower;
                                }
                                std::cout << lowerRight;
                            }
                        }
                    }

                    // print the arrow
                    if (blockIdx < blockNum - 1)
                        std::cout << "»";
                    else
                        std::cout << space;

                    // print the end of the line
                    std::cout << std::endl;
                }
            }
            std::cout << std::endl;
        }

    }

private:
    Node *root;
    int size;
    int index_width; // width of the index, 2 for qubit
};


#endif //TDD_C_CONTRACTIONTREE_HPP
