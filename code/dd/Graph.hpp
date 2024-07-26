//
// Created by Jason Fu on 24-7-16.
//

#ifndef TDD_C_GRAPH_HPP
#define TDD_C_GRAPH_HPP

#include <set>
#include <vector>
#include <queue>

namespace dd {

    /**
     * Graph class for the Girvan-Newman algorithm
     * The datastructure is adjacency array, which should be created externally.
     * An undirected graph is represented by two directed edges.
     */
    class Graph {
    public:
        // node structure in a graph
        struct GraphNode {
            int offset;    // offset in the adjacency array
            // ====================================
            // variables useful in the GN algorithm
            int distance;
            long long edgeCount;

            // ====================================
            explicit GraphNode(int offset = -1) : offset(offset), edgeCount(0), distance(-1) {}

            // reset predecessor and edgeCount
            inline void reset() {
                edgeCount = 0;
                distance = -1;
            }
        };


        // edge structure in a graph
        struct GraphEdge {
            int u; // index in a GraphNode vector
            int v; // index in a GraphNode vector
            // ====================================
            // variables useful in the GN algorithm
            double betweenness;
            bool removed;

            // ====================================
            GraphEdge(int u, int v) : u(u), v(v), betweenness(0), removed(false) {}

            // reset betweenness and removed
            inline void reset() {
                betweenness = 0;
                removed = false;
            }

            bool operator<(const GraphEdge &other) const {
                return betweenness < other.betweenness;
            }

            inline int neighbor(int self) const {
                return self == u ? v : u;
            }
        };

        std::vector<GraphNode> nodes; // The last element should be the length of the adjacency list (aka `edges`)
        std::vector<GraphEdge> edges;


        // constructor
        explicit Graph(int size) {
            nodes.reserve(size + 1);
            edges.reserve(size);
        }

        /**
         * start the GN algorithm O(m^2 n)
         * @return an order vector of edges by their removal order
         * The first element in the vector is the first edge being removed.
         */
        std::vector<GraphEdge> girvanNewman() {
            // start the GN algorithm
            std::vector<GraphEdge> result;
            for (int i = 0; i < edges.size() / 2; ++i) {
                // calculate all the edge betweenness
                calculateEdgeBetweenness();
                // find the edge with the highest betweenness
                GraphEdge &maxEdge = *std::max_element(edges.begin(), edges.end());
                // remove the edge
                maxEdge.removed = true;
                // remove the corresponding edge in the adjacency list
                for (int j = nodes[maxEdge.v].offset; j < nodes[maxEdge.v + 1].offset; ++j) {
                    if (edges[j].v == maxEdge.u) {
                        edges[j].removed = true;
                        break;
                    }
                }
                // push to the result
                result.push_back(maxEdge);
            }
            return result;
        }

        /**
         * get the connected components of the graph
         * @return a vector of connected components
         */
        std::vector<Graph *> getConnectedComponents() {
            int *ccIndex = new int[nodes.size() - 1]; // index of the connected component for each node
            for (int i = 0; i < nodes.size() - 1; ++i) {
                ccIndex[i] = -1;
            }
            // calculate connected components
            // O(m)
            int index = 0;
            for (int j = 0; j < nodes.size() - 1; ++j) {
                if (ccIndex[j] == -1) {
                    std::queue<int> q;
                    q.push(j);
                    while (!q.empty()) {
                        int u = q.front();
                        q.pop();
                        ccIndex[u] = index;
                        for (int i = nodes[u].offset; i < nodes[u + 1].offset; ++i) {
                            int v = edges[i].v;
                            if (ccIndex[v] == -1) {
                                q.push(v);
                            }
                        }
                    }
                    ++index;
                }
            }
            // create connected components
            std::vector<Graph *> result(index);
            for (int i = 0; i < index; ++i) {
                // find the indices of the nodes in the connected component
                std::vector<int> nodeIndex;
                for (int j = 0; j < nodes.size() - 1; ++j) {
                    if (ccIndex[j] == i) {
                        nodeIndex.push_back(j);
                    }
                }
                // create a new graph
                result[i] = new Graph(nodeIndex.size());
                // create the nodes
                for (int j: nodeIndex) {
                    result[i]->nodes.push_back(nodes[j]);
                    result[i]->nodes.back().offset = result[i]->edges.size();
                    // copy the edges
                    for (int k = nodes[j].offset; k < nodes[j + 1].offset; ++k) {
                        auto &edge = edges[k];
                        result[i]->edges.push_back(edge);
                        result[i]->edges.back().u = edge.u;
                        result[i]->edges.back().v = edge.v;
                    }
                }
                // update the end offset
                result[i]->nodes.emplace_back(result[i]->edges.size());
            }
            delete[] ccIndex;
            return result;
        }

    private:

        // calculate edge betweenness for each edge. O(mn)
        void calculateEdgeBetweenness() {
            for (auto &edge: edges) {
                edge.betweenness = 0;
            }
            for (int i = 0; i < nodes.size() - 1; ++i) {
                for (auto &node: nodes) {
                    node.reset();
                }
                calculateEdgeBetweenness(i);
            }
        }

        /**
         * calculate edge betweenness from source s to all the connected nodes. O(m)
         * @param s index of the source node
         */
        void calculateEdgeBetweenness(int s) {
            std::queue<int> q;
            nodes[s].distance = 0;
            nodes[s].edgeCount = 1;
            q.push(s);
            while (!q.empty()) {
                // pop the first element
                auto u = q.front();
                q.pop();
                // iterate through its neighbors
                for (int i = nodes[u].offset; i < nodes[u + 1].offset; ++i) {
                    auto &edge = edges[i];
                    if (!edge.removed) {
                        int v = edge.v;
                        // undiscovered
                        if (nodes[v].distance == -1) {
                            nodes[v].distance = nodes[u].distance + 1;
                            nodes[v].edgeCount = nodes[u].edgeCount;
                            q.push(v);
                        }
                            // discovered
                        else if (nodes[v].distance == nodes[u].distance + 1) {
                            nodes[v].edgeCount += nodes[u].edgeCount;
                        }
                    }
                }
            }
            // calculate edge betweenness
            for (auto &edge: edges) {
                if (!edge.removed) {
                    auto &u = nodes[edge.u];
                    auto &v = nodes[edge.v];
                    if (u.distance == v.distance + 1) {
                        // v -> u
                        edge.betweenness += (double) v.edgeCount / (double) u.edgeCount;
                    } else if (u.distance == v.distance - 1) {
                        // u -> v
                        edge.betweenness += (double) u.edgeCount / (double) v.edgeCount;
                    }
                }
            }
        }

    };
}; // namespace dd


#endif //TDD_C_GRAPH_HPP
