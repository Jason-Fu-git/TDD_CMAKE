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
     * Graph class
     * @tparam T type of data stored in a node
     * @author Jason Fu
     */
    template<typename T>
    class Graph {
    public:
        // node structure in a graph
        struct GraphNode {
            T data;    // data stored in the node
            std::set<int> neighbors; // index in a GraphEdge vector
            // ====================================
            // variables useful in the GN algorithm
            long long edgeCount;
            int distance;

            // ====================================
            explicit GraphNode(T data) : data(data), edgeCount(0), distance(-1) {}

            // reset predecessor and edgeCount
            void reset() {
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
            void reset() {
                betweenness = 0;
                removed = false;
            }

            bool operator<(const GraphEdge &other) const {
                return betweenness < other.betweenness;
            }

            int neighbor(int self) const {
                return self == u ? v : u;
            }
        };


        // constructor
        explicit Graph(const std::vector<T> &data) {
            nodes.reserve(data.size());
            for (const auto &item: data) {
                nodes.emplace_back(item);
            }
        }

        /**
         * add an edge to the graph
         * @param u index for node u (begin in 0)
         * @param v index for node v (begin in 0)
         */
        void addEdge(int u, int v) {
            edges.emplace_back(u, v);
            nodes[u].neighbors.insert(edges.size() - 1);
            nodes[v].neighbors.insert(edges.size() - 1);
        }

        std::vector<GraphNode> nodes;
        std::vector<GraphEdge> edges;

        /**
         * start the GN algorithm O(m^2 n)
         * @return an order vector of edges by their removal order
         * The first element in the vector is the first edge being removed.
         */
        std::vector<GraphEdge> girvanNewman(){
            // start the GN algorithm
            for (auto &edge: edges) {
                edge.reset();
            }
            std::vector<GraphEdge> result;
            for (int i = 0; i < edges.size(); ++i) {
                // calculate all the edge betweenness
                calculateEdgeBetweenness();
                // find the edge with the highest betweenness
                GraphEdge &maxEdge = *std::max_element(edges.begin(), edges.end());
                // remove the edge
                maxEdge.removed = true;
                // push to the result
                result.push_back(maxEdge);
            }
            return result;
        }

    private:

        // calculate edge betweenness for each edge. O(mn)
        void calculateEdgeBetweenness(){
            for (auto &edge: edges) {
                edge.betweenness = 0;
            }
            for (int i = 0 ; i < nodes.size(); ++i) {
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
        void calculateEdgeBetweenness(int s){
            std::queue<int> q;
            nodes[s].distance = 0;
            nodes[s].edgeCount = 1;
            q.push(s);
            while (!q.empty()) {
                // pop the first element
                auto u = q.front();
                q.pop();
                // iterate through its neighbors
                for (auto &neighbor: nodes[u].neighbors) {
                    auto &edge = edges[neighbor];
                    if (!edge.removed) {
                        int v = edge.neighbor(u);
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
                    }else if (u.distance == v.distance - 1) {
                        // u -> v
                        edge.betweenness += (double) u.edgeCount / (double) v.edgeCount;
                    }
                }
            }
        }

    };
}; // namespace dd


#endif //TDD_C_GRAPH_HPP
