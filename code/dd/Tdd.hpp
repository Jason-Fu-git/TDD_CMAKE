#pragma once


#include "Definitions.hpp"
#include "Edge.hpp"
#include "Node.hpp"

#include <map>

namespace dd {

    /**
     * @brief In and out index for a tensor.
     * @var key - String representation of the index.
     * @var idx - Numbered in short.
     */
    struct Index {
        std::string key;
        short idx;
    };

    inline bool operator==(const Index &lhs, const Index &rhs) {
        return lhs.idx == rhs.idx && lhs.key == rhs.key;
    }

    inline bool operator!=(const Index &lhs, const Index &rhs) {
        return !(lhs == rhs);
    }

    inline bool operator<(const Index &lhs, const Index &rhs) {
        return lhs.key < rhs.key || (lhs.key == rhs.key && lhs.idx < rhs.idx);
    }

    inline bool operator<=(const Index &lhs, const Index &rhs) {
        return lhs < rhs || lhs == rhs;
    }

    inline bool operator>(const Index &lhs, const Index &rhs) {
        return !(lhs <= rhs);
    }

    inline bool operator>=(const Index &lhs, const Index &rhs) {
        return !(lhs < rhs);
    }

    /**
     * @brief Tensor Decision Diagram
     * @var e - The root edge for the TDD.
     * @var index_set - The set of indices for the TDD.
     * @var key_2_index - The map from key to index. (int -> string)
     */
    struct TDD {

        Edge<mNode> e;
        std::vector<Index> index_set;
        std::vector<std::string> key_2_index;

    };


    struct key_2_new_key_node {
        short level;
        float new_key;
        std::map<float, key_2_new_key_node *> next;
        key_2_new_key_node *father;
    };

}