#pragma once

#include "Complex.hpp"
#include "ComplexValue.hpp"
#include "Definitions.hpp"
#include "Edge.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace dd {


	// NOLINTNEXTLINE(readability-identifier-naming)
	/**
	 * \brief The node in TDD. It wraps its out-coming edges, succeeding nodes and a Qubit index.
	 * \var e - out-coming edges, organized in a vector
	 * \var next - used to link nodes in unique table
	 * \var ref - reference count
	 * \var v - variable index (nonterminal) value (-1 for terminal)
	 */
	struct mNode {
		//std::array<Edge<mNode>, RADIX> e{}; // edges out of this node
		std::vector<Edge<mNode>> e;          // edges out of this node
		mNode* next{};                      // used to link nodes in unique table
		RefCount ref{};                     // reference count
		Qubit v{}; // variable index (nonterminal) value (-1 for terminal)
		std::uint8_t flags = 0;

		// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
		static mNode terminal;

		static constexpr bool isTerminal(const mNode* p) { return p == &terminal; }
		static constexpr mNode* getTerminal() { return &terminal; }

	};
	using mEdge = Edge<mNode>;
	using mCachedEdge = CachedEdge<mNode>;

} // namespace dd


