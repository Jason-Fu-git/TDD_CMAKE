cmake_minimum_required(VERSION 2.8)
project(TDD_C CXX)

set(CMAKE_CXX_STANDARD 17)

# replace this with your actual library path where kahypar is installed
link_directories(/usr/local/lib)

add_subdirectory(code/deps/xtl)
add_subdirectory(code/deps/xtensor)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories(code)
include_directories(code/dd)
include_directories(code/operations)
include_directories(code/parsers/qasm_parser)

set(TDD_C_INCLUDES
        code/dd/Complex.hpp
        code/dd/ComplexCache.hpp
        code/dd/ComplexNumbers.hpp
        code/dd/ComplexTable.hpp
        code/dd/ComplexValue.hpp
        code/dd/ComputeTable.hpp
        code/dd/Control.hpp
        code/dd/DDpackageConfig.hpp
        code/dd/Definitions.hpp
        code/dd/Edge.hpp
        code/dd/Export.hpp
        code/dd/GateMatrixDefinitions.hpp
        code/dd/Node.hpp
        code/dd/Operations.hpp
        code/dd/Package.hpp
        code/dd/Package_fwd.hpp
        code/dd/Tdd.hpp
        code/dd/Tensor.hpp
        code/dd/UniqueTable.hpp
        code/operations/ClassicControlledOperation.hpp
        code/operations/CompoundOperation.hpp
        code/operations/Control.hpp
        code/operations/Expression.hpp
        code/operations/NonUnitaryOperation.hpp
        code/operations/Operation.hpp
        code/operations/OpType.hpp
        code/operations/StandardOperation.hpp
        code/operations/SymbolicOperation.hpp
        code/parsers/qasm_parser/Parser.hpp
        code/parsers/qasm_parser/Scanner.hpp
        code/parsers/qasm_parser/Token.hpp
        code/CircuitOptimizer.hpp
        code/Definitions.hpp
        code/Permutation.hpp
        code/QuantumComputation.hpp
        code/Cir_import.h
        code/dd/ContractionTree.hpp
        code/dd/ContractionOptimizer.hpp
        code/dd/Graph.hpp
)

set(TDD_C_SOURCES
        code/dd/Edge.cpp
        code/dd/Node.cpp
        code/dd/Operations.cpp
        code/operations/Expression.cpp
        code/operations/NonUnitaryOperation.cpp
        code/operations/Operation.cpp
        code/operations/StandardOperation.cpp
        code/operations/SymbolicOperation.cpp
        code/parsers/qasm_parser/Parser.cpp
        code/parsers/qasm_parser/Scanner.cpp
        code/parsers/GRCSParser.cpp
        code/parsers/QASMParser.cpp
        code/parsers/QCParser.cpp
        code/parsers/RealParser.cpp
        code/parsers/TFCParser.cpp
        code/CircuitOptimizer.cpp
        code/main.cpp
        code/QuantumComputation.cpp
        code/dd/ContractionOptimizer.cpp
)

add_executable(${PROJECT_NAME} ${TDD_C_SOURCES} ${TDD_C_INCLUDES})

target_link_libraries(${PROJECT_NAME} xtl)
target_link_libraries(${PROJECT_NAME} xtensor)
target_link_libraries(${PROJECT_NAME} kahypar)