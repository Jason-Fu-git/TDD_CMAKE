
#include <iostream>
#include <string>
#include <ctime>
#include "QuantumComputation.hpp"

#include "dd/Export.hpp"
#include "Cir_import.h"
#include "libkahypar.h"

using namespace std;

int save_data();


// parameters: [output mode] (output path) [input files]
// [output mode]: 0: print to console, 1: print to file
int main(int argc, char *argv[]) {

    // parse arguments
    assert(argc > 2);

    string inputPath = "Benchmarks/";
    string dataPath;

    string mode = argv[1];
    bool debug = false;

    int s = 2;
    if (mode == "0") {
        std::cout << "Output mode: print to console" << std::endl;
    } else if (mode == "1") {
        std::cout << "Output mode: print to file" << std::endl;
        std::cout << "Benchmark data will be saved into path " << argv[2] << std::endl;
        s = 3;
        dataPath = argv[2];
        if (dataPath.back() != '/') {
            dataPath += "/";
        }
    } else if (mode == "2") {
        debug = true;
    } else {
        std::cout << "Invalid output mode" << std::endl;
        return 1;
    }

    if (debug) {
        // iterate through benchmark files
        for (int i = s; i < argc; i++) {
            std::string file_name = argv[i];
            int n = get_qubits_num(inputPath + file_name);

            // print the circuit
            auto gateSet = import_circuit(inputPath + file_name);
            ContractionTree::printQuantumCircuit(gateSet, n, true);
        }
    }

    std::map<OptimizingMethod, std::string> method_map = {
            {OptimizingMethod::PARTITION_1,  "PartitionScheme1"},
            {OptimizingMethod::PARTITION_2,  "PartitionScheme2"},
            {OptimizingMethod::GN_COMMUNITY, "GNCommunity"},
            {OptimizingMethod::GREEDY,       "GREEDY"},
            {OptimizingMethod::KAHYPAR,      "KaHyPar"}
    };

    // iterate through methods
    for (auto const &[method, method_name]: method_map) {
        std::cout << "Optimization method: " << method_name << std::endl;
        std::ofstream *ofile = nullptr;
        if (mode == "1") {
            ofile = new std::ofstream(dataPath + method_name + ".csv", ios::app);
            assert(ofile->is_open());
            *ofile << "Benchmark" << "," << "Gate num" << ","
                   << "Time" << "," << "Node num max" << "," << "Node num final" << "," << "Contraction cost (log)"
                   << std::endl;
        }
        std::cout << "===================================" << std::endl;
        // iterate through benchmark files
        for (int i = s; i < argc; i++) {
            std::string file_name = argv[i];
            int n = get_qubits_num(inputPath + file_name);
            auto dd = std::make_unique<dd::Package<>>(10 * n);
            int *nodes;

            std::cout << "File name:" << file_name << std::endl;

            nodes = Simulate_with_ContractionOptimizer(inputPath, file_name, dd, method, ofile, debug);

            std::cout << "Nodes max:" << *nodes << std::endl;
            std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
            std::cout << "===================================" << std::endl;

            delete[] nodes;
        }
        if (mode == "1") {
            ofile->close();
            delete ofile;
        }

    }


    return 0;
}
