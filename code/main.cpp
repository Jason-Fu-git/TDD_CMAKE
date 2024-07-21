
#include <iostream>
#include <string>
#include <ctime>
#include "QuantumComputation.hpp"

#include "dd/Export.hpp"
#include "Cir_import.h"
#include "libkahypar.h"

using namespace std;

int save_data();


int main(int argc, char *argv[]) {


//    kahypar_context_t *context = kahypar_context_new();
//    kahypar_configure_context_from_file(context, "config/kahypar_config.ini");
//
//    kahypar_set_seed(context, 42);
//
//    const kahypar_hypernode_id_t num_vertices = 4;
//    const kahypar_hyperedge_id_t num_hyperedges = 2;
//
//    std::unique_ptr<kahypar_hyperedge_weight_t[]> hyperedge_weights = std::make_unique<kahypar_hyperedge_weight_t[]>(2);
//
//    // force the cut to contain hyperedge 0 and 2
//    hyperedge_weights[0] = 1;
//    hyperedge_weights[1] = 1;
//
//    std::unique_ptr<size_t[]> hyperedge_indices = std::make_unique<size_t[]>(3);
//
//    hyperedge_indices[0] = 0;
//    hyperedge_indices[1] = 2;
//    hyperedge_indices[2] = 4;
//
//    std::unique_ptr<kahypar_hyperedge_id_t[]> hyperedges = std::make_unique<kahypar_hyperedge_id_t[]>(4);
//
//    // hypergraph from hMetis manual page 14
//    hyperedges[0] = 0;
//    hyperedges[1] = 1;
//    hyperedges[2] = 2;
//    hyperedges[3] = 3;
//
//    const double imbalance = 0.03;
//    const kahypar_partition_id_t k = 2;
//
//    kahypar_hyperedge_weight_t objective = 0;
//
//    std::vector<kahypar_partition_id_t> partition(num_vertices, -1);
//
//    kahypar_partition(num_vertices, num_hyperedges,
//                      imbalance, k,
//            /*vertex_weights */ nullptr, hyperedge_weights.get(),
//                      hyperedge_indices.get(), hyperedges.get(),
//                      &objective, context, partition.data());
//
//    for (int i = 0; i != num_vertices; ++i) {
//        std::cout << i << ":" << partition[i] << std::endl;
//    }
//
//    kahypar_context_free(context);


//    qc::QuantumComputation qc1{};
//
//    string path = "Benchmarks/bv_10.qasm";
//    qc1.import(path, qc::Format::OpenQASM);
//    const qc::MatrixDD dd1 = buildFunctionality(&qc1, dd);
//    dd->printInformation();
//    serialize(dd1, "output.ser");

//    dd::Graph<int> g(v);
//    g.addEdge(1, 2);
//    g.addEdge(2, 3);
//    g.addEdge(1, 3);
//
//    g.addEdge(3, 4);
//
//    g.addEdge(4, 5);
//    g.addEdge(5, 6);
//    g.addEdge(6, 4);
//
//    g.addEdge(2, 6);
//    g.addEdge(1, 5);
//
//    auto order = g.girvanNewman();
//    for (int i = 0; i < order.size(); ++i) {
//        printf("Edge %d -> %d\n", order[i].u, order[i].v);
//    }
//    std::vector<int> v = {0, 1, 2, 3, 4, 5, 6};


    assert(argc == 2);

    string path2 = "Benchmarks/";

    string file_name = argv[1];
    int *nodes;
    int n = get_qubits_num(path2 + file_name);
    auto dd = std::make_unique<dd::Package<>>(4 * n);
    std::cout << "File name:" << file_name << std::endl;
    nodes = Simulate_with_tdd(path2, file_name, dd);
    std::cout << "Naive Method" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;


    std::cout << "File name:" << file_name << std::endl;
    auto dd2 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_partition1(path2, file_name, dd2);

    std::cout << "Partition 1" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    auto dd3 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_partition2(path2, file_name, dd3);

    std::cout << "Partition 2" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

//    std::cout << "File name:" << file_name << std::endl;
//    std::cout << "Exhaustive Search" << std::endl;
//
//    auto dd4 = std::make_unique<dd::Package<>>(4 * n);
//    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd4, OptimizingMethod::EXHAUSTIVE_SEARCH);
//
//
//    std::cout << "Nodes max:" << *nodes << std::endl;
//    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
//    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    std::cout << "Partition Scheme 1" << std::endl;
    auto dd5 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd5, OptimizingMethod::PARTITION_1);

    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    std::cout << "Partition Scheme 2" << std::endl;

    auto dd6 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd6, OptimizingMethod::PARTITION_2);
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    std::cout << "GN Community" << std::endl;

    auto dd7 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd7, OptimizingMethod::GN_COMMUNITY);
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    std::cout << "GREEDY" << std::endl;

    auto dd8 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd8, OptimizingMethod::GREEDY);
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    auto dd9 = std::make_unique<dd::Package<>>(4 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd9, OptimizingMethod::KAHYPAR);
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

//    save_data();
//    system("pause");
    return 0;
}

int save_data() {

    std::ofstream ofile;

    string path2 = "Benchmarks/";
    std::string file_list_txt = "test.txt";
    std::ifstream file_list2;
    std::string line2;
    clock_t start2, finish2;
    double time2;
    int *nodes2;


    ofile.open("data.csv", ios::app);
    ofile << "Simulate_with_tdd" << endl;
    ofile << "benchmarks" << "," << "time" << "," << "node max" << "," << "node final" << endl;
    file_list2.open(file_list_txt);
    while (std::getline(file_list2, line2)) {
        std::cout << "file name:" << line2 << std::endl;
        int n = get_qubits_num(path2 + line2);
        auto dd = std::make_unique<dd::Package<>>(4 * n);
        start2 = clock();
        nodes2 = Simulate_with_tdd(path2, line2, dd);
        finish2 = clock();
        time2 = (double) (finish2 - start2) / CLOCKS_PER_SEC;
        std::cout << "time:" << time2 << std::endl;
        std::cout << "nodes max:" << *nodes2 << std::endl;
        std::cout << "nodes final:" << *(nodes2 + 1) << std::endl;
        ofile << line2 << "," << time2 << "," << *nodes2 << "," << *(nodes2 + 1) << endl;
    }
    file_list2.close();
    ofile.close();


    ofile.open("data.csv", ios::app);
    ofile << "Simulate_with_partition1" << endl;
    ofile << "benchmarks" << "," << "time" << "," << "node max" << "," << "node final" << endl;
    file_list2.open(file_list_txt);
    while (std::getline(file_list2, line2)) {
        std::cout << "file name:" << line2 << std::endl;
        int n = get_qubits_num(path2 + line2);
        auto dd = std::make_unique<dd::Package<>>(4 * n);
        start2 = clock();
        nodes2 = Simulate_with_partition1(path2, line2, dd);
        finish2 = clock();
        time2 = (double) (finish2 - start2) / CLOCKS_PER_SEC;
        std::cout << "time:" << time2 << std::endl;
        std::cout << "nodes max:" << *nodes2 << std::endl;
        std::cout << "nodes final:" << *(nodes2 + 1) << std::endl;
        ofile << line2 << "," << time2 << "," << *nodes2 << "," << *(nodes2 + 1) << endl;
    }
    file_list2.close();
    ofile.close();

    ofile.open("data.csv", ios::app);
    ofile << "Simulate_with_partition2" << endl;
    ofile << "benchmarks" << "," << "time" << "," << "node max" << "," << "node final" << endl;
    file_list2.open(file_list_txt);
    while (std::getline(file_list2, line2)) {
        std::cout << "file name:" << line2 << std::endl;
        int n = get_qubits_num(path2 + line2);
        auto dd = std::make_unique<dd::Package<>>(4 * n);
        start2 = clock();
        nodes2 = Simulate_with_partition2(path2, line2, dd);
        finish2 = clock();
        time2 = (double) (finish2 - start2) / CLOCKS_PER_SEC;
        std::cout << "time:" << time2 << std::endl;
        std::cout << "nodes max:" << *nodes2 << std::endl;
        std::cout << "nodes final:" << *(nodes2 + 1) << std::endl;
        ofile << line2 << "," << time2 << "," << *nodes2 << "," << *(nodes2 + 1) << endl;
    }
    file_list2.close();
    ofile.close();

    system("pause");
    return 0;
}