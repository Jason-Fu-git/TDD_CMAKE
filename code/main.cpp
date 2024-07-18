
#include <iostream>
#include <string>
#include <ctime>
#include "QuantumComputation.hpp"

#include "dd/Export.hpp"
#include "Cir_import.h"
#include "dd/Tensor.hpp"
#include "dd/Graph.hpp"

using namespace std;

int save_data();


int main(int argc, char *argv[]) {

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
//
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