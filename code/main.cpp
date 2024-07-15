
#include <iostream>
#include <string>
#include <ctime>
#include "QuantumComputation.hpp"

#include "dd/Export.hpp"
#include "Cir_import.h"
#include "dd/Tensor.hpp"

using namespace std;

int save_data();


#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>




int main(int argc, char *argv[]) {

//    qc::QuantumComputation qc1{};
//
//    string path = "Benchmarks/bv_10.qasm";
//    qc1.import(path, qc::Format::OpenQASM);
//    const qc::MatrixDD dd1 = buildFunctionality(&qc1, dd);
//    dd->printInformation();
//    serialize(dd1, "output.ser");
    assert(argc == 2);

    string path2 = "Benchmarks/";

    string file_name = argv[1];
    int *nodes;
    int n = get_qubits_num(path2 + file_name);
    auto dd = std::make_unique<dd::Package<>>(3 * n);
    std::cout << "File name:" << file_name << std::endl;
    nodes = Simulate_with_tdd(path2, file_name, dd);
    std::cout << "Naive Method" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;


    std::cout << "File name:" << file_name << std::endl;
    auto dd2 = std::make_unique<dd::Package<>>(3 * n);
    nodes = Simulate_with_partition1(path2, file_name, dd2);

    std::cout << "Partition 1" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    auto dd3 = std::make_unique<dd::Package<>>(3 * n);
    nodes = Simulate_with_partition2(path2, file_name, dd3);

    std::cout << "Partition 2" << std::endl;
    std::cout << "Nodes max:" << *nodes << std::endl;
    std::cout << "Nodes Final:" << *(nodes + 1) << std::endl;
    std::cout << "===================================" << std::endl;

    std::cout << "File name:" << file_name << std::endl;
    auto dd4 = std::make_unique<dd::Package<>>(3 * n);
    nodes = Simulate_with_ContractionOptimizer(path2, file_name, dd3, 0);

    std::cout << "Exhaustive Search" << std::endl;
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
        auto dd = std::make_unique<dd::Package<>>(3 * n);
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
        auto dd = std::make_unique<dd::Package<>>(3 * n);
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
        auto dd = std::make_unique<dd::Package<>>(3 * n);
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