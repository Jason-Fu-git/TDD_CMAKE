#pragma once


#include <xtensor/xio.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include "Tdd.hpp"

namespace dd {

    /**
     * @brief Wraps the tensor's data and index set.
     * @var data - The data of the tensor. Using xtensor::xarray<ComplexValue> to store the data.
     * @var index_set - The index set for the tensor. Using std::vector<Index> to store the index set.
     *
     */
    struct Tensor {

        xt::xarray<ComplexValue> data;

        std::vector<Index> index_set;

        std::string name=NULL;

    };

    /**
     * @brief A collection of tensors.
     * @var tensors - The collection of tensors. Using std::vector<Tensor> to store the tensors.
     */
    struct TensorNetwork {
        std::vector<Tensor> tensors;
    };


}