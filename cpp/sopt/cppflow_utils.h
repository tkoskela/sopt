#ifndef SOPT_CPPFLOWUTILS_H
#define SOPT_CPPFLOWUTILS_H

#include <Eigen/Core>
#include "sopt/types.h"
#include <cppflow/cppflow.h>
#include "cppflow/ops.h"

namespace sopt {
namespace cppflowutils {
//! Converts a sopt::Image to a cppflow::tensor
cppflow::tensor convert_image_to_tensor(sopt::Image<double> const &image, int image_rows, int image_cols);

//! Convert a cppflow:tensor to an Eigen::Array
Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>> convert_tensor_to_image(std::vector<float> model_output, int image_rows, int image_cols);

}  // namespace cppflowutils
}  // namespace sopt
#endif
