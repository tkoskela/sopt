#include <algorithm>
#include <exception>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <ctime>
#include <catch.hpp>

#include <sopt/imaging_forward_backward.h>
#include <sopt/logging.h>
#include <sopt/maths.h>
#include <sopt/relative_variation.h>
#include <sopt/sampling.h>
#include <sopt/types.h>
#include <sopt/utilities.h>
#include <sopt/wavelets.h>
#include <cppflow/cppflow.h>
#include "cppflow/ops.h"
#include "cppflow/model.h"

// This header is not part of the installed sopt interface
// It is only present in tests
#include <tools_for_tests/directories.h>
#include <tools_for_tests/tiffwrappers.h>

// \min_{x} ||\Psi^Tx||_1 \quad \mbox{s.t.} \quad ||y - Ax||_2 < \epsilon and x \geq 0

typedef double Scalar;
typedef sopt::Vector<Scalar> Vector;
typedef sopt::Matrix<Scalar> Matrix;
typedef sopt::Image<Scalar> Image;

TEST_CASE("Cppflow"){

  std::string const input_image = "cameraman256";
  Image const image = sopt::notinstalled::read_standard_tiff(input_image);

  //TODO use image.shape nad image.size instead of hardcoding values
  //TODO check data type and allow float or double
  
  // create a vector
  std::cout << "============Create vector" << std::endl;
  int const image_rows = image.rows();
  int const image_cols = image.cols();
  int const image_size = image.size();

  std::vector<int64_t> tensor_shape = {image_rows, image_cols};
  std::vector<float> values(image_size, 1);


  // Initialize all elements to image values.
  for (int i = 0; i < image.rows(); ++i) {
    for (int j = 0; j < image.cols(); ++j) {
      values[i*image_rows+j] = image(i,j);
    }
  }
  
  // create a tensor from vector
  std::cout << "============Create tensor" << std::endl;
  cppflow::tensor cf_tensor(values, tensor_shape);
  
  auto input = cppflow::cast(cf_tensor, TF_UINT8, TF_FLOAT);
  // Add batch dimension at start
  input = cppflow::expand_dims(input, 0);
  // add extra spatial dimension at end??
  // cppflow::decode_jpeg results in a shape (256, 256, 1) so we assume this is needed
  input = cppflow::expand_dims(input, -1);

  // Read in model
  std::cout << "============Reading model file" << std::endl;
  cppflow::model model(std::string("/home/sarah/Projects/LEXCI/sopt/cppflow/examples/lexci_model/model"));

  // Run model on image
  std::cout << "============Run model on tensor" << std::endl;
  auto output = model({{"serving_default_input0:0", input}}, {"StatefulPartitionedCall:0"});

  // Get values from output
  auto results = output[0].get_data<float>();
  std::vector<double> doubleResults(results.begin(), results.end());

  //TODO this has the image size hardcoded but it doesn't like using image_rows (int declared above) or image.rows()
  Eigen::Map<Eigen::Array<double, 256, 256>> model_output(doubleResults.data());
  // Map transposes the image so we transpose it back
  // This only works on square images, can't modify shape if it is not square
  model_output.transposeInPlace();

  sopt::utilities::write_tiff(model_output, "./cameraman_output.tiff");

  // compare input image to cleaned output image
  // calculate mean squared error sum_i ( ( x_true(i) - x_est(i) ) **2 ) 
  // check this is less than the number of pixels * 0.01

  auto mse = (image - model_output).square().sum() / image.size();
  CAPTURE(mse);
  CHECK(mse < 0.01);

  

}
