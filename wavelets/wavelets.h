#ifndef SOPT_WAVELETS_WAVELETS_H
#define SOPT_WAVELETS_WAVELETS_H

#include <Eigen/Core>
#include "traits.h"

namespace sopt { namespace wavelets {

  //! Holds wavelets coefficients
  struct WaveletData {
    //! Type of the underlying scalar
    typedef t_real t_scalar;
    //! Type of the underlying vector
    typedef Eigen::Matrix<t_scalar, Eigen::Dynamic, 1> t_vector;
    //! Wavelet coefficient per-se
    const t_vector coefficients;

    //! Holds filters for direct transform
    struct {
      //! Low-pass filter for direct transform
      const t_vector low;
      //! High-pass filter for direct transform
      const t_vector high;
    } direct_filter;

    //! Holds filters for indirect transform
    struct {
      //! High-pass filter for direct transform
      const t_vector low_even;
      const t_vector low_odd;
      const t_vector high_even;
      const t_vector high_odd;
    } indirect_filter;

    //! Constructs from initializers
    WaveletData(std::initializer_list<t_scalar> const &coefs);
    //! Constructs from vector
    WaveletData(t_vector const &coefs);
  };
  extern const WaveletData Dirac;
  extern const WaveletData Daubechies1;
  extern const WaveletData Daubechies2;
  extern const WaveletData Daubechies3;
  extern const WaveletData Daubechies4;
  extern const WaveletData Daubechies5;
  extern const WaveletData Daubechies6;
  extern const WaveletData Daubechies7;
  extern const WaveletData Daubechies8;
  extern const WaveletData Daubechies9;
  extern const WaveletData Daubechies10;
}}
#endif
