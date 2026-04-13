#pragma once

#include <vector>
#include <utility>
#include <cstdint>
#include <functional>
#include <cmath>
#include <type_traits>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <map>

namespace Interpolation
{
using f64 = double;
using f32 = float;
using i32 = int32_t;
using u32 = uint32_t;

using vector_d = std::vector<double>;
struct matrix_d {
public:
   matrix_d(size_t n, size_t m, double c = 0.) : _data(n * m, c), _rows(n), _cols(m)
   {
   }

   double &operator()(size_t i, size_t j)
   {
      return _data[j + i * _rows];
   }

   size_t rows()
   {
      return _rows;
   }
   size_t cols()
   {
      return _cols;
   }

private:
   vector_d _data;
   size_t _rows, _cols;
};

} // namespace Interpolation
