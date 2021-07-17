// ====================================================================
// This file is part of Polylogarithm.
//
// Polylogarithm is licenced under the MIT License.
// ====================================================================

#pragma once
#include <complex>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

namespace polylogarithm {
namespace test {

/**
 * Reads complex numbers and corresponding polylogarithms from a file.
 *
 * @param filename file name
 * @tparam T data type for real and imaginary parts
 *
 * @return vector of pairs, where the first element in the pair is the
 * complex number and the second element is the corresponding
 * polylogarithm.
 */
template <typename T>
std::vector<std::pair<std::complex<T>, std::complex<T>>>
read_from_file(const std::string& filename)
{
   using Cmpl_t = std::complex<T>;

   std::vector<std::pair<Cmpl_t, Cmpl_t>> data;
   std::string line;
   std::ifstream fstr(filename);

   while (std::getline(fstr, line)) {
      T re_z{}, im_z{}, re_li{}, im_li{};

      std::istringstream isstr(line);
      isstr >> re_z >> im_z >> re_li >> im_li;

      data.push_back(std::make_pair(Cmpl_t(re_z, im_z), Cmpl_t(re_li, im_li)));
   }

   return data;
}

/**
 * Reads real numbers and corresponding results from a file.
 *
 * @param filename file name
 * @tparam T data type for real numbers
 *
 * @return vector of pairs, where the first element in the pair is the
 * number and the second element is the corresponding result.
 */
template <typename T>
std::vector<std::pair<T, T>>
read_reals_from_file(const std::string& filename)
{
   std::vector<std::pair<T, T>> data;
   std::string line;
   std::ifstream fstr(filename);

   while (std::getline(fstr, line)) {
      T x{}, y{};

      std::istringstream isstr(line);
      isstr >> x >> y;

      data.push_back(std::make_pair(x, y));
   }

   return data;
}

} // namespace test
} // namespace polylogarithm
