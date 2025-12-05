
#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/to_string.hpp"

namespace Bmad {
std::string to_string(const Bmad::FortranCharArray1D& arr) {
  auto vec = arr.to_vector();
  if (vec.empty())
    return "[]";

  std::ostringstream oss;
  oss << "[";
  auto it = vec.begin();
  oss << "'" << *it << "'";
  for (++it; it != vec.end(); ++it) {
    oss << ", '" << *it << "'";
  }
  oss << "]";
  return oss.str();
}

} // namespace Bmad
