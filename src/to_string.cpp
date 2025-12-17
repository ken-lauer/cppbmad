
#include <string>

#include "bmad/fortran_arrays.hpp"
#include "bmad/to_string.hpp"

namespace Bmad {

std::string to_string(const FCharArray1D& arr) {
  if (!arr.is_valid())
    return "[]";
  std::ostringstream oss;
  oss << "[";
  for (int i = 0; i < arr.size(); ++i) {
    if (i > 0)
      oss << ", ";
    // implicit conversion to std::string
    oss << "\"" << static_cast<std::string>(arr[i]) << "\"";
  }
  oss << "]";
  return oss.str();
}

std::string repr(
    const void* obj,
    const std::string& class_name,
    const std::initializer_list<std::pair<std::string, std::string>>& fields)

{
  std::ostringstream oss;
  oss << class_name << "(0x" << std::hex
      << reinterpret_cast<std::uintptr_t>(obj);

  if (fields.size()) {
    for (const auto& field : fields) {
      oss << ", ";
      oss << field.first << "=" << field.second;
    }
  }

  oss << ")";
  return oss.str();
}
} // namespace Bmad
