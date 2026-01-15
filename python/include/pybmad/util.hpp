#pragma once

#include <functional>
#include <optional>
#include <string>

template <typename T>
std::optional<std::reference_wrapper<T>> make_opt_ref(std::optional<T> &opt_val) {
  if (opt_val.has_value()) {
    return std::ref(opt_val.value());
  }
  return std::nullopt;
}

// using json = nlohmann::json;
// template <typename T>
// std::string instance_to_json(const T& self) {
//   json j;
//   to_json(j, self);
//   return to_string(j);
// }
