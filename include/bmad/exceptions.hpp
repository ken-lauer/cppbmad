#pragma once

#include <stdexcept>

namespace Bmad {

class BmadException : public std::runtime_error {
public:
  explicit BmadException(const std::string &message)
      : std::runtime_error(message) {}
};

class InvalidIndexException : public BmadException {
public:
  InvalidIndexException(const std::string &index_type, int index, int max_value)
      : BmadException(
            "Invalid " + index_type + " index " + std::to_string(index) + " (valid range: 0-" +
            std::to_string(max_value - 1) + ")"
        ) {}
};

class NullPointerException : public BmadException {
public:
  NullPointerException(const std::string &context)
      : BmadException("Null pointer encountered in " + context) {}
};
}; // namespace Bmad
