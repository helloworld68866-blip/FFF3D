#pragma once

#include <sstream>
#include <stdexcept>
#include <string>

namespace dec3d::test {

[[noreturn]] inline void Fail(const char* expression, const char* file, int line, const std::string& detail = {}) {
  std::ostringstream message;
  message << file << ":" << line << " assertion failed: " << expression;

  if (!detail.empty()) {
    message << " (" << detail << ")";
  }

  throw std::runtime_error(message.str());
}

template <typename TLeft, typename TRight>
inline void FailEqual(
    const char* left_expr,
    const char* right_expr,
    const TLeft& left,
    const TRight& right,
    const char* file,
    int line) {
  std::ostringstream detail;
  detail << left_expr << " != " << right_expr << " [" << left << " vs " << right << "]";
  Fail("equality", file, line, detail.str());
}

}  // namespace dec3d::test

#define DEC3D_CHECK(expr) \
  do { \
    if (!(expr)) { \
      ::dec3d::test::Fail(#expr, __FILE__, __LINE__); \
    } \
  } while (false)

#define DEC3D_CHECK_EQ(left, right) \
  do { \
    const auto& dec3d_left_eval = (left); \
    const auto& dec3d_right_eval = (right); \
    if (!(dec3d_left_eval == dec3d_right_eval)) { \
      ::dec3d::test::FailEqual(#left, #right, dec3d_left_eval, dec3d_right_eval, __FILE__, __LINE__); \
    } \
  } while (false)
