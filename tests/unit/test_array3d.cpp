#include "core/array/array3d.hpp"
#include "test_assert.hpp"

#include <iostream>

int main() {
  try {
    dec3d::core::Array3D<double> field(2, 3, 4, -1.0);

    DEC3D_CHECK_EQ(field.extent_r(), static_cast<std::size_t>(2));
    DEC3D_CHECK_EQ(field.extent_theta(), static_cast<std::size_t>(3));
    DEC3D_CHECK_EQ(field.extent_phi(), static_cast<std::size_t>(4));
    DEC3D_CHECK_EQ(field.size(), static_cast<std::size_t>(24));
    DEC3D_CHECK(!field.empty());

    field(1, 2, 3) = 42.5;
    DEC3D_CHECK_EQ(field(1, 2, 3), 42.5);
    DEC3D_CHECK_EQ(field(0, 0, 0), -1.0);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
