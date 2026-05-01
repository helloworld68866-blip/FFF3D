#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace dec3d::core {

template <typename T>
class Array3D {
 public:
  Array3D() = default;

  Array3D(std::size_t radial, std::size_t polar, std::size_t azimuthal, const T& initial_value = T{})
      : radial_(radial),
        polar_(polar),
        azimuthal_(azimuthal),
        storage_(radial * polar * azimuthal, initial_value) {}

  [[nodiscard]] std::size_t extent_r() const noexcept { return radial_; }
  [[nodiscard]] std::size_t extent_theta() const noexcept { return polar_; }
  [[nodiscard]] std::size_t extent_phi() const noexcept { return azimuthal_; }
  [[nodiscard]] std::size_t size() const noexcept { return storage_.size(); }
  [[nodiscard]] bool empty() const noexcept { return storage_.empty(); }

  T& operator()(std::size_t radial, std::size_t polar, std::size_t azimuthal) {
    return storage_.at(linear_index(radial, polar, azimuthal));
  }

  const T& operator()(std::size_t radial, std::size_t polar, std::size_t azimuthal) const {
    return storage_.at(linear_index(radial, polar, azimuthal));
  }

  [[nodiscard]] const std::vector<T>& storage() const noexcept { return storage_; }
  [[nodiscard]] std::vector<T>& storage() noexcept { return storage_; }

 private:
  [[nodiscard]] std::size_t linear_index(std::size_t radial, std::size_t polar, std::size_t azimuthal) const {
    if (radial >= radial_ || polar >= polar_ || azimuthal >= azimuthal_) {
      throw std::out_of_range("Array3D index out of range");
    }

    return ((radial * polar_) + polar) * azimuthal_ + azimuthal;
  }

  std::size_t radial_{0};
  std::size_t polar_{0};
  std::size_t azimuthal_{0};
  std::vector<T> storage_;
};

}  // namespace dec3d::core
