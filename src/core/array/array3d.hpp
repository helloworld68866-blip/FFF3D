#pragma once

#include <cassert>
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

  T& operator()(std::size_t radial, std::size_t polar, std::size_t azimuthal) noexcept {
    assert(in_bounds(radial, polar, azimuthal));
    return storage_[linear_index_unchecked(radial, polar, azimuthal)];
  }

  const T& operator()(std::size_t radial, std::size_t polar, std::size_t azimuthal) const noexcept {
    assert(in_bounds(radial, polar, azimuthal));
    return storage_[linear_index_unchecked(radial, polar, azimuthal)];
  }

  T& checked(std::size_t radial, std::size_t polar, std::size_t azimuthal) {
    return storage_.at(linear_index_checked(radial, polar, azimuthal));
  }

  const T& checked(std::size_t radial, std::size_t polar, std::size_t azimuthal) const {
    return storage_.at(linear_index_checked(radial, polar, azimuthal));
  }

  [[nodiscard]] std::size_t linear_index_unchecked(
      std::size_t radial,
      std::size_t polar,
      std::size_t azimuthal) const noexcept {
    return ((radial * polar_) + polar) * azimuthal_ + azimuthal;
  }

  [[nodiscard]] T* data() noexcept { return storage_.data(); }
  [[nodiscard]] const T* data() const noexcept { return storage_.data(); }
  [[nodiscard]] const std::vector<T>& storage() const noexcept { return storage_; }
  [[nodiscard]] std::vector<T>& storage() noexcept { return storage_; }

 private:
  [[nodiscard]] bool in_bounds(std::size_t radial, std::size_t polar, std::size_t azimuthal) const noexcept {
    return radial < radial_ && polar < polar_ && azimuthal < azimuthal_;
  }

  [[nodiscard]] std::size_t linear_index_checked(
      std::size_t radial,
      std::size_t polar,
      std::size_t azimuthal) const {
    if (radial >= radial_ || polar >= polar_ || azimuthal >= azimuthal_) {
      throw std::out_of_range("Array3D index out of range");
    }

    return linear_index_unchecked(radial, polar, azimuthal);
  }

  std::size_t radial_{0};
  std::size_t polar_{0};
  std::size_t azimuthal_{0};
  std::vector<T> storage_;
};

}  // namespace dec3d::core
