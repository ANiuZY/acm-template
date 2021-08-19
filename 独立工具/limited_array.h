//
// Created by aniuzhao on 2021/8/19.
//

#ifndef MEMORY_MONITOR_LIMITED_ARRAY_H
#define MEMORY_MONITOR_LIMITED_ARRAY_H

#include <stdint.h>
#include <vector>

template <typename T>
class limited_array {
private:
  int32_t size_ = 0;
  int32_t max_size_ = 10;
  int32_t begin_ = 0;

public:
  std::vector<T> data_;

public:
  limited_array() = default ;
  limited_array(const int32_t max_size) {
    max_size_ = max_size;
  }

private:
  void check(){
    if (begin_ == max_size_) {
      for (int i=0; i < size_; ++i){
        data_[i] = data_[begin_ + i];
      }
      begin_ = 0;
    }
  }

public:
  void push_back(T item) {
    check();
    if (data_.size() <= begin_ + size_) {
      data_.push_back(item);
    } else {
      data_[begin_ + size_] = item;
    }
    if (size_ < max_size_) {
      size_ ++;
    } else {
      begin_ ++;
    }
  }
  void pop_front() {
    check();
    if (size_ > 0) {
      size_ --;
      begin_ ++;
    }
  }
  void pop_back() {
    check();
    if (size_ > 0) {
      size_ --;
    }
  }
  T& operator [] (const int32_t index) {
    return data_[begin_ + index];
  }
  void setMaxSize(const int32_t max_size) {
    max_size_ = max_size;
  }
  int32_t size() const {
    return size_;
  }
};

#endif //MEMORY_MONITOR_LIMITED_ARRAY_H
