//
// Created by Paul Botros on 9/13/20.
//

#ifndef PARENT_TCROSSER_SRC_TCROSSER_H_
#define PARENT_TCROSSER_SRC_TCROSSER_H_

#include <vector>
#include <exception>
#include <sstream>
#include <string>
#include <algorithm>
#include <queue>
#include <memory>

namespace tcrosser {

/**
 * Generic exception used in ThresholdCrossingCalculator.
 */
class ThresholdCrossingException : public std::exception {
 public:
  explicit ThresholdCrossingException(const std::string& message) {
    std::stringstream s;
    s << "[ThresholdCrossingException] " << message;
    _message = s.str();
  }

  const char *what() const noexcept override {
    return _message.c_str();
  }

 private:
  std::string _message;
};

/**
 * Detects threshold crossings on a single channel's stream of data. This class is responsible for tracking its
 * absolute index into the data and for detecting threshold crossings in each chunk of data.
 */
template <class DataT>
class Thresholder {
 public:
  virtual ~Thresholder() = default;

  virtual void set_threshold(DataT threshold) = 0;
  virtual DataT get_threshold() const = 0;
  virtual int get_channel() const = 0;

  /**
   *
   * @param data A chunk of data for this channel.
   * @param num_samples  How many samples are in this chunk of data.
   * @return a vector of *absolute* data indices in which a threshold crossing occurred.
   */
  virtual std::vector<int> threshold(const DataT* data, const int num_samples) = 0;
};

/**
 * Simple implementation of thresholding where the threshold is a static value, and a crossing occurs whenever the
 * previous value does not exceed the threshold while the current value does. Allows for a "dead time", which is the
 * minimum number of samples in between detected threshold crossings.
 *
 * @tparam DataT the type of the multi-channel data stream.
 */
template <class DataT>
class StaticThresholder : public Thresholder<DataT> {
 public:
  StaticThresholder(
      DataT threshold,
      int channel_index,
      int deadtime_samples) {
    threshold_ = threshold;
    channel_index_ = channel_index;
    deadtime_samples_ = deadtime_samples;
    cur_idx = 0;
    next_okay_index_ = 0;
    last_sample_ = 0;
  }

  void set_threshold(DataT threshold) override {
    threshold_ = threshold;
  }

  DataT get_threshold() const override {
    return threshold_;
  }

  int get_channel() const override {
    return channel_index_;
  }

  std::vector<int> threshold(const DataT* data, const int num_samples) override {
    std::vector<int> ret;

    // Cache locally in case the class-level var gets modified during the execution of this branch.
    DataT threshold = threshold_;

    for (int i = 0; i < num_samples; i++) {
      int64_t global_idx = cur_idx + i;

      bool is_tc;
      DataT prev_sample;
      if (i == 0) {
        if (global_idx == 0) {
          // last_sample_ doesn't mean anything on the very first sample, so skip entirely
          continue;
        }
        // If it's the first sample of this batch, have to look back at the previous batch:
        prev_sample = last_sample_;
      } else {
        prev_sample = data[i - 1];
      }

      DataT sample = data[i];
      if (threshold > 0) {
        is_tc = sample >= threshold && prev_sample < threshold;
      } else {
        is_tc = sample <= threshold && prev_sample > threshold;
      }
      if (is_tc && global_idx >= next_okay_index_) {
        next_okay_index_ = global_idx + deadtime_samples_ + 1;
        ret.push_back(global_idx);
      }

      if (i == num_samples - 1) {
        last_sample_ = sample;
      }
    }
    cur_idx += num_samples;
    return ret;
  }

 private:
  DataT threshold_;
  DataT last_sample_;
  int deadtime_samples_;
  int64_t cur_idx;
  int64_t next_okay_index_;
  int channel_index_;
};

/**
 * How to align threshold crossings waveforms.
 */
enum class AlignmentDirection {
  GLOBAL_MAXIMA,
  GLOBAL_MINIMA,
};

template <class DataT>
class ThresholdCrossing {
 public:
  int64_t data_index;
  int aligned_channel_index;
  int alignment_samples;

  // of dimensions n_channels (first dim) x n_samples (second dim)
  std::vector<std::vector<DataT>> multichannel_waveform;

  ThresholdCrossing(
      int64_t data_index,
      int aligned_channel_index,
      int alignment_samples,
      std::vector<std::vector<DataT>> multichannel_waveform) :
      data_index(data_index),
      aligned_channel_index(aligned_channel_index),
      alignment_samples(alignment_samples),
      multichannel_waveform(multichannel_waveform) {}

  struct less : std::binary_function<
      ThresholdCrossing<DataT>, ThresholdCrossing<DataT>, bool> {
    bool operator() (ThresholdCrossing<DataT> a, ThresholdCrossing<DataT> b) const {
      return a.data_index < b.data_index;
    }
  };
};

template <class DataT>
class ThresholdCrossingCalculator {
 public:
  ThresholdCrossingCalculator(
      std::vector<std::shared_ptr<Thresholder<DataT>>> thresholders,
      int prepeak_samples,
      int postpeak_samples,
      int cross_channel_deadtime_samples,
      int max_alignment_window_samples,
      AlignmentDirection alignment_direction,
      int buffer_size_per_channel = -1) {
    thresholders_ = thresholders;
    num_channels_ = thresholders.size();
    prepeak_samples_ = prepeak_samples;
    postpeak_samples_ = postpeak_samples;
    max_alignment_window_samples_ = max_alignment_window_samples;
    alignment_direction_ = alignment_direction;
    buffer_idx_ = 0;
    next_okay_index_ = 0;
    cross_channel_deadtime_samples_ = cross_channel_deadtime_samples;

    // Some arbitrary multiple of the maximum number of samples we'll need to remember
    if (buffer_size_per_channel <= 0) {
      buffer_size_per_channel_ = 16 * (
          prepeak_samples_ +
              postpeak_samples_ +
              max_alignment_window_samples_);
    } else {
      buffer_size_per_channel_ = buffer_size_per_channel;
    }

    buffer_ = std::vector<DataT>(buffer_size_per_channel_ * num_channels_);

    if (num_channels_ <= 0) {
      throw ThresholdCrossingException("Number of channels must be positive.");
    }
    if (prepeak_samples_ <= 0) {
      throw ThresholdCrossingException("prepeak_samples must be positive.");
    }
    if (postpeak_samples_ <= 0) {
      throw ThresholdCrossingException("postpeak_samples must be positive.");
    }
    if (max_alignment_window_samples_ < 0) {
      throw ThresholdCrossingException("max_alignment_window_samples must be nonnegative.");
    }
  };

  /**
   * Consume a new chunk of multi-channel data, and return any new threshold crossings. The returned threshold crossings
   * don't necessarily need to entirely lie in the new chunk of data.
   *
   * @param data_by_channel Data by "channel", passed in the same order as the Thresholders's in the constructor.
   * @param num_samples_per_channel How many samples *per channel* are in this chunk of data.
   * @return
   */
  std::vector<ThresholdCrossing<DataT>> consume(const DataT **data_by_channel, const int num_samples_per_channel) {
    // Copy over the new data into our circular buffer, and detect threshold crossings while we're at it
    for (int chidx = 0; chidx < num_channels_; chidx++) {
      const DataT *channel_data = data_by_channel[chidx];
      for (int j = 0; j < num_samples_per_channel; j++) {
        *buffer_at(chidx, buffer_idx_ + j) = channel_data[j];
      }

      auto thresholder = thresholders_[chidx];
      std::vector<int> channel_crossings = thresholder->threshold(channel_data, num_samples_per_channel);
      // Push our pending threshold crossings to a priority queue, so that we always process our threshold crossings in
      // order.
      for (const auto &crossing_index : channel_crossings) {
        pending_crossings_ch_time_.push(std::make_pair(chidx, crossing_index));
      }
    }

    int64_t max_alignable_index =
        buffer_idx_ + num_samples_per_channel - max_alignment_window_samples_ - postpeak_samples_;
    buffer_idx_ += num_samples_per_channel;

    std::vector<ThresholdCrossing<DataT>> ret;

    // Process any threshold crossings that we can, i.e. crossings that occurred enough samples ago. Since we process
    // in ascending order, we can stop once we hit the first crossing that occurred too recently.
    while (!pending_crossings_ch_time_.empty()) {
      auto crossing_ch_time = pending_crossings_ch_time_.top();
      int crossing_chidx = crossing_ch_time.first;
      int crossing_time = crossing_ch_time.second;

      if (crossing_time >= max_alignable_index) {
        break;
      }

      pending_crossings_ch_time_.pop();

      if (crossing_time < next_okay_index_) {
        // A threshold crossing that's too close to the last one - ignore!
        continue;
      }

      next_okay_index_ = crossing_time + cross_channel_deadtime_samples_ + 1;

      // Compute the "alignment offset", i.e. how  many samples from the threshold crossing to the aligned peak
      int offset = 0;
      int aligned_channel_index = 0;
      switch (alignment_direction_) {
        case AlignmentDirection::GLOBAL_MAXIMA: {
          DataT max_value = -1e10;
          // TODO(multichannel): sweep across *all* channels, not just the one on which there
          // was a crossing
          for (int channel_idx = 0; channel_idx < num_channels_; channel_idx++) {
            for (int sample_idx = 0; sample_idx < max_alignment_window_samples_; sample_idx++) {
              DataT sample = *buffer_at(channel_idx, crossing_time + sample_idx);
              if (sample > max_value) {
                max_value = sample;
                offset = sample_idx;
                aligned_channel_index = channel_idx;
              }
            }
          }
        }
          break;
        case AlignmentDirection::GLOBAL_MINIMA: {
          DataT min_value = 1e10;
          for (int channel_idx = 0; channel_idx < num_channels_; channel_idx++) {
            for (int sample_idx = 0; sample_idx < max_alignment_window_samples_; sample_idx++) {
              DataT sample = *buffer_at(channel_idx, crossing_time + sample_idx);
              if (sample < min_value) {
                min_value = sample;
                offset = sample_idx;
                aligned_channel_index = channel_idx;
              }
            }
          }
        }
          break;
        default: {
          throw ThresholdCrossingException("Unhandled alignment direction!");
        }
      }

      int64_t peak_index = offset + crossing_time;

      // Grab a full window from [-prepeak_samples, postpeak_samples) according to the peak_index
      int waveform_samples = prepeak_samples_ + postpeak_samples_;

      // Invalidate any waveforms that are "too early"
      int64_t waveform_start_samples = peak_index - prepeak_samples_;
      if (waveform_start_samples < 0) {
        continue;
      }

      // Copy over the aligned waveform
      std::vector<std::vector<DataT>> multichannel_waveform;
      for (int channel_index = 0; channel_index < num_channels_; channel_index++) {
        std::vector<DataT> waveform(waveform_samples);
        for (int i = 0; i < waveform_samples; i++) {
          int64_t sample_idx = waveform_start_samples + i;
          waveform[i] = *buffer_at(channel_index, sample_idx);
        }
        multichannel_waveform.push_back(waveform);
      }

      // Translate our local "channel index" into an externally-facing channel index
      int external_chidx = thresholders_[aligned_channel_index]->get_channel();
      ret.push_back(ThresholdCrossing<DataT>(
          waveform_start_samples,
          external_chidx,
          offset,
          multichannel_waveform));
    }

    // Sort *after* alignment so that we always guarantee that we emit TCs in ascending order.
    std::sort(ret.begin(), ret.end(), typename ThresholdCrossing<DataT>::less());
    return ret;
  }

  std::shared_ptr<Thresholder<DataT>> thresholder(int channel_index) {
    for (const auto& thresholder : thresholders_) {
      if (thresholder->get_channel() == channel_index) {
        return thresholder;
      }
    }
    return std::shared_ptr<Thresholder<DataT>>(nullptr);
  }

 private:
  std::vector<std::shared_ptr<Thresholder<DataT>>> thresholders_;
  int num_channels_;
  int prepeak_samples_;
  int postpeak_samples_;
  int max_alignment_window_samples_;
  int buffer_size_per_channel_;
  int next_okay_index_;
  int cross_channel_deadtime_samples_;
  AlignmentDirection alignment_direction_;

  std::vector<DataT> buffer_;
  int64_t buffer_idx_;

  struct ch_time_greater : std::binary_function<
      std::pair<int, int>, std::pair<int, int>, bool> {
    bool operator() (std::pair<int, int> a, std::pair<int, int> b) const {
      if (a.second == b.second) {
        return a.first > b.first;
      }
      return a.second > b.second;
    }
  };

  std::priority_queue<std::pair<int, int>,
                      std::vector<std::pair<int, int>>,
                      ch_time_greater> pending_crossings_ch_time_;

  inline DataT* buffer_at(int channel_index, int64_t global_sample_index) {
    // Buffer is indexed as <n_channels, buffer_size_per_channel>, i.e., each channel is its own row
    return &buffer_.front()
        + channel_index * buffer_size_per_channel_
        + (global_sample_index % buffer_size_per_channel_);
  }
};

}

#endif //PARENT_TCROSSER_SRC_TCROSSER_H_
