#include <utility>

#include "gtest/gtest.h"
#include "tcrosser.h"

using namespace tcrosser;

class ThresholdCrossingTest : public ::testing::Test {
 protected:
  void SetUp() override {
  }
};

// Constants fixed for this test.
const static int _BUFFER_SIZE = 4096;
const static int _WAVEFORM_PREPEAK_SAMPLES = 16;
const static int _WAVEFORM_POSTPEAK_SAMPLES = 48;
const static int _MAX_ALIGNMENT_WINDOW = _WAVEFORM_POSTPEAK_SAMPLES - 1;

class FakeThresholder : public Thresholder<double> {
 public:
  explicit FakeThresholder(std::queue<int> tc_indices, int channel_index)
      : data_index_(0), channel_index_(channel_index), tc_indices_(std::move(tc_indices)) {}

  std::vector<int> threshold(const double *, const int num_samples) override {
    std::vector<int> ret;
    while (!tc_indices_.empty()) {
      int tc_index = tc_indices_.front();
      if (tc_index < data_index_ + num_samples) {
        tc_indices_.pop();
        ret.push_back(tc_index);
      } else {
        break;
      }
    }
    data_index_ += num_samples;
    return ret;
  }

  void set_threshold(double threshold) override {}

  double get_threshold() const override {
    return 0;
  }

  int get_channel() const override {
    return channel_index_;
  }

 private:
  int data_index_;
  int channel_index_;
  std::queue<int> tc_indices_;
};

static std::vector<ThresholdCrossing<double>> do_threshold_crossings(std::vector<std::vector<int>> tc_indices, double **fake_data, int num_samples, int num_channels) {
  std::vector<std::shared_ptr<Thresholder<double>>> thresholders(num_channels);
  for (int chidx = 0; chidx < num_channels; chidx++) {
    std::queue<int> tc_indices_for_channel;

    std::sort(tc_indices[chidx].begin(), tc_indices[chidx].end());
    for (auto tc_index : tc_indices[chidx]) {
      tc_indices_for_channel.push(tc_index);
    }

    thresholders[chidx] = std::make_shared<FakeThresholder>(tc_indices_for_channel, chidx);
  }

  ThresholdCrossingCalculator<double> calculator(
      thresholders,
      _WAVEFORM_PREPEAK_SAMPLES,
      _WAVEFORM_POSTPEAK_SAMPLES,
      _WAVEFORM_POSTPEAK_SAMPLES,
      _MAX_ALIGNMENT_WINDOW,
      AlignmentDirection::GLOBAL_MAXIMA,
      _BUFFER_SIZE);

  std::vector<ThresholdCrossing<double>> ret;

  // Consume in arbitrarily small batches of 256 samples at a time to test our streaming logic.
  int idx = 0;
  while (idx < num_samples) {
    auto **fake_data_chunk = new const double*[num_channels];
    for (int chidx = 0; chidx < num_channels; chidx++) {
      fake_data_chunk[chidx] = fake_data[chidx] + idx;
    }

    int num_to_align = std::min(num_samples - idx, 256);
    auto crossings = calculator.consume(fake_data_chunk, num_to_align);
    ret.insert(ret.end(), crossings.begin(), crossings.end());

    // Advance each channel pointer in our fake data
    idx += num_to_align;
    delete[] fake_data_chunk;
  }

  return ret;
}

TEST_F(ThresholdCrossingTest, TestWorks) {
  int total_samples = 10000;
  int num_channels = 1;

  auto **fake_data = new double*[num_channels];
  double value = 0.0;
  for (int chidx = 0; chidx < num_channels; chidx++) {
    fake_data[chidx] = new double[total_samples];
    for (int j = 0; j < total_samples; j++) {
      fake_data[chidx][j] = value;
      value += 1.0;
    }
  }

  std::vector<std::pair<int, int>> spike_test_cases_ch_time(
      {{
           // 1. Spike on 1st and 5th channel in third batch of reading
           {0, 80},
           {0, 160},
           // 2. Spike on 16th channel that is placed right before, at, and after "wraparound" of the buffer
           {0, _BUFFER_SIZE - _WAVEFORM_POSTPEAK_SAMPLES - 5},
           {0, _BUFFER_SIZE},
           {0, _BUFFER_SIZE + _WAVEFORM_POSTPEAK_SAMPLES + 5},
       }});

  std::vector<std::vector<int>> tc_indices(num_channels);
  for (const auto& ch_time : spike_test_cases_ch_time) {
    tc_indices[ch_time.first].push_back(ch_time.second);
    // Make it a global maximum at the start of the waveform window
    fake_data[ch_time.first][ch_time.second] += 500;
  }

  auto crossings = do_threshold_crossings(tc_indices, fake_data, total_samples, num_channels);
  ASSERT_EQ(crossings.size(), spike_test_cases_ch_time.size());
  for (int i = 0; i < spike_test_cases_ch_time.size(); i++) {
    const auto &ch_time = spike_test_cases_ch_time[i];
    ThresholdCrossing<double> &crossing = crossings[i];
    ASSERT_EQ(crossing.aligned_channel_index, ch_time.first);
    ASSERT_EQ(crossing.data_index, ch_time.second - _WAVEFORM_PREPEAK_SAMPLES);
    ASSERT_EQ(crossing.alignment_samples, 0);
    ASSERT_EQ(crossing.multichannel_waveform.size(), num_channels);
    auto waveform = crossing.multichannel_waveform[0];
    for (int j = 0; j < waveform.size(); j++) {
      int expected_sample_idx = ch_time.second - _WAVEFORM_PREPEAK_SAMPLES + j;
      ASSERT_DOUBLE_EQ(waveform[j], fake_data[0][expected_sample_idx]);
    }
  }

  for (int chidx = 0; chidx < num_channels; chidx++) {
    delete[] fake_data[chidx];
  }
  delete[] fake_data;
}

TEST_F(ThresholdCrossingTest, TestDedupesSpikesAcrossChannels) {
  int total_samples = 10000;
  int num_channels = 16;
  int alignment_samples_to_test = 3;

  auto **fake_data = new double*[num_channels];
  for (int chidx = 0; chidx < num_channels; chidx++) {
    double value = 0.0;
    fake_data[chidx] = new double[total_samples];
    for (int j = 0; j < total_samples; j++) {
      fake_data[chidx][j] = value;
      value += 1.0;
    }
  }

  std::vector<std::pair<int, int>> spike_test_cases_ch_time(
      {{
           // 1. A bunch of spikes that should all be collapsed into 1
           {0, 80},
           {4, 80},
           {5, 80 + 1},
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES},
           // 2. Just after the threshold:
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES + 1},
       }});

  std::vector<std::pair<int, int>> spike_test_cases_expected_ch_time(
      {{
           {0, 80},
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES + 1},
       }});

  std::vector<std::vector<int>> tc_indices(num_channels);
  for (const auto& ch_time : spike_test_cases_expected_ch_time) {
    tc_indices[ch_time.first].push_back(ch_time.second);
    // Make it a global maximum at the start of the waveform window
    fake_data[ch_time.first][ch_time.second + alignment_samples_to_test] += float(5e6);
  }

  auto crossings = do_threshold_crossings(tc_indices, fake_data, total_samples, num_channels);
  ASSERT_EQ(crossings.size(), spike_test_cases_expected_ch_time.size());
  for (int i = 0; i < spike_test_cases_expected_ch_time.size(); i++) {
    const auto &ch_time = spike_test_cases_expected_ch_time[i];
    ThresholdCrossing<double> &crossing = crossings[i];
    ASSERT_EQ(crossing.aligned_channel_index, ch_time.first);
    ASSERT_EQ(crossing.data_index, ch_time.second - _WAVEFORM_PREPEAK_SAMPLES + alignment_samples_to_test);
    ASSERT_EQ(crossing.alignment_samples, alignment_samples_to_test);
    ASSERT_EQ(crossing.multichannel_waveform.size(), num_channels);
    for (int chidx = 0; chidx < num_channels; chidx++) {
      auto waveform = crossing.multichannel_waveform[chidx];
      for (int j = 0; j < waveform.size(); j++) {
        auto expected_sample_idx = crossing.data_index + j;
        ASSERT_DOUBLE_EQ(waveform[j], fake_data[chidx][expected_sample_idx]);
      }
    }
  }

  for (int chidx = 0; chidx < num_channels; chidx++) {
    delete[] fake_data[chidx];
  }
  delete[] fake_data;
}


TEST_F(ThresholdCrossingTest, TestWorksIntegration) {
  int total_samples = 10000;
  int num_channels = 1;

  auto **fake_data = new double*[num_channels];
  double value = 0.0;
  for (int chidx = 0; chidx < num_channels; chidx++) {
    fake_data[chidx] = new double[total_samples];
    for (int j = 0; j < total_samples; j++) {
      fake_data[chidx][j] = value;
      value += 1.0;
    }
  }

  std::vector<std::pair<int, int>> spike_test_cases_ch_time(
      {{
           // 1. Spike on 1st and 5th channel in third batch of reading
           {0, 80},
           {0, 160},
           // 2. Spike on 16th channel that is placed right before, at, and after "wraparound" of the buffer
           {0, _BUFFER_SIZE - _WAVEFORM_POSTPEAK_SAMPLES - 5},
           {0, _BUFFER_SIZE},
           {0, _BUFFER_SIZE + _WAVEFORM_POSTPEAK_SAMPLES + 5},
       }});

  std::vector<std::vector<int>> tc_indices(num_channels);
  for (const auto& ch_time : spike_test_cases_ch_time) {
    tc_indices[ch_time.first].push_back(ch_time.second);
    // Make it a global maximum at the start of the waveform window
    fake_data[ch_time.first][ch_time.second] += 500;
  }

  auto crossings = do_threshold_crossings(tc_indices, fake_data, total_samples, num_channels);
  ASSERT_EQ(crossings.size(), spike_test_cases_ch_time.size());
  for (int i = 0; i < spike_test_cases_ch_time.size(); i++) {
    const auto &ch_time = spike_test_cases_ch_time[i];
    ThresholdCrossing<double> &crossing = crossings[i];
    ASSERT_EQ(crossing.aligned_channel_index, ch_time.first);
    ASSERT_EQ(crossing.data_index, ch_time.second - _WAVEFORM_PREPEAK_SAMPLES);
    ASSERT_EQ(crossing.alignment_samples, 0);
    ASSERT_EQ(crossing.multichannel_waveform.size(), num_channels);
    auto waveform = crossing.multichannel_waveform[0];
    for (int j = 0; j < waveform.size(); j++) {
      int expected_sample_idx = ch_time.second - _WAVEFORM_PREPEAK_SAMPLES + j;
      ASSERT_DOUBLE_EQ(waveform[j], fake_data[0][expected_sample_idx]);
    }
  }

  for (int chidx = 0; chidx < num_channels; chidx++) {
    delete[] fake_data[chidx];
  }
  delete[] fake_data;
}

TEST_F(ThresholdCrossingTest, TestDedupeSpikesAcrossChannels) {
  int total_samples = 10000;
  int num_channels = 16;
  int alignment_samples_to_test = 3;

  auto **fake_data = new double*[num_channels];
  for (int chidx = 0; chidx < num_channels; chidx++) {
    double value = 0.0;
    fake_data[chidx] = new double[total_samples];
    for (int j = 0; j < total_samples; j++) {
      fake_data[chidx][j] = value;
      value += 1.0;
    }
  }

  std::vector<std::pair<int, int>> spike_test_cases_ch_time(
      {{
           // 1. A bunch of spikes that should all be collapsed into 1
           {0, 80},
           {4, 80},
           {5, 80 + 1},
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES},
           // 2. Just after the threshold:
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES + 1},
       }});

  std::vector<std::pair<int, int>> spike_test_cases_expected_ch_time(
      {{
           {0, 80},
           {0, 80 + _WAVEFORM_POSTPEAK_SAMPLES + 1},
       }});

  std::vector<std::vector<int>> tc_indices(num_channels);
  for (const auto& ch_time : spike_test_cases_expected_ch_time) {
    tc_indices[ch_time.first].push_back(ch_time.second);
    // Make it a global maximum at the start of the waveform window
    fake_data[ch_time.first][ch_time.second + alignment_samples_to_test] += float(5e6);
  }

  auto crossings = do_threshold_crossings(tc_indices, fake_data, total_samples, num_channels);
  ASSERT_EQ(crossings.size(), spike_test_cases_expected_ch_time.size());
  for (int i = 0; i < spike_test_cases_expected_ch_time.size(); i++) {
    const auto &ch_time = spike_test_cases_expected_ch_time[i];
    ThresholdCrossing<double> &crossing = crossings[i];
    ASSERT_EQ(crossing.aligned_channel_index, ch_time.first);
    ASSERT_EQ(crossing.data_index, ch_time.second - _WAVEFORM_PREPEAK_SAMPLES + alignment_samples_to_test);
    ASSERT_EQ(crossing.alignment_samples, alignment_samples_to_test);
    ASSERT_EQ(crossing.multichannel_waveform.size(), num_channels);
    for (int chidx = 0; chidx < num_channels; chidx++) {
      auto waveform = crossing.multichannel_waveform[chidx];
      for (int j = 0; j < waveform.size(); j++) {
        auto expected_sample_idx = crossing.data_index + j;
        ASSERT_DOUBLE_EQ(waveform[j], fake_data[chidx][expected_sample_idx]);
      }
    }
  }

  for (int chidx = 0; chidx < num_channels; chidx++) {
    delete[] fake_data[chidx];
  }
  delete[] fake_data;
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
