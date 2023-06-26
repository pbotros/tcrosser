#include <utility>
#include <iostream>
#include <fstream>
#include "include/tcrosser.h"
#include <cassert>

int main(int argc, char **argv) {
    if (argc != 8) {
        std::cout << argv[0] << " <prepeak samples> <postpeak samples> <n channels> <threshold_filename> <threshold multiplier> <data_filename> <output filename>" << std::endl;
        exit(1);
    }

    int prepeak_samples = std::stoi(argv[1]);
    int postpeak_samples = std::stoi(argv[2]);
    int n_channels = std::stoi(argv[3]);

    std::string threshold_filename = argv[4];
    double threshold_multiplier = std::stod(argv[5]);
    std::string data_filename = argv[6];
    std::string output_filename = argv[7];

    // Open and read 1 float at a time:
    std::ifstream threshold_stream(threshold_filename, std::ios::binary | std::ios::in);
    std::vector<float> thresholds;
    for (int i = 0; i < n_channels; i++) {
        float value;
        threshold_stream.read(reinterpret_cast<char *>(&value), sizeof(float));
        thresholds.push_back(value * threshold_multiplier);
    }
    threshold_stream.close();

    std::vector<std::shared_ptr<tcrosser::Thresholder<float>>> thresholders;
    for (int i = 0; i < n_channels; i++) {
        thresholders.push_back(std::make_shared<tcrosser::StaticThresholder<float>>(
            thresholds[i],
            i,
            postpeak_samples));
    }

    tcrosser::ThresholdCrossingCalculator<float> calculator(
        thresholders,
        prepeak_samples,
        postpeak_samples,
        postpeak_samples, // cross channel dead time
        postpeak_samples - 1,
        threshold_multiplier > 0 ? tcrosser::AlignmentDirection::GLOBAL_MAXIMA : tcrosser::AlignmentDirection::GLOBAL_MINIMA
    );

    std::ifstream data_stream(data_filename, std::ios::binary | std::ios::in);

    int64_t chunk_size_samples = 1024;
    std::vector<float> data_chunk_flattened(n_channels * chunk_size_samples);

    std::vector<std::vector<float>> data_chunk(n_channels);
    for (int i = 0; i < n_channels; i++) {
        data_chunk[i].resize(chunk_size_samples);
    }

    int64_t n_bytes_read_total = 0;
    std::ofstream output_stream(output_filename, std::ios::binary | std::ios::out);
    while (true) {
        if (!data_stream) {
            break;
        }
        data_stream.read(
            reinterpret_cast<char *>(data_chunk_flattened.data()),
            sizeof(float) * n_channels * chunk_size_samples);
        auto bytes_read = data_stream.gcount();
        if (bytes_read <= 0) {
            break;
        }
        auto n_samples_read = bytes_read / (sizeof(float) * n_channels);
        if (n_samples_read == 0) {
            break;
        }
        n_bytes_read_total += bytes_read;

        int64_t flattened_sample_index = 0;
        for (int64_t sample_idx = 0; sample_idx < n_samples_read; sample_idx++) {
            for (int64_t channel_idx = 0; channel_idx < n_channels; channel_idx++) {
                data_chunk[channel_idx][sample_idx] = data_chunk_flattened[flattened_sample_index];
                flattened_sample_index++;
            }
        }

        std::vector<const float *> data_by_channel(n_channels);
        for (int i = 0; i < n_channels; i++) {
            data_by_channel[i] = &data_chunk[i][0];
        }
        auto crossings = calculator.consume(data_by_channel.data(), (int) n_samples_read);
        if (!crossings.empty()) {
            for (const auto &item : crossings) {
                output_stream.write(reinterpret_cast<const char *>(&item.data_index),
                                    sizeof(item.data_index));
                output_stream.write(reinterpret_cast<const char *>(&item.aligned_channel_index),
                                    sizeof(item.aligned_channel_index));
            }
        }
    }

    output_stream.close();
    std::cout << "Successfully read " << n_bytes_read_total << " bytes." << std::endl;
}
