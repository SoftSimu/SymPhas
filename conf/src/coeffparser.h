#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>

// Parses the new PFC coefficient file format
// Returns a flattened vector of coefficients
inline std::vector<double> parse_pfc_coeff_file(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Could not open coefficients file: " + filename);

    size_t num_fields = 0;
    double default_value = 1.0;
    std::vector<std::vector<double>> non_coupling;
    std::vector<std::vector<double>> coupling;
    enum Section { NONE, DEFAULT, NONCOUPLING, COUPLING } section = NONE;
    std::string line;
    while (std::getline(infile, line)) {
        // Remove comments
        auto hash = line.find('#');
        if (hash != std::string::npos) line = line.substr(0, hash);
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        if (line.empty()) continue;
        if (line.find("fields:") == 0) {
            num_fields = std::stoul(line.substr(7));
            continue;
        }
        if (line == "[default]") { section = DEFAULT; continue; }
        if (line == "[non_coupling]") { section = NONCOUPLING; continue; }
        if (line == "[coupling]") { section = COUPLING; continue; }
        if (section == DEFAULT) {
            default_value = std::stod(line);
            section = NONE;
            continue;
        }
        std::vector<double> values;
        std::stringstream ss(line);
        std::string item;
        while (std::getline(ss, item, ',')) {
            item.erase(0, item.find_first_not_of(" \t"));
            item.erase(item.find_last_not_of(" \t") + 1);
            if (!item.empty()) values.push_back(std::stod(item));
        }
        if (section == NONCOUPLING) non_coupling.push_back(values);
        else if (section == COUPLING) coupling.push_back(values);
    }
    // Validation
    if (num_fields == 0) throw std::runtime_error("fields: not specified in coefficients file");
    if (non_coupling.size() != num_fields)
        throw std::runtime_error("Number of [non_coupling] rows does not match fields count");
    size_t expected_coupling = num_fields * (num_fields - 1) / 2;
    if (coupling.size() != expected_coupling)
        throw std::runtime_error("Number of [coupling] rows does not match fields*(fields-1)/2");
    // Flatten
    std::vector<double> flat;
    for (const auto& row : non_coupling) flat.insert(flat.end(), row.begin(), row.end());
    for (const auto& row : coupling) flat.insert(flat.end(), row.begin(), row.end());
    // Fill missing with default
    for (auto& v : flat) if (std::isnan(v)) v = default_value;
    return flat;
}
