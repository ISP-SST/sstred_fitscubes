#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <sstream>

#if defined(__linux__)
#include <sys/sysinfo.h>
#elif defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

using namespace std;

// Helper to clean up FITS strings and strip out comments
std::string clean_fits_value(const std::string& str) {
    size_t comment_pos = str.find('/');
    std::string val = (comment_pos != std::string::npos) ? str.substr(0, comment_pos) : str;
    std::string clean_val = "";
    for (char c : val) {
        if (std::isdigit(c) || c == '-' || c == '+') clean_val += c;
    }
    return clean_val;
}

// Automatically fetch available system memory in bytes
int64_t get_available_memory() {
#if defined(__linux__)
    struct sysinfo info;
    if (sysinfo(&info) == 0) return static_cast<int64_t>(info.freeram) * info.mem_unit;
#elif defined(__APPLE__)
    int64_t page_size = 0;
    size_t len = sizeof(page_size);
    sysctlbyname("hw.pagesize", &page_size, &len, NULL, 0);
    uint32_t page_count = 0;
    size_t len_pages = sizeof(page_count);
    sysctlbyname("vm.page_pageable_internal_count", &page_count, &len_pages, NULL, 0);
    if (page_size > 0 && page_count > 0) return static_cast<int64_t>(page_count) * page_size;
#endif
    return 8ULL * 1024 * 1024 * 1024;
}

// Robust FITS structural parsing directly from file header
bool parse_fits_header(const std::string& filename, int& bitpix, std::vector<int64_t>& dims, int64_t& offset_in) {
    std::ifstream in(filename, std::ios::binary);
    if (!in) return false;
    std::vector<char> record_buf(80);
    dims.clear();
    int naxis = 0;
    int64_t bytes_read = 0;
    bool end_found = false;
    while (!end_found && in.read(record_buf.data(), 80)) {
        bytes_read += 80;
        std::string line(record_buf.data(), 80);
        if (line.rfind("SIMPLE  =", 0) == 0) continue;
        if (line.rfind("BITPIX  =", 0) == 0) {
            std::string val = clean_fits_value(line.substr(10, 20));
            if (!val.empty()) bitpix = std::stoi(val);
        } else if (line.rfind("NAXIS   =", 0) == 0) {
            std::string val = clean_fits_value(line.substr(10, 20));
            if (!val.empty()) {
                naxis = std::stoi(val);
                dims.resize(naxis, 1);
            }
        } else if (line.rfind("NAXIS", 0) == 0 && naxis > 0) {
            size_t idx = line.find("=");
            if (idx != std::string::npos) {
                std::string key = line.substr(0, idx);
                key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
                if (key.size() > 5 && key.substr(0, 5) == "NAXIS") {
                    std::string axis_str = key.substr(5);
                    if (std::all_of(axis_str.begin(), axis_str.end(), ::isdigit)) {
                        int axis_num = std::stoi(axis_str);
                        if (axis_num <= naxis) {
                            std::string val = clean_fits_value(line.substr(idx + 1));
                            if (!val.empty()) dims.at(axis_num - 1) = std::stoll(val);
                        }
                    }
                }
            }
        } else if (line.rfind("END     ", 0) == 0) end_found = true;
    }
    offset_in = ((bytes_read - 1) / 2880 + 1) * 2880;
    return (naxis >= 5 && dims.size() >= 5);
}
// (Detta fortsätter direkt efter Del 1 i fitscube_flip_cli.cpp)

template <typename T>
bool do_flip(const std::string& infile, const std::string& outfile, 
             int64_t Nx, int64_t Ny, int64_t Ntuning, int64_t Nstokes, int64_t Nscans,
             int64_t offset_in, int64_t offset_out, int64_t max_memory_bytes) 
{
    std::ifstream in(infile, std::ios::binary);
    std::ofstream out(outfile, std::ios::binary | std::ios::in | std::ios::out);
    if (!in || !out) return false;
    int64_t typesize = sizeof(T);
    int64_t bytes_per_Y_row = Nx * Ntuning * Nstokes * Nscans * typesize;
    int64_t Ny1 = max_memory_bytes / bytes_per_Y_row;
    if (Ny1 < 1) Ny1 = 1;
    if (Ny1 > Ny) Ny1 = Ny;
    int64_t Nchunks = std::ceil(static_cast<double>(Ny) / Ny1);
    int64_t frame_elements = Nx * Ny;
    int64_t total_frames_per_chunk = Nscans * Nstokes * Ntuning;
    auto start_time = std::chrono::steady_clock::now();
    int64_t total_frames_processed = 0;
    int64_t total_frames_to_do = Nchunks * total_frames_per_chunk;
    for (int64_t ichunk = 0; ichunk < Nchunks; ++ichunk) {
        int64_t i0 = ichunk * Ny1;
        int64_t i1 = std::min(i0 + Ny1 - 1, Ny - 1);
        int64_t current_Ny1 = i1 - i0 + 1;
        std::vector<T> subcube(Nx * current_Ny1 * Ntuning * Nstokes * Nscans);
        for (int64_t iscan = 0; iscan < Nscans; ++iscan) {
            for (int64_t istokes = 0; istokes < Nstokes; ++istokes) {
                for (int64_t ituning = 0; ituning < Ntuning; ++ituning) {
                    int64_t iframe = ituning + Ntuning * (istokes + Nstokes * iscan);
                    int64_t file_frame_offset = offset_in + (iframe * frame_elements * typesize) + (i0 * Nx * typesize);
                    in.seekg(file_frame_offset, std::ios::beg);
                    int64_t target_index = (((iscan * Nstokes + istokes) * Ntuning + ituning) * current_Ny1) * Nx;
                    T* dest_ptr = subcube.data() + target_index;
                    in.read(reinterpret_cast<char*>(dest_ptr), Nx * current_Ny1 * typesize);
                    total_frames_processed++;
                    if (total_frames_processed % 10 == 0 || total_frames_processed == total_frames_to_do) {
                        auto current_time = std::chrono::steady_clock::now();
                        double elapsed = std::chrono::duration<double>(current_time - start_time).count();
                        double progress = static_cast<double>(total_frames_processed) / total_frames_to_do;
                        double remaining = (elapsed / progress) - elapsed;
                        int rem_min = static_cast<int>(remaining) / 60;
                        int rem_sec = static_cast<int>(remaining) % 60;
                        printf("\rTransposing data: %.1f%% | Frame %ld/%ld | Remaining: %dm %ds       ", progress * 100.0, total_frames_processed, total_frames_to_do, rem_min, rem_sec);
                        fflush(stdout);
                    }
                }
            }
        }
        std::vector<T> transposed(subcube.size());
        int64_t out_idx = 0;
        for (int64_t iy = 0; iy < current_Ny1; ++iy) {
            for (int64_t ix = 0; ix < Nx; ++ix) {
                for (int64_t istokes = 0; istokes < Nstokes; ++istokes) {
                    for (int64_t iscan = 0; iscan < Nscans; ++iscan) {
                        for (int64_t ituning = 0; ituning < Ntuning; ++ituning) {
                            int64_t in_idx = ix + Nx * (iy + current_Ny1 * (ituning + Ntuning * (istokes + Nstokes * iscan)));
                            *(transposed.data() + out_idx) = subcube.at(in_idx);
                            out_idx++;
                        }
                    }
                }
            }
        }
        int64_t chunk_byte_offset = offset_out + (i0 * Nx * Ntuning * Nstokes * Nscans * typesize);
        out.seekp(chunk_byte_offset, std::ios::beg);
        out.write(reinterpret_cast<const char*>(transposed.data()), transposed.size() * typesize);
    }
    printf("\n"); fflush(stdout);
    return true;
}
// (Detta fortsätter direkt efter Del 2 i fitscube_flip_cli.cpp)

std::string remap_fits_card_axis(const std::string& line, const std::vector<int>& perm) {
    if (line.size() < 8) return line;
    std::string key = line.substr(0, 8);
    std::string clean_key = key;
    clean_key.erase(std::remove_if(clean_key.begin(), clean_key.end(), ::isspace), clean_key.end());
    std::vector<std::string> prefixes = {"NAXIS", "CTYPE", "CUNIT", "CNAME", "CRPIX", "CRVAL", "CDELT", "CSYER", "CWERR", "CWDIS"};
    for (const auto& pre : prefixes) {
        if (clean_key.size() == pre.size() + 1 && clean_key.substr(0, pre.size()) == pre) {
            char axis_char = clean_key.at(pre.size());
            if (std::isdigit(axis_char)) {
                int old_axis = axis_char - '1'; 
                auto it = std::find(perm.begin(), perm.end(), old_axis);
                if (it != perm.end()) {
                    int new_axis = std::distance(perm.begin(), it) + 1; 
                    std::ostringstream oss;
                    oss << std::left << std::setw(5) << pre << new_axis << line.substr(6);
                    std::string new_line = oss.str();
                    if (new_line.size() < 80) new_line.resize(80, ' ');
                    return new_line.substr(0, 80);
                }
            }
        }
    }
    if (clean_key.size() >= 5 && (clean_key.substr(0, 2) == "PS" || clean_key.substr(0, 2) == "PV")) {
        char axis_char = clean_key.at(2);
        if (std::isdigit(axis_char) && clean_key.at(3) == '_') {
            int old_axis = axis_char - '1';
            auto it = std::find(perm.begin(), perm.end(), old_axis);
            if (it != perm.end()) {
                int new_axis = std::distance(perm.begin(), it) + 1;
                std::ostringstream oss;
                oss << clean_key.substr(0, 2) << new_axis << line.substr(3);
                std::string new_line = oss.str();
                if (new_line.size() < 80) new_line.resize(80, ' ');
                return new_line.substr(0, 80);
            }
        }
    }
    return line;
}

void create_mapped_header(const std::string& infile, const std::string& outfile, int64_t offset_out) {
    std::ifstream src(infile, std::ios::binary);
    std::ofstream dst(outfile, std::ios::binary);
    std::vector<char> output_block(offset_out);
    src.read(output_block.data(), offset_out);
    std::vector<int> perm = {2, 4, 3, 0, 1}; 
    std::vector<std::string> output_cards;
    
    // Fixed: Now correctly uses output_block references for the buffer loop
    for (size_t i = 0; i < output_block.size(); i += 80) {
        std::string line(output_block.data() + i, 80);
        std::string processed_line = remap_fits_card_axis(line, perm);
        if (processed_line.rfind("FILENAME=", 0) == 0) {
            std::string out_base = outfile.substr(outfile.find_last_of("/\\") + 1);
            std::ostringstream oss;
            oss << "FILENAME= '" << out_base << "'";
            processed_line = oss.str();
            if (processed_line.size() < 80) processed_line.resize(80, ' ');
            processed_line = processed_line.substr(0, 80);
        }
        output_cards.push_back(processed_line);
    }
    size_t write_offset = 0;
    for (const auto& card : output_cards) {
        if (write_offset + 80 <= static_cast<size_t>(offset_out)) {
            std::copy(card.begin(), card.begin() + 80, output_block.data() + write_offset);
            write_offset += 80;
        }
    }
    dst.write(output_block.data(), offset_out);
}

void append_fits_extensions(const std::string& infile, const std::string& outfile, int64_t data_end_in) {
    std::ifstream src(infile, std::ios::binary | std::ios::ate);
    int64_t file_size = src.tellg();
    int64_t extension_bytes = file_size - data_end_in;
    if (extension_bytes > 0) {
        src.seekg(data_end_in, std::ios::beg);
        std::vector<char> ext_buffer(extension_bytes);
        src.read(ext_buffer.data(), extension_bytes);
        std::ofstream dst(outfile, std::ios::binary | std::ios::app); 
        dst.write(ext_buffer.data(), extension_bytes);
    }
}

bool run_transposition(const std::string& infile, const std::string& outfile, int bitpix, const std::vector<int64_t>& dims, int64_t offset_in, int64_t offset_out, int64_t max_mem) {
    int64_t d_0 = dims.at(0); int64_t d_1 = dims.at(1); int64_t d_2 = dims.at(2); int64_t d_3 = dims.at(3); int64_t d_4 = dims.at(4);
    if (bitpix == 16) return do_flip<int16_t>(infile, outfile, d_0, d_1, d_2, d_3, d_4, offset_in, offset_out, max_mem);
    else if (bitpix == -32) return do_flip<float>(infile, outfile, d_0, d_1, d_2, d_3, d_4, offset_in, offset_out, max_mem);
    return false;
}

bool execute_flip_pipeline(const std::string& infile, const std::string& outfile, int64_t custom_max_mem) {
    int bitpix = 0; std::vector<int64_t> dims; int64_t offset_in = 0;
    if (!parse_fits_header(infile, bitpix, dims, offset_in)) return false;
    int64_t max_mem = (custom_max_mem > 0) ? custom_max_mem : (get_available_memory() / 2);
    int64_t offset_out = offset_in; 
    int64_t total_elements = dims.at(0) * dims.at(1) * dims.at(2) * dims.at(3) * dims.at(4);
    int64_t total_data_bytes = total_elements * (std::abs(bitpix) / 8);
    create_mapped_header(infile, outfile, offset_out);
    { std::ofstream target_dst(outfile, std::ios::binary | std::ios::in | std::ios::out); target_dst.seekp(offset_out + total_data_bytes - 1); char zero = 0; target_dst.write(&zero, 1); }
    bool success = run_transposition(infile, outfile, bitpix, dims, offset_in, offset_out, max_mem);
    if (success) append_fits_extensions(infile, outfile, offset_in + total_data_bytes);
    return success;
}

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cout << "Usage: fitscube_flip <infile.fits> [max_mem_gb]\n"; return 1; }
    
    // Fixed: Native standard C++ string fetching from argv array strings
    std::string infile = argv[1]; 
    std::string outfile = "";
    
    if (infile.size() > 8 && infile.substr(infile.size() - 8) == "_im.fits") outfile = infile.substr(0, infile.size() - 8) + "_sp.fits";
    else if (infile.size() > 5 && infile.substr(infile.size() - 5) == ".fits") outfile = infile.substr(0, infile.size() - 5) + "_sp.fits";
    else outfile = infile + "_sp.fits";
    
    // Fixed: Native std::stoll reading from terminal arg parameters
    int64_t custom_max_mem = (argc > 2) ? (std::stoll(argv[2]) * 1024 * 1024 * 1024) : 0;
    return execute_flip_pipeline(infile, outfile, custom_max_mem) ? 0 : 1;
}
