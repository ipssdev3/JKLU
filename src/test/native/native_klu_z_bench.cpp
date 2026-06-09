#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "klu.h"
}

struct Csc {
    int n;
    std::vector<int> Ap;
    std::vector<int> Ai;
    std::vector<double> Ax;
};

static bool read_data_line(std::istream& in, std::string& line) {
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] != '%') return true;
    }
    return false;
}

static Csc read_matrix(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open matrix");
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("missing matrix header");
    std::istringstream header(line);
    std::string banner, object, format, field, symmetry;
    header >> banner >> object >> format >> field >> symmetry;
    std::transform(banner.begin(), banner.end(), banner.begin(), ::tolower);
    std::transform(object.begin(), object.end(), object.begin(), ::tolower);
    std::transform(format.begin(), format.end(), format.begin(), ::tolower);
    std::transform(field.begin(), field.end(), field.begin(), ::tolower);
    std::transform(symmetry.begin(), symmetry.end(), symmetry.begin(), ::tolower);
    if (banner != "%%matrixmarket" || object != "matrix" || format != "coordinate") {
        throw std::runtime_error("unsupported matrix header");
    }
    if (!read_data_line(in, line)) throw std::runtime_error("missing matrix size");
    int nrow, ncol, nz;
    {
        std::istringstream ss(line);
        ss >> nrow >> ncol >> nz;
    }
    if (nrow != ncol) throw std::runtime_error("matrix must be square");
    struct Entry { int row; int col; double re; double im; };
    std::vector<Entry> entries;
    entries.reserve(nz);
    for (int p = 0; p < nz; p++) {
        if (!read_data_line(in, line)) throw std::runtime_error("missing matrix entry");
        Entry e;
        std::istringstream ss(line);
        ss >> e.row >> e.col;
        if (field == "pattern") {
            e.re = 1.0;
            e.im = 0.0;
        } else if (field == "complex") {
            ss >> e.re >> e.im;
        } else if (field == "real" || field == "integer") {
            ss >> e.re;
            e.im = 0.0;
        } else {
            throw std::runtime_error("unsupported matrix field");
        }
        e.row--;
        e.col--;
        entries.push_back(e);
        if (symmetry != "general" && e.row != e.col) {
            Entry m;
            m.row = e.col;
            m.col = e.row;
            if (symmetry == "symmetric") {
                m.re = e.re;
                m.im = e.im;
            } else if (symmetry == "skew-symmetric") {
                m.re = -e.re;
                m.im = -e.im;
            } else if (symmetry == "hermitian") {
                m.re = e.re;
                m.im = -e.im;
            } else {
                throw std::runtime_error("unsupported matrix symmetry");
            }
            entries.push_back(m);
        }
    }
    std::stable_sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) {
        return a.col == b.col ? a.row < b.row : a.col < b.col;
    });
    Csc csc;
    csc.n = nrow;
    csc.Ap.assign(nrow + 1, 0);
    csc.Ai.resize(entries.size());
    csc.Ax.resize(2 * entries.size());
    for (const auto& e : entries) csc.Ap[e.col + 1]++;
    for (int col = 0; col < nrow; col++) csc.Ap[col + 1] += csc.Ap[col];
    std::vector<int> next = csc.Ap;
    for (const auto& e : entries) {
        int p = next[e.col]++;
        csc.Ai[p] = e.row;
        csc.Ax[2 * p] = e.re;
        csc.Ax[2 * p + 1] = e.im;
    }
    return csc;
}

static std::vector<double> read_rhs(const std::string& path, int n) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open rhs");
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("missing rhs header");
    std::istringstream header(line);
    std::string banner, object, format, field, symmetry;
    header >> banner >> object >> format >> field >> symmetry;
    std::transform(field.begin(), field.end(), field.begin(), ::tolower);
    if (!read_data_line(in, line)) throw std::runtime_error("missing rhs size");
    int rows, cols;
    {
        std::istringstream ss(line);
        ss >> rows >> cols;
    }
    if (rows != n || cols != 1) throw std::runtime_error("unexpected rhs size");
    std::vector<double> rhs(2 * n);
    for (int i = 0; i < n; i++) {
        if (!read_data_line(in, line)) throw std::runtime_error("missing rhs entry");
        std::istringstream ss(line);
        ss >> rhs[2 * i];
        if (field == "complex") {
            ss >> rhs[2 * i + 1];
        } else {
            rhs[2 * i + 1] = 0.0;
        }
    }
    return rhs;
}

static std::vector<double> ones_rhs(const Csc& a) {
    std::vector<double> rhs(2 * a.n);
    for (int col = 0; col < a.n; col++) {
        for (int p = a.Ap[col]; p < a.Ap[col + 1]; p++) {
            int row = a.Ai[p];
            rhs[2 * row] += a.Ax[2 * p];
            rhs[2 * row + 1] += a.Ax[2 * p + 1];
        }
    }
    return rhs;
}

static double ms_since(std::chrono::steady_clock::time_point start,
        std::chrono::steady_clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: native_klu_z_bench matrix.mtx [rhs.mtx|-] [warmup] [iterations]\n";
        return 2;
    }
    bool hasRhs = argc > 2 && std::string(argv[2]) != "-";
    int warmup = argc > 3 ? std::atoi(argv[3]) : 0;
    int iterations = argc > 4 ? std::atoi(argv[4]) : 1;
    Csc a = read_matrix(argv[1]);
    std::vector<double> rhs0 = hasRhs ? read_rhs(argv[2], a.n) : ones_rhs(a);
    bool reuseSymbolic = false;
    if (const char* reuse = std::getenv("NATIVE_KLU_REUSE_SYMBOLIC")) {
        reuseSymbolic = std::atoi(reuse) != 0;
    }

    double analyzeMs = 0.0, factorMs = 0.0, refactorMs = 0.0, solveMs = 0.0;
    double luEntries = 0.0;
    klu_common reusedCommon;
    klu_symbolic* reusedSymbolic = nullptr;
    if (reuseSymbolic) {
        klu_defaults(&reusedCommon);
        if (const char* btf = std::getenv("NATIVE_KLU_BTF")) {
            reusedCommon.btf = std::atoi(btf);
        }
        auto t0 = std::chrono::steady_clock::now();
        reusedSymbolic = klu_analyze(a.n, a.Ap.data(), a.Ai.data(), &reusedCommon);
        auto t1 = std::chrono::steady_clock::now();
        analyzeMs = ms_since(t0, t1);
        if (!reusedSymbolic) {
            std::cerr << "KLU analyze failed status=" << reusedCommon.status << "\n";
            return 1;
        }
    }
    for (int iter = -warmup; iter < iterations; iter++) {
        klu_common common;
        klu_common* commonp = &common;
        if (reuseSymbolic) {
            commonp = &reusedCommon;
        } else {
            klu_defaults(&common);
            if (const char* btf = std::getenv("NATIVE_KLU_BTF")) {
                common.btf = std::atoi(btf);
            }
        }
        auto t0 = std::chrono::steady_clock::now();
        klu_symbolic* symbolic = reusedSymbolic;
        if (!reuseSymbolic) {
            symbolic = klu_analyze(a.n, a.Ap.data(), a.Ai.data(), commonp);
        }
        auto t1 = std::chrono::steady_clock::now();
        klu_numeric* numeric = klu_z_factor(a.Ap.data(), a.Ai.data(), a.Ax.data(), symbolic, commonp);
        auto t2 = std::chrono::steady_clock::now();
        int refactorOk = klu_z_refactor(a.Ap.data(), a.Ai.data(), a.Ax.data(), symbolic, numeric, commonp);
        auto t3 = std::chrono::steady_clock::now();
        std::vector<double> rhs = rhs0;
        int solveOk = klu_z_solve(symbolic, numeric, a.n, 1, rhs.data(), commonp);
        auto t4 = std::chrono::steady_clock::now();
        if (!symbolic || !numeric || refactorOk != 1 || solveOk != 1) {
            std::cerr << "KLU failed status=" << commonp->status << "\n";
            return 1;
        }
        if (iter >= 0) {
            if (!reuseSymbolic) {
                analyzeMs += ms_since(t0, t1);
            }
            factorMs += ms_since(t1, t2);
            refactorMs += ms_since(t2, t3);
            solveMs += ms_since(t3, t4);
            luEntries += numeric->lnz + numeric->unz - a.n +
                    (numeric->Offp ? numeric->Offp[a.n] : 0);
        }
        klu_z_free_numeric(&numeric, commonp);
        if (!reuseSymbolic) {
            klu_free_symbolic(&symbolic, commonp);
        }
    }
    if (reuseSymbolic) {
        klu_free_symbolic(&reusedSymbolic, &reusedCommon);
    }
    double denom = std::max(1, iterations);
    std::cout << "case=native-klu-z,file=" << argv[1]
              << ",n=" << a.n
              << ",nnz=" << a.Ap[a.n]
              << ",iterations=" << iterations
              << ",analyzeMs=" << analyzeMs / denom
              << ",factorMs=" << factorMs / denom
              << ",refactorMs=" << refactorMs / denom
              << ",solveMs=" << solveMs / denom
              << ",luEntries=" << luEntries / denom
              << "\n";
    return 0;
}
