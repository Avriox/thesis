//
// Created by Jakob Goldmann on 11.03.25.
// Enhanced version with decorator-like macros for ergonomic timing
//

#ifndef TIMER_H
#define TIMER_H

#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <type_traits>
#include <utility>

// Macro to get function name (compiler-specific)
#if defined(__GNUC__) || defined(__clang__)
    #define FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
    #define FUNCTION_NAME __FUNCSIG__
#else
    #define FUNCTION_NAME __func__
#endif

// Define TIMER_VERBOSE to enable console output for each timed operation
// #define TIMER_VERBOSE

// Forward declaration for the global timer
class Timer;
Timer& get_global_timer();

// Global timer macros (no timer instance needed)
#define TIME(name, func_call) \
    [&]() { \
        auto start = std::chrono::high_resolution_clock::now(); \
        auto result = func_call; \
        auto end = std::chrono::high_resolution_clock::now(); \
        get_global_timer().record(name, start, end); \
        return result; \
    }()

#define TIME_AUTO(func_call) \
    [&]() { \
        auto start = std::chrono::high_resolution_clock::now(); \
        auto result = func_call; \
        auto end = std::chrono::high_resolution_clock::now(); \
        get_global_timer().record(#func_call, start, end); \
        return result; \
    }()

#define TIME_SCOPE(name) \
    auto _timer_scope_guard_##__LINE__ = get_global_timer().make_block_guard(name)

class Timer {
private:
    struct FunctionData {
        std::string name;
        std::chrono::nanoseconds duration;

        FunctionData(const std::string &n) : name(n), duration(0) {}
    };

    // For scoped timing blocks
    class BlockGuard {
    private:
        Timer& timer_;
        std::string name_;
        std::chrono::high_resolution_clock::time_point start_;

    public:
        BlockGuard(Timer& timer, const std::string& name)
            : timer_(timer), name_(name), start_(std::chrono::high_resolution_clock::now()) {}

        ~BlockGuard() {
            auto end = std::chrono::high_resolution_clock::now();
            timer_.record(name_, start_, end);
        }
    };

    std::vector<FunctionData> functions;

public:
    // Record timing directly with start/end times
    void record(const std::string& name,
                std::chrono::high_resolution_clock::time_point start,
                std::chrono::high_resolution_clock::time_point end) {
        functions.emplace_back(name);
        auto& current_func = functions.back();
        current_func.duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

#ifdef TIMER_VERBOSE
        std::cout << "Timed execution of: " << name << " ("
                  << format_duration(current_func.duration) << ")" << std::endl << std::flush;
#endif
    }

    // Create a block guard for scoped timing
    BlockGuard make_block_guard(const std::string& name) {
        return BlockGuard(*this, name);
    }

    // Report methods
    void report_functions() const {
        if (functions.empty()) {
            std::cout << "No functions were timed." << std::endl;
            return;
        }

        size_t max_name_length = std::max_element(functions.begin(), functions.end(),
                                                  [](const FunctionData &a, const FunctionData &b) {
                                                      return a.name.length() < b.name.length();
                                                  }
        )->name.length();

        std::cout << std::endl << std::left << std::setw(max_name_length + 2) << "Function Name" << "| Execution Time"
                  << std::endl;
        std::cout << std::string(max_name_length + 2, '-') << "+" << std::string(20, '-') << std::endl;

        for (const auto &func: functions) {
            std::cout << std::left << std::setw(max_name_length + 2) << func.name << "| "
                      << format_duration(func.duration) << std::endl;
        }
    }

    // One-line report method
    void report() const {
        report_functions();
    }

private:
    static std::string format_duration(const std::chrono::nanoseconds &duration) {
        const auto ns = duration.count();

        if (ns >= 60'000'000'000) {  // More than a minute
            const auto minutes = ns / 60'000'000'000;
            const auto seconds = (ns % 60'000'000'000) / 1'000'000'000;
            return std::to_string(minutes) + "m " + std::to_string(seconds) + "s";
        } else if (ns >= 1'000'000'000) {  // More than a second
            const auto seconds = ns / 1'000'000'000;
            const auto milliseconds = (ns % 1'000'000'000) / 1'000'000;
            return std::to_string(seconds) + "." + std::to_string(milliseconds) + "s";
        } else if (ns >= 1'000'000) {  // More than a millisecond
            return std::to_string(ns / 1'000'000) + "ms";
        } else if (ns >= 1'000) {  // More than a microsecond
            return std::to_string(ns / 1'000) + "µs";
        } else {
            return std::to_string(ns) + "ns";
        }
    }
};

// Global timer implementation as a singleton
Timer& get_global_timer() {
    static Timer instance;
    return instance;
}

// Optional global functions for simpler API
inline void timer_report() {
    get_global_timer().report();
}

#endif //TIMER_H