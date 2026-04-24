//
// Created by khara on 10/23/25.
//

#pragma once

#include <libconfig.h++>
#include <string>
#include <sstream>
#include <stdexcept>
#include "mfem.hpp"

namespace mycustom {
    class ConfigReader {
    public:
        explicit ConfigReader(bool verbose = true)
            : verbose_(verbose), rank_(mfem::Mpi::WorldRank()) {
        }

        void load(const std::string& filename) {
            try {
                cfg_.readFile(filename.c_str());
            }
            catch (const libconfig::FileIOException&) {
                throw std::runtime_error("I/O error while reading config file: " + filename);
            }
            catch (const libconfig::ParseException& e) {
                std::ostringstream oss;
                oss << "Parse error in " << e.getFile() << ":" << e.getLine()
                    << " - " << e.getError();
                throw std::runtime_error(oss.str());
            }
        }

        // In ConfigReader
        template <typename T>
        void getRequired(const std::string& path, T& out) const {
            if (!lookupCoerce(path, out)) {
                throw std::runtime_error("Required setting '" + path + "' not found or wrong type.");
            }
            printKV(path, out, ValueSource::Config);
        }

        template <typename T>
        void get(const std::string& path, T& out, const T& default_value) const {
            out = default_value;
            bool found = lookupCoerce(path, out); // stays default if missing
            printKV(path, out, found ? ValueSource::Config : ValueSource::Default);
        }

        bool has(const std::string& path) const {
            try {
                (void)cfg_.lookup(path);
                return true;
            }
            catch (const libconfig::SettingNotFoundException&) {
                return false;
            }
        }

        // Optional: toggle printing at runtime
        void setVerbose(bool v) { verbose_ = v; }

        libconfig::Config& getCfg() {
            return cfg_;
        }

        void printCfg() {
            if (!rank_) {
                cfg_.writeFile("dump.cfg");
            }
        }

    private:
        libconfig::Config cfg_;
        mutable bool verbose_ = true;
        int rank_;

        enum class ValueSource { Default, Config };

        template <typename T>
        void printKV(const std::string& key, const T& value, ValueSource src) const {
            if (!verbose_ || rank_ != 0) return;
            std::cout << key << " = " << toString(value)
                << " (" << (src == ValueSource::Config ? "config" : "default") << ")"
                << std::endl;
        }

        template <typename T>
        static std::string toString(const T& value)
        {
            using U = std::decay_t<T>; // remove cv/ref and decay arrays

            if constexpr (std::is_same_v<U, std::string>) {
                return "\"" + value + "\"";
            }
            else if constexpr (std::is_same_v<U, const char*> || std::is_same_v<U, char*>) {
                return std::string("\"") + (value ? value : "") + "\"";
            }
            else if constexpr (std::is_same_v<U, bool>) {
                return value ? "true" : "false";
            }
            else if constexpr (std::is_integral_v<U>) { // excludes bool due to earlier branch
                std::ostringstream oss;
                oss << value;
                return oss.str();
            }
            else if constexpr (std::is_floating_point_v<U>) {
                std::ostringstream oss;
                oss << std::setprecision(17) << value; // tweak precision as desired
                return oss.str();
            }
            else {
                return "<unprintable>";
            }
        }

        template <typename T>
        bool lookupCoerce(const std::string& path, T& out) const {
            if constexpr (std::is_floating_point_v<T>) {
                // Accept float/double OR ints
                if (cfg_.lookupValue(path, out)) return true;
                long long i64;
                int i;
                if (cfg_.lookupValue(path, i64)) {
                    out = static_cast<T>(i64);
                    return true;
                }
                if (cfg_.lookupValue(path, i)) {
                    out = static_cast<T>(i);
                    return true;
                }
                return false;
            }
            else if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) {
                // Accept int/int64; optionally accept float (narrowing). Remove the double branch to forbid narrowing.
                if (cfg_.lookupValue(path, out)) return true;
                long long i64;
                if (cfg_.lookupValue(path, i64)) {
                    out = static_cast<T>(i64);
                    return true;
                }
                double d;
                if (cfg_.lookupValue(path, d)) {
                    out = static_cast<T>(d);
                    return true;
                }
                return false;
            }
            else if constexpr (std::is_same_v<T, bool> || std::is_same_v<T, std::string>) {
                return cfg_.lookupValue(path, out);
            }
            else {
                // Fallback: exact type only
                return cfg_.lookupValue(path, out);
            }
        }

    };
}
