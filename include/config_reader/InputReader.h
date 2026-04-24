#pragma once

#include <map>  // for std::map
#include <string>  // for std::string
#include <utility>  // for std::pair
#include <vector>  // for std::vector

#include <libconfig.h++>
#include "mfem.hpp"


namespace mylib {
    /**
     * This class reads in values from a config file.
     *
     * It uses libconfig for the file parsing. LibConfig supports a nested
     * syntax structure where values are referenced by a path. For example,
     * in a config file with the following:
     *
     * mesh_params = {
     *   bounds = {
     *     x_min = 0.0;
     *     x_max = 10.0;
     *     y_min = 0.0;
     *     y_max = 10.0;
     *     z_min = 0.0;
     *     z_max = 10.0;
     *   }
     * }
     *
     * The mesh values would be referenced by the paths "mesh_params.bounds.x_min",
     * "mesh_params.bounds.x_max", etc...
     *
     * To simplify the reading process and to make it easier to change paths,
     * this class can store a path prefix to read nested values.
     *
     * Example usage:
     * double xmin, xmax, ymin, ymax;
     * libconfig::Config cfg;
     * InputReader reader("mesh_params.bounds", &cfg);
     * reader.ReadConfigFile("config.txt")  // must be called before reading values
     * if (reader.ReadValue("x_min", xmin)) { ** process value ** }
     * if (reader.ReadValue("x_max", xmax)) { ** process value ** }
     * ...
     *
     * This can also be used without a prefix by specifying the full path in
     * the ReadValue command:
     * InputReader reader_no_prefix("", &cfg);
     * reader_no_prefix.ReadConfigFile("config.txt")
     * if (reader_no_prefix.ReadValue("mesh_params.bounds.x_min", xmin)) {
     *   ** process value **
     * }
     */
    class InputReader {
    public:
        /**
         * Create object.
         *
         * This takes a pointer to a libconfig Config object. The object does not
         * need to be initialized prior to being passed here. The pointer's target
         * is not copied by this class and the memoery is not freed here. If the
         * object is deleted, this class will not work properly.
         *
         * @param prefix path prefix within the config file
         * @param config pointer to config object
         */
        InputReader(const std::string &prefix, libconfig::Config *config)
            : InputReader(prefix, config, nullptr) {
        }

        /**
         * Create an object from an existing reader
         *
         * This constructor copies the pointer to the libconfig Config object. The
         * intent of this version is to allow easy reading of struct within the
         * config file by creating a reader with a new path.
         *
         * @param prefix path prefix within the config file
         * @param reader InputReader to copy config pointer from
         */
        InputReader(const std::string &prefix, InputReader *reader)
            : InputReader(prefix, reader->cfg, nullptr) {
        }

        /**
         * Create object with a mapping fro tracking read objects.
         *
         * This takes a pointer to a libconfig Config object. The object does not
         * need to be initialized prior to being passed here. The pointer's target
         * is not copied by this class and the memoery is not freed here. If the
         * object is deleted, this class will not work properly.
         *
         * @param prefix path prefix within the config file
         * @param config pointer to config object
         * @param mapping
         */
        InputReader(const std::string &prefix, InputReader *reader,
                    std::map<std::string, bool> *mapping)
            : InputReader(prefix, reader->cfg, mapping) {
        }

        /**
         * Create object with a mapping fro tracking read objects.
         *
         * This takes a pointer to a libconfig Config object. The object does not
         * need to be initialized prior to being passed here. The pointer's target
         * is not copied by this class and the memoery is not freed here. If the
         * object is deleted, this class will not work properly.
         *
         * @param prefix path prefix within the config file
         * @param config pointer to config object
         * @param mapping
         */
        InputReader(const std::string &prefix, libconfig::Config *config,
                    std::map<std::string, bool> *mapping): cfg(config),
                                                           mapping_(mapping),
                                                           rank_(mfem::Mpi::WorldRank()) {
            // If a prefix is given, the path will have a '.' before the variable name.
            // Add this now so it does not need to be added for each read.
            if (!prefix.empty()) {
                prefix_path_ = prefix + ".";
            }
        }

        virtual ~InputReader() = default;

        /**
         * Prints the full libconfig dictionary in a given file
         */
        void PrintConfig(const std::string& filename = "dump.cfg") const {
            if (!rank_) {
                cfg->writeFile(filename.c_str());
            }
        }

        /**
         * Read values from the config file into the file object.
         *
         * This must be called prior to trying to read specific config values.
         *
         * @param filename Name of config file
         */
        void ReadConfigFile(const std::string &filename) {
            try {
                cfg->setAutoConvert(true); // deprecated
                //    this is the correct autoconvert call for the current libconfig version:
                //    cfg.setOptions(libconfig::Config::OptionAutoConvert);
                cfg->readFile(filename.c_str());
            } catch (const libconfig::FileIOException &fioex) {
                MFEM_ABORT("I/O error while reading file.");
            } catch (const libconfig::ParseException &pex) {
                std::stringstream ss;
                ss << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError();
                MFEM_ABORT(ss.str().c_str());
            }
        }

        /**
         * Reads the value of the given key from the configuration file and
         * stores it in the given value variable. The variable and its value is
         * output if the appropriate flags are set.
         *
         * @param key_name The name of the parameter to fetch
         * @param value Variable that will be filled with the value fetched
         * @return true if read is successful
         */
        template<typename T>
        bool ReadValue(const std::string &key_name, T &value) {
            bool result = cfg->lookupValue(prefix_path_ + key_name, value);
            if (result)
                mark_found(key_name);
            if (!rank_) {
                std::cout << "" << prefix_path_ + key_name << ": " << value << " (" << (result ? "config" : "default") << ")\n";
                //PrintStatusStream(std::cerr, "", prefix_path_ + key_name, ": ", value);
            }
            return result;
        }

        /**
         * Reads an array from the input file and stores it in the given variable.
         *
         * The array length is validated against an expected value, thorwing an
         * exception if they don't match.
         *
         * The user is required to pass a pointer to store the the array date. The
         * pointer must be unitialized and the memory will be allocated within the
         * function. The user is responsible for freeing the memory. Memory will only
         * be allocated if the input file has an array of proper length associated
         * with the given key.
         *
         * @param key_name the key name of the configuation value to retrieve
         * @param[out] arr array to store the data (will be allocated here)
         * @param expected_count expected length of the array to read
         * @throw TALYException if key points to a value that is not an array
         * @throw TALYException if array length is incorrect
         * @return true if the array was found in the input file
         */
        template<typename T>
        bool ReadArray(const std::string &key_name, T * &arr, int expected_count) {
            if (!ValidateIsArray(key_name, expected_count)) { return false; }
            libconfig::Setting &setting = cfg->lookup(prefix_path_ + key_name);

            arr = new T[expected_count];
            for (int i = 0; i < expected_count; i++) {
                T value;

                std::string item_path = key_name + ".[" + std::to_string(i) + "]";
                if (!ReadValue(item_path, value)) {
                    return false;
                }

                arr[i] = value;
            }
            mark_found(key_name);
            return true;
        }

        /**
         * Reads an array from the input file and stores it in the given vector.
         *
         * The array length is validated against an expected value, throwing an
         * exception if they don't match.
         *
         * The given vector is resized to fit the data. Any existing data in the
         * vector is lost. Resizing will only happen if the input file has an array
         * of proper length associated with the given key.
         *
         * @param key_name the key name of the configuation value to retrieve
         * @param[out] arr vector to store the data (will be allocated here)
         * @param expected_count expected length of the array to read
         * @throw TALYException if key points to a value that is not an array
         * @throw TALYException if array length is incorrect
         * @return true if the array was found in the input file
         */
        template<typename T>
        bool ReadArray(const std::string &key_name, std::vector<T> &arr,
                       int expected_count) {
            if (!ValidateIsArray(key_name, expected_count)) { return false; }
            libconfig::Setting &setting = cfg->lookup(prefix_path_ + key_name);

            arr.resize(expected_count);
            for (int i = 0; i < expected_count; i++) {
                T value;

                std::string item_path = key_name + ".[" + std::to_string(i) + "]";
                if (!ReadValue(item_path, value)) {
                    return false;
                }

                arr[i] = value;
            }
            mark_found(key_name);
            return true;
        }

        /**
         * Reads an array from the input file and stores it in the given vector.
         *
         * The given vector is resized to fit the data. Any existing data in the
         * vector is lost. Resizing will only happen if the input file has an array
         * of proper length associated with the given key.
         *
         * @param key_name the key name of the configuation value to retrieve
         * @param[out] arr vector to store the data (will be allocated here)
         * @throw TALYException if key points to a value that is not an array
         * @return true if the array was found in the input file
         */
        template<typename T>
        bool ReadArray(const std::string &key_name, std::vector<T> &arr) {
            if (!ValidateIsArray(key_name)) { return false; }
            libconfig::Setting &setting = cfg->lookup(prefix_path_ + key_name);

            int expected_count = setting.getLength();

            arr.clear();
            for (int i = 0; i < expected_count; i++) {
                T value;

                std::string item_path = key_name + ".[" + std::to_string(i) + "]";
                if (!ReadValue(item_path, value)) {
                    return false;
                }

                arr.push_back(value);
            }
            mark_found(key_name);
            return true;
        }

        /**
         * Return the prefix path used by the reader
         *
         * @return prefix path used by the reader
         */
        std::string prefix_path() { return prefix_path_; }

    private:
        /**
         * Mark the given key as found in the reading process
         *
         * @param key key to mark as found
         */
        void mark_found(const std::string &key) {
            if (mapping_) {
                mapping_->insert(std::pair<std::string, bool>(key, true));
            }
        }

        /**
         * Confirms that the given key points to an array.
         *
         * An exception is thrown if this is not an array.
         *
         * @param key_name the key name of the configuation value to retrieve
         * @throw TALYException if key points to a value that is not an array
         * @return true if the key points to an array in the input file
         */
        bool ValidateIsArray(const std::string &key_name) {
            // this is here because lookup throws an exception if not found
            if (!cfg->exists(prefix_path_ + key_name)) {
                return false;
            }

            libconfig::Setting &setting = cfg->lookup(prefix_path_ + key_name);
            if (!setting.isArray()) {
                // confirm this is an array
                std::stringstream ss;
                ss << "expected array input but found single value" << " (key = " << prefix_path_ + key_name << ")";
                MFEM_ABORT(ss.str().c_str());
            }
            return true;
        }

        /**
         * Confirms that the given key points to an array of given length.
         *
         * The array length is validated against an expected value, throwing an
         * exception if they don't match.
         *
         * @param key_name the key name of the configuation value to retrieve
         * @param expected_count expected length of the array to read
         * @throw TALYException if key points to a value that is not an array
         * @throw TALYException if array length is incorrect
         * @return true if a proper length array was found in the input file
         */
        bool ValidateIsArray(const std::string &key_name, int expected_count) {
            // first make sure this is an array
            if (!ValidateIsArray(key_name)) { return false; }

            libconfig::Setting &setting = cfg->lookup(prefix_path_ + key_name);
            if (setting.getLength() != expected_count) {
                // validate length
                std::stringstream ss;
                ss << "input array is wrong length"
                      << " (key = " << prefix_path_ + key_name << ")\n"
                      << "    expected length = " << expected_count
                      << " actual length = " << setting.getLength();
                MFEM_ABORT(ss.str().c_str());
            }
            return true;
        }

        libconfig::Config *cfg; ///< libconfig config object
        std::string prefix_path_; ///< prefix of path in config file
        std::map<std::string, bool> *mapping_ = nullptr; ///< mapping for tracking what was read
        int rank_;
    };
} // namespace TALYFEMLIB
