#ifndef SECTIONEDCONFIG_H
#define SECTIONEDCONFIG_H

#include <fnmatch.h>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <sstream>
#include <stdlib.h>

/** Class to read a configuration file organized in sections.
 */
class SectionedConfig
{
  public:
    typedef std::map<std::string, std::string> rec_t;

    std::map<std::string, rec_t> configMap;

    mutable std::map<std::string, rec_t> retrievedConfig;

    mutable bool read_ = false;

    /** Reads configuration file
     * @parm fname file path
     */
    bool read(const char *fname, bool check = false)
    {

        bool rc = true;

        if (check && read_) return rc;

        if (fname == 0 || fname[0] == 0) return rc;

        const bool debug = false;

        std::ifstream f(fname);
        if (!f.good()) {
            std::cerr << "Error. Configuration file '" << fname << "' was not found\n";
            abort();
        }

        read_ = true;

        configMap.clear();

        std::string line;
        std::regex comment_or_empty_line_pat("[[:space:]]*(#.*)?");
        std::regex section_pat("[[:space:]]*\\[([^\\]]*)\\][[:space:]]*(#.*)?");
        std::regex param_pat("[[:space:]]*([^=[:space:]]*)[[:space:]]*=[[:space:]]*([^=[:space:]]*)"
                             "[[:space:]]*(#.*)?");
        std::string section;
        std::map<std::string, rec_t>::iterator rcd =
            configMap.insert(std::pair<std::string, rec_t>("*", rec_t())).first;
        int iline = 0;
        while (getline(f, line).good()) {
            ++iline;
            if (std::regex_match(line, comment_or_empty_line_pat)) continue;
            std::smatch sm;
            if (debug) std::cout << "Interpreting line " << line << "\n";
            if (std::regex_match(line, sm, section_pat)) {
                if (debug) std::cout << "\tIt's a section header\n";
                section = sm[1];
                if (debug) std::cout << "\tSection name: " << section << "\n";
                rcd = configMap.insert(std::pair<std::string, rec_t>(section, rec_t())).first;
            } else if (std::regex_match(line, sm, param_pat)) {
                if (debug)
                    std::cout << "\tIt's a parameter line\n"
                              << "\tParameter name: " << sm[1] << "\tParameter value: " << sm[2]
                              << "\n";
                rcd->second.operator[](sm[1]) = sm[2];
            } else {
                std::cerr << "Syntax error line " << iline << " of configuration file " << fname
                          << "\n";
            }
        }
        return rc;
    }

    /** Outputs stored configuration to the standard ouput.
     */
    void ls() const
    {
        std::string c = "";
        for (auto r : configMap) {
            std::cout << c << "[" << r.first << "]\n";
            c = "\n";
            for (auto c : r.second) {
                std::cout << c.first << " = " << c.second << "\n";
            }
        }
    }

    /** Outputs retrieved configuration
     * @param stream to write the retrieve configuration to.
     */
    void dumpRetrieved(std::ostream &o) const
    {
        std::map<std::string, std::pair<std::string, bool>> common;
#if 1
        if (retrievedConfig.size() > 1) {
            // search for common paramaters
            // the boolean is set to true if the parameter value is common
            // to all the sections
            for (auto r : retrievedConfig) {
                for (auto p : r.second) {
                    std::pair<std::string, bool> val_and_flag(p.second, true);
                    auto ir = common.insert(std::pair<std::string, std::pair<std::string, bool>>(
                        p.first, val_and_flag));
                    bool inserted = ir.second;
                    std::pair<std::string, bool> &common_rec = ir.first->second;
                    if (!inserted &&
                        common_rec.first !=
                            val_and_flag.first) { // common contains a different value
                        common_rec.second =
                            false; // the value is not common to all sections, invalidate it
                    }
                }
            }

            // keep only common values in the "common" collection
            for (auto it = common.begin(); it != common.end(); ++it) {
                if (!it->second.second) common.erase(it);
            }
        }

        std::string c;
        if (common.size() > 0) {
            o << "[*]\n\n";
            for (auto r : common) {
                o << r.first << " = " << r.second.first << "\n";
            }
            c = "\n";
        }
#else
        std::string c;
#endif

        for (auto r : retrievedConfig) {
            o << c << "[" << r.first << "]\n\n";
            for (auto p : r.second) {
                if (common.find(p.first) == common.end()) {
                    o << p.first << " = " << p.second << "\n";
                }
            }
            c = "\n";
        }
    }

    /** Retrieves a parameter value.
     * @param section section the parameter belong to.
     * @param param parameter name
     * @param [out] parameter value
     * @return true iff the parameter was found
     */
    template <typename T> bool get1(const char *section, const char *param, T &val) const
    {
        const bool debug = false;
        bool found = false;
        std::string sval = "";
        for (auto r : configMap) {
            if (debug) std::cout << "Compare " << r.first.c_str() << " to " << section << "\n";
            if (0 == fnmatch(r.first.c_str(), section, 0)) {
                if (debug) std::cout << "\tThey match.";
                rec_t::iterator it = r.second.find(param);
                if (it != r.second.end()) {
                    found = true;
                    if (debug) std::cout << "\t" << param << " <- " << it->second << "\n";
                    sval = it->second;
                }
            }
        }
        if (found) {
            std::stringstream buf(sval);
            buf >> val;
        } else {
            std::stringstream buf;
            buf << val;
            sval = buf.str();

            // check for consistency of default values used within the code
            auto it1 = retrievedConfig.find(section);
            if (it1 != retrievedConfig.end()) {
                auto it2 = it1->second.find(param);
                if (it2 != it1->second.end()) {
                    std::stringstream buf;
                    if (sval != it2->second) {
                        std::cerr
                            << "Bug found: code uses inconsistent default value for parameter "
                            << param << " of section " << section << ". Aborts.\n";
                        abort();
                    }
                }
            }
        }
        //    std::cout << "adding " << section << ":" << param << " to retrieved config. Size: "
        //<< retrievedConfig.size() << "\n";
        retrievedConfig[std::string(section)][std::string(param)] = sval;
        return found;
    }

    template <typename T> T get(const char *section, const char *param, const T &defaultVal) const
    {
        T val = defaultVal;
        get1(section, param, val);
        return val;
    }

    template <typename T> T get(const char *section, const char *param) const
    {
        T val;
        return get(section, param, val);
    }
};

#endif // SECTIONEDCONFIG_H not defined
