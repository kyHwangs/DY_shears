#include "triggers.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

#include <boost/tokenizer.hpp>

#include <TTree.h>

namespace physics
{
trigger_mask::trigger_mask(util::job::info &info) : _accepts_any_trigger(false)
{}

trigger_mask::trigger_mask(util::job::info &info, const std::string &definition, bool verbose)
    : trigger_mask(info)
{
    using boost::escaped_list_separator;
    using boost::tokenizer;

    int accepted_count = 0;

    // Retrieve tokens
    escaped_list_separator<char> sep("\\", "\t ,\n\r", "\"\'");
    tokenizer<escaped_list_separator<char>> tok(definition, sep);
    for (auto it = tok.begin(); it != tok.end(); ++it) {
        const std::string &token = *it;

        if (token.empty()) {
            continue;
        }

        if (token[0] == '^') { // Trigger is a veto
            if (!veto(info.reader, token.substr(1))) {
                std::string msg = "Could not find trigger ";
                msg += token.substr(1);
                throw std::invalid_argument(msg);
            }
            if (verbose) std::cout << "\tVETO\t" << token.substr(1) << std::endl;
        } else { // Regular trigger
            accepted_count++;
            if (!accept(info.reader, token)) {
                throw std::invalid_argument("Could not find trigger " + token);
            }
            if (verbose) std::cout << "\tACCEPT\t" << token << std::endl;
        }
    }
    if (accepted_count == 0) {
        set_accepts_any_trigger(true);
        if (verbose) std::cout << "\tACCEPT ANY TRIGGER" << std::endl;
    }
}

bool trigger_mask::accept(TTreeReader &tr, const std::string &name)
{
    _accepted_triggers.push_back (TTreeReaderValue<bool> (tr, name.data()));
    return true;
}

bool trigger_mask::accept(TTreeReader &tr, const std::vector<std::string> &names)
{
    bool ok = true;
    for (auto &name : names) {
        ok &= accept(tr, name);
    }
    return ok;
}

bool trigger_mask::veto(TTreeReader &tr, const std::string &name)
{
    _vetoed_triggers.push_back (TTreeReaderValue<bool> (tr, name.data()));
    return true;
}

bool trigger_mask::veto(TTreeReader &tr, const std::vector<std::string> &names)
{
    bool ok = true;
    for (auto &name : names) {
        ok &= veto(tr, name);
    }
    return ok;
}

bool trigger_mask::passes()
{
    if (std::any_of(_vetoed_triggers.begin(),
                    _vetoed_triggers.end(),
                    [] (auto &reader) { return *reader; })) {
        return false;
    }
    if (std::any_of(_accepted_triggers.begin(),
                    _accepted_triggers.end(),
                    [] (auto &reader) { return *reader; })) {
        return true;
    }
    bool pass = _accepts_any_trigger;
    return pass;

}
} // namespace physics
