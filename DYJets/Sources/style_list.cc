#include "style_list.h"

#include <stdexcept>

namespace util
{

style_list::style_list(const YAML::Node &specification)
    : _rules(specification.as<std::vector<style_list::rule>>())
{
}

std::string style_list::get_formatted(const std::string &name,
                                      const std::string &rule_match,
                                      const std::string &default_value) const
{
    for (auto it = _rules.crbegin(); it != _rules.crend(); ++it) {
        if (it->node[name] && std::regex_match(rule_match, it->regex)) {
            return std::regex_replace(rule_match,
                                      it->regex,
                                      it->node[name].as<std::string>());
        }
    }
    return default_value;
}

std::vector<std::string> style_list::get_formatted_all(const std::string &name,
                                                       const std::string &rule_match) const
{
    std::vector<std::string> vec;
    for (auto it = _rules.cbegin(); it != _rules.cend(); ++it) {
        if (it->node[name] && std::regex_match(rule_match, it->regex)) {
            auto repl = std::regex_replace(rule_match,
                                           it->regex,
                                           it->node[name].as<std::string>());
            vec.emplace_back(repl);
        }
    }
    return vec;
}
}

/// \cond
namespace YAML
{

template <> struct convert<util::style_list::rule>
{
    static bool decode(const Node &node, util::style_list::rule &rule)
    {
        if (!node["selector"]) {
            throw std::runtime_error("Style list item without selector");
        }

        rule.selector = node["selector"].as<std::string>();
        rule.regex = std::regex(rule.selector, std::regex_constants::extended);
        rule.node = node;

        return true;
    }
};
} // namespace YAML
/// \endcond
