#include "event_counter.h"

#include <iomanip>
#include <sstream>

#include "logging.h"

namespace util {

void event_counter::print() const
{
    // Find optimal column length
    std::size_t namelength = 0, countlength = 0, effcountlength = 0;
    for (const std::string &name : _names) {
        namelength = std::max(name.size(), namelength);

        const counter &c = _counts.at(name);
        countlength = std::max(std::to_string(c.count).size(), countlength);

        double effcount = c.effcount / c.effabscount * c.count;
        std::stringstream ss;
        ss << std::fixed << std::setprecision(0) << effcount;
        effcountlength = std::max(ss.str().size(), effcountlength);
    }
    // Table header labels
    const std::string namelabel = "Counter";
    namelength = std::max(namelabel.size(), namelength);
    const std::string countlabel = "Events";
    countlength = std::max(countlabel.size(), countlength);
    const std::string effcountlabel = "Eff.";
    effcountlength = std::max(effcountlabel.size(), effcountlength);
    // Print table header
    logging::info << std::setw(namelength + countlength + effcountlength + 4)
                  << std::setfill('=')
                  << ""
                  << std::setfill(' ')
                  << std::endl;
    logging::info << std::setw(namelength)
                  << std::setiosflags(std::ios_base::left)
                  << namelabel
                  << std::resetiosflags(std::ios_base::left)
                  << "  "
                  << std::setw(countlength)
                  << countlabel
                  << "  "
                  << std::setw(effcountlength)
                  << effcountlabel
                  << std::endl;
    logging::info << std::setw(namelength + countlength + effcountlength + 4)
                  << std::setfill('-')
                  << ""
                  << std::setfill(' ')
                  << std::endl;
    // Print table contents
    for (const std::string &name : _names) {
        const counter &c = _counts.at(name);
        logging::info << std::setw(namelength)
                      << std::setiosflags(std::ios_base::left)
                      << name
                      << "  "
                      << std::resetiosflags(std::ios_base::left)
                      << std::setiosflags(std::ios_base::right)
                      << std::setw(countlength)
                      << c.count
                      << "  "
                      << std::setw(effcountlength)
                      << std::fixed
                      << std::setprecision(0)
                      << (c.effcount / c.effabscount * c.count)
                      << std::resetiosflags(std::ios_base::right)
                      << std::endl;
    }
    // Table bottom
    logging::info << std::setw(namelength + countlength + effcountlength + 4)
                  << std::setfill('=')
                  << ""
                  << std::setfill(' ')
                  << std::endl;
}

void event_counter::declare(const std::string &name)
{
    _names.push_back(name);
    _counts[name] = {0, 0, 0};
}

void event_counter::count(const std::string &name, double weight)
{
    if (_counts.count(name) == 0) {
        throw std::logic_error("Counter not declared: " + name);
    }
    counter &c = _counts[name];
    c.count++;
    c.effabscount += std::abs(weight);
    c.effcount += weight;
}
} // namespace util
