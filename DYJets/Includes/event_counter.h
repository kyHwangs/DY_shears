#ifndef EVENT_COUNTER_H
#define EVENT_COUNTER_H

#include <map>
#include <string>
#include <vector>

namespace util
{

/**
 * \brief Counts events
 *
 * This class can be used to count the number of events passing certain criteria. It is used in a
 * way similar to \ref histo_set.
 */
class event_counter final
{
    struct counter
    {
        long long count;
        double effabscount, effcount;
    };

    std::vector<std::string> _names;
    std::map<std::string, counter> _counts;

  public:
    /// \brief Constructor
    explicit event_counter() = default;

    /// \brief Declares a new counter
    void declare(const std::string &name);

    /**
     * \brief Increment the counter given by \c name.
     * \throws std::logic_error if the counter wasn't declared
     */
    void count(const std::string &name, double weight);

    /**
     * \brief Prints a summary table of declared counters.
     */
    void print() const;
};
} // namespace util

#endif // EVENT_COUNTER_H
