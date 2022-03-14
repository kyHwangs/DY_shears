#ifndef TABLES_H
#define TABLES_H

#include <map>
#include <string>
#include <vector>

namespace util
{

/**
 * \brief Handles table of scale factors.
 *
 * The class reads a table from the filesystem, and provides access to its entries.
 *
 * YAML string nodes can be converted to tables.
 */
class table
{
    /// \brief A bin
    struct record
    {
        double ptLow, ptHi, etaLow, etaHi, effi, effiErrorLow, effiErrorHigh;

        bool belongToEta(double) const;
        bool belongTo(double, double) const;
        bool equalTo(int num) const;
    };

    std::vector<record> _recd;

  public:
    /// \brief Creates an invalid table.
    explicit table() = default;

    /**
     * \brief Reads a table from the file provided in argument.
     * \throws std::invalid_argument if the file doesn't exist
     * \warning The reader isn't very robust. Be careful when adding tables!
     */
    explicit table(const std::string &filename);

    /// \brief Returns the efficiency in the given (\c pt, \c eta) bin.
    double getEfficiency(double pt, double eta) const;

    /// \brief Returns the lower bound for the efficiency in the given (\c pt, \c eta) bin.
    double getEfficiencyLow(double pt, double eta) const;

    /// \brief Returns the upper bound for the efficiency in the given (\c pt, \c eta) bin.
    double getEfficiencyHigh(double pt, double eta) const;
};

/// \brief A set of tables identified by their name.
using tables = std::map<std::string, table>;
} // namespace util

#ifdef DYJETS_NEW_API

/// \cond
namespace YAML
{

class Node;
template<class T> struct convert;

template <> struct convert<util::table>
{
    static bool decode(const Node &node, util::table &table);
};
} // namespace YAML
/// \endcond

#endif // DYJETS_NEW_API

#endif // TABLES_H
