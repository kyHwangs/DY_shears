#ifndef CATALOG_H
#define CATALOG_H

#include <limits>
#include <memory>
#include <string>
#include <vector>

class TChain;

/// \brief Classes describing data and datasets
namespace data
{

/**
 * \brief Represents the contents of a catalog.
 */
class catalog
{
    std::shared_ptr<TChain> _event_chain;
    bool _chains_initialized;

    std::vector<std::string> _files;
    double _lumi, _xsec;
    long long _primary_events = -1;
    bool _isdata = false;

  public:
    /**
     * \brief Constructor.
     *
     * The constructor reads the catalog. If `filename` is relative, it is
     * interpreted as being relative to `bonzaiDir`. Paths beginning with
     * "/store/" are interpreted as pointing to an `eos`-based directory.
     *
     * \param filename The location of the catalog.
     * \param bonzaiDir The `bonzaiDir` parameter from the config file, used
     *                  when paths in the catalog aren't absolute (Run I catalog
     *                  support).
     * \param maxFiles The maximum number of files to be used (-1 for no limit).
     */
    explicit catalog(const std::string &filename,
                     const std::string &nanoDir,
                     std::size_t maxFiles = std::numeric_limits<int>::max());

    /// \brief Destructor
    virtual ~catalog();

    /// \brief Returns a \c TChain pointing to event data.
    [[deprecated]] std::shared_ptr<TChain> event_chain()
    {
        if (!_chains_initialized) {
            initialize_chains();
        }
        return _event_chain;
    }


    /// \brief Returns the list of files read from the catalog.
    std::vector<std::string> files() const { return _files; }

    /// \brief Returns the integrated luminosity read from the catalog.
    double lumi() const { return _lumi; }

    /// \brief Returns the integrated cross section read from the catalog.
    double xsec() const { return _xsec; }

    /// \brief Returns true if the catalog contains data.
    bool isdata() const { return _isdata; }

    /// \brief Returns true if the catalog contains an MC dataset.
    bool ismc() const { return !_isdata; }

    /// \brief Returns the number of primary events in the dataset.
    long long primary_events() const { return _primary_events; }
  private:
    void initialize_chains();
};
} // namespace data

#endif // CATALOG_H
