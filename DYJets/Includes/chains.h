#ifndef CHAINS_H
#define CHAINS_H

#include <memory>
#include <string>
#include <vector>

class TChain;

namespace util
{

/**
 * \brief Utilitiy class to get \c TChain from Tupel-generated ROOT trees.
 */
class chains
{
    std::shared_ptr<TChain> _events;
    
  public:
    /// \brief Constructor.
    explicit chains(const std::vector<std::string> &files);

    /// \brief Returns a \c TChain pointing to event data.
    std::shared_ptr<TChain> events() { return _events; }

};
} // namespace util

#endif // CHAINS_H
