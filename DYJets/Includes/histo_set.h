#ifndef HISTO_SET_H
#define HISTO_SET_H

#include <map>
#include <stdexcept>

#include <TH1D.h>

#include "GenH1D.h"
#include "style_list.h"

namespace util
{

/**
 * \brief Manages a collection of histograms.
 *
 * This class can be used to manage histograms, simplifying the declare-fill-write paradigm. It is
 * used as follows:
 *
 * 1. \ref declare "Declare" your histograms as if you were creating them
 * 1. \ref fill "Fill" them as usual
 * 1. \ref write "Write" them all at once
 *
 * This class also features *tags* that allow one to use the same declaration for similar
 * histograms (ie when the selection is different).
 *
 * Histogram definitions can also come from a stylesheet set using \ref set_style. No call to
 * \ref declare is needed in this case. Information in the stylesheet take precedence over calls to
 * \c declare.
 */
class histo_set
{
  public:
    using histogram_type = GenH1D;

  private:
    std::map<std::string, histogram_type> _models;
    style_list _style;
    std::map<std::pair<std::string, std::string>, histogram_type> _histograms;

  public:
    /// \brief Constructor.
    explicit histo_set() = default;

    /// \brief Destructor.
    virtual ~histo_set() = default;

    /**
     * \brief Declares a new histogram.
     *
     * Histograms have to be declared before being filled. Histogram names must be unique.
     *
     * \param name The name of the histogram, forwarded to the histogram constructor as first
     *             argument
     * \param args Further arguments to forward to the histogram constructor.
     */
    template <class... Args> void declare(const std::string &name, const Args &... args)
    {
        if (_models.count(name) > 0) {
            throw std::logic_error("Histogram declared twice: " + name);
        }
        _models.emplace(name, histogram_type(name.c_str(), args...));
    }

    /**
     * \brief Fills an histogram.
     *
     * Fills an histogram named \c name, previously \ref declare "declared". The \c tag can be used
     * to distinguish between different version of the same histogram (eg using different
     * selections). Other arguments are forwarded to the \c Fill function of the histogram.
     */
    template <class... Args>
    void fill(const std::string &name, const std::string &tag, const Args &... args)
    {
        get(name, tag).Fill(args...);
    }

    /// \brief Sets the stylesheet.
    void set_style(const style_list &style) { _style = style; }

    /**
     * \brief Writes all histograms in the current (\c ROOT) directory.
     *
     * \note Empty histograms are never written.
     */
    void write();

    /**
     * \brief Retrieves a string that combines the given name and tag.
     * \returns \c name if \c tag is empty, \c name_tag otherwise.
     */
    static std::string combined_name(const std::string &name, const std::string &tag);

    /**
     * \brief Retrieves the histogram with the given \c name and \c tag.
     *
     * The histogram is created if it doesn't exist.
     */
    histogram_type &get(const std::string &name, const std::string &tag = "");

  private:
    /**
     * \brief Creates an histogram from the style.
     *
     * \returns An iterator pointing to the end of \c _histograms if the histogram couldn't be
     *          created, an iterator pointing to its position in \c _histograms otherwise.
     */
    auto create_from_style(const std::string &name, const std::string &tag)
        -> decltype(_histograms)::iterator;
};
} // namespace util

#endif // HISTO_SET_H
