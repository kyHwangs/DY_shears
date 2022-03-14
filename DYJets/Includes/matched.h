#ifndef MATCHED_H
#define MATCHED_H

#include <boost/optional.hpp>

namespace util
{

/**
 * \brief Groups gen/rec variables together.
 */
template <class T> class matched
{
  public:
    using type_t = T;

    boost::optional<type_t> gen; ///< \brief The gen-level value, if any
    boost::optional<type_t> rec; ///< \brief The reco-level value, if any

    /**
     * \brief Compares two \c matched objects for equality
     */
    bool operator==(const matched<type_t> &other) const
    {
        return gen == other.gen && rec == gen.rec;
    }

    /**
     * \brief Returns a \c matched object with the result of applying \c f on
     *        \ref gen and \ref rec.
     *
     * The generated function calls are:
     *
     *     f(*gen, args...)
     *     f(*rec, args...)
     *
     * If \c gen and/or \c rec is empty, the corresponding value(s) in the
     * returned will also be empty.
     */
    template <class Functor, class... Args>
    auto apply(Functor &&f, Args &&... args)
        -> matched<typename std::decay<decltype(f(*gen, std::forward<Args>(args)...))>::type>
    {
        if (gen && rec) {
            return {f(*gen, std::forward<Args>(args)...), f(*rec, std::forward<Args>(args)...)};
        } else if (gen) {
            return {f(*gen, std::forward<Args>(args)...), boost::none};
        } else {
            return {boost::none, f(*rec, std::forward<Args>(args)...)};
        }
    }

    /**
     * \brief Returns a \c matched object with the result of applying \c f on
     *        \ref gen and \ref rec.
     *
     * The generated function calls are:
     *
     *     f(*gen, args...)
     *     f(*rec, args...)
     *
     * If \c gen and/or \c rec is empty, the corresponding value(s) in the
     * returned will also be empty.
     */
    template <class Functor, class... Args>
    auto apply(Functor &&f, Args &&... args) const
        -> matched<typename std::decay<decltype(f(*gen, std::forward<Args>(args)...))>::type>
    {
        if (gen && rec) {
            return {f(*gen, std::forward<Args>(args)...), f(*rec, std::forward<Args>(args)...)};
        } else if (gen) {
            return {f(*gen, std::forward<Args>(args)...), boost::none};
        } else {
            return {boost::none, f(*rec, std::forward<Args>(args)...)};
        }
    }

    /**
     * \brief Returns a \c matched object with the result of applying the member
     *        function \c f on \ref gen and \ref rec.
     *
     * The generated function calls are:
     *
     *     (*gen.*f)(args...)
     *     (*rec.*f)(args...)
     *
     * If \c gen and/or \c rec is empty, the corresponding value(s) in the
     * returned will also be empty.
     */
    template <class MemberFunctor, class... Args>
    auto apply(MemberFunctor &&f, Args &&... args)
        -> matched<typename std::decay<decltype((*gen.*f)(std::forward<Args>(args)...))>::type>
    {
        if (gen && rec) {
            return {(*gen.*f)(std::forward<Args>(args)...), (*rec.*f)(std::forward<Args>(args)...)};
        } else if (gen) {
            return {(*gen.*f)(std::forward<Args>(args)...), boost::none};
        } else {
            return {boost::none, (*rec.*f)(std::forward<Args>(args)...)};
        }
    }

    /**
     * \brief Returns a \c matched object with the result of applying the member
     *        function \c f on \ref gen and \ref rec.
     *
     * The generated function calls are:
     *
     *     (*gen.*f)(args...)
     *     (*rec.*f)(args...)
     *
     * If \c gen and/or \c rec is empty, the corresponding value(s) in the
     * returned will also be empty.
     */
    template <class MemberFunctor, class... Args>
    auto apply(MemberFunctor &&f, Args &&... args) const
        -> matched<typename std::decay<decltype((*gen.*f)(std::forward<Args>(args)...))>::type>
    {
        if (gen && rec) {
            return {(*gen.*f)(std::forward<Args>(args)...), (*rec.*f)(std::forward<Args>(args)...)};
        } else if (gen) {
            return {(*gen.*f)(std::forward<Args>(args)...), boost::none};
        } else {
            return {boost::none, (*rec.*f)(std::forward<Args>(args)...)};
        }
    }
};

} // namespace util

#endif // MATCHED_H
