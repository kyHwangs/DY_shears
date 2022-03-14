\page coding-style Coding Style

This page describes several aspects that you should keep in mind while coding. *Please*
try to apply the recommandations.

## Indentation, braces and stuff

The coding style is defined in a `clang-format` configuration file located at the root of the
source tree. `clang-format` is a tool that can automatically modify any C++ file to use a given
coding style. It is used as:

~~~
clang-format -i <file>
~~~

Run the command above on any modified (C++) file before committing changes to git.

## Documentation

It is very important that you document your code. Add a comment starting with `///` or `/**` in
front of every function or class in a `.h` file. They will be picked by Doxygen and used to
generate a nice description of your code. Example (without syntax highlighting, sorry):

~~~
/**
 * \\brief Computes the foo of \c bar
 *  ^ Use only one \ above
 *
 * The foo of bar is defined as ...
 */
void foo(double bar);
~~~

A well-written documentation can be very, very useful to newcomers. It's also a good way to check
that you understand what your function is doing: if you can't describe it in one or two short
sentences, your function is probably too complicated.

Documentation checklist:

- Are all your functions documented?
- Did you use \\brief to provide a short description?
- Did you describe what every argument does?

I encourage that you educate yourself about Doxygen commands; your users will love you. You can
even use \f$\mathrm\LaTeX\f$ formulas in the documentation!

## Modern C++ features

The code requires C++ 14, and you are thus allowed to use any C++ 14 feature. This includes among
others the [range-based for loop](http://en.cppreference.com/w/cpp/language/range-for), the
[`auto`](http://en.cppreference.com/w/cpp/language/auto) and
[`nullptr`](http://en.cppreference.com/w/cpp/language/nullptr) keywords, and
[list initialization](http://en.cppreference.com/w/cpp/language/list_initialization). A special
mention should also go to [lambda functions](http://en.cppreference.com/w/cpp/language/lambda),
which are very useful when doing physics.

Modern C++ isn't limited to new keywords. It also advocates against the use of raw pointers, which
proved to be very error-prone (who doesn't know about segmentation violations?). Try to avoid
pointers in the first place, and use
[shared pointers](http://en.cppreference.com/w/cpp/memory/shared_ptr) when needed. You shouldn't
need raw pointers at all, except when interacting with ROOT or doing very low-level stuff.

C++ provides excellent support for many things. Since ROOT is more likely to be buggy/have unclear
interfaces/break than the Standard Library, prefer standard containers and strings over their ROOT
equivalents.

Did you know that C++ had built-in support for many common
[algorithms](http://en.cppreference.com/w/cpp/algorithm) such as sorting, for
[regular expressions](http://en.cppreference.com/w/cpp/regex), for
[numeric algorithms](http://en.cppreference.com/w/cpp/numeric#Numeric_algorithms) such as inner
product of vectors and for
[pseudo-random number generation](http://en.cppreference.com/w/cpp/numeric/random)? If the
Standard Library doesn't fulfill your needs, you might still find an existing algorithm in one of
the very well-designed and well-documented [Boost C++ Libraries](http://boost.org): there are even
numerical integrators for multi-dimensional differential equations!

## Error reporting

Errors can be divided in two categories: the ones that can safely be ignored and the ones that
can't. In general, the former should result in a warning, and the latter in an exception being
thrown. Example code:

~~~{.cc}
#include <stdexcept>
#include "logging.h"

void foo(double bar)
{
    if (bar == 0) {
        // Let's say this can be ignored without breaking the program
        util::logging::warn << "foo got bar = 0" << std::endl;
    } else if (bar < 0) {
        // This cannot happen for whatever reason
        throw std::range_error("foo got bar < 0");
    }
}
~~~
