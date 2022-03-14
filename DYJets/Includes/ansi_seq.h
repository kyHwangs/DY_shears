#ifndef ANSI_SEQ_H
#define ANSI_SEQ_H

namespace util
{

/// \brief ANSI terminal control sequences
namespace ansi
{

/// \brief Represents a color
struct color
{
    const int code;
};

/// \brief Represents a color type (background and foreground)
struct color_type
{
    const int code;
};

const color_type foreground = color_type{30}; ///< \brief A foreground color.
const color_type background = color_type{40}; ///< \brief A background color.

const color black = color{0};   ///< \brief Black
const color red = color{1};     ///< \brief Red
const color green = color{2};   ///< \brief Green
const color yellow = color{3};  ///< \brief Yellow
const color blue = color{4};    ///< \brief Blue
const color magenta = color{5}; ///< \brief Magenta
const color cyan = color{6};    ///< \brief Cyan
const color white = color{7};   ///< \brief White

const color bright_red = color{61};     ///< \brief Bright red
const color bright_green = color{62};   ///< \brief Bright green
const color bright_yellow = color{63};  ///< \brief Bright yellow
const color bright_blue = color{64};    ///< \brief Bright blue
const color bright_magenta = color{65}; ///< \brief Bright magenta
const color bright_cyan = color{66};    ///< \brief Bright cyan
const color bright_white = color{67};   ///< \brief Bright white

/// Returns the sequence one has to print in order to use the given color.
inline std::string setcolor(color color, color_type type = foreground)
{
    return "\33[" + std::to_string(type.code + color.code) + "m";
}

/// Returns the sequence one has to print to get the terminal back to normal.
inline std::string reset() { return "\33[0m"; }

/// Returns the sequence one has to print to clear the current line.
inline std::string clear_line() { return "\033[K"; }
} // namespace ansi
} // namespace util

#endif // ANSI_SEQ_H
