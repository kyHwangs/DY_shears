#include <fstream>
#include <string>
#include <vector>

bool printTable(std::ostream &o,
                const std::vector<std::string> &colHeaders,
                const std::vector<std::string> &lineHeaders,
                const std::vector<std::vector<std::string>> &vals,
                const char *caption,
                const char *label,
                const char *format = "latex");

bool printTable(std::ostream &o,
                const std::vector<std::string> &colHeaders,
                const std::vector<std::string> &lineHeaders,
                const std::vector<std::vector<double>> &vals,
                const char *numFormat,
                const char *caption,
                const char *label = "",
                const char *format = "latex");
