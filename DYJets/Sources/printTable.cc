#include "printTable.h"
#include "TString.h"
#include <algorithm>
#include <iostream>

bool printTable(std::ostream &o,
                const std::vector<std::string> &colHeaders,
                const std::vector<std::string> &lineHeaders,
                std::vector<std::vector<double>> &vals,
                const char *numFormat,
                const char *caption,
                const char *label,
                const char *format = "latex")
{
    std::vector<std::vector<std::string>> svals(vals.size());
    for (unsigned irow = 0; irow < vals.size(); ++irow) {
        for (unsigned icol = 0; icol < vals[irow].size(); ++icol) {
            svals[irow].push_back(TString::Format(numFormat, vals[irow][icol]).Data());
        }
    }
    return printTable(o, colHeaders, lineHeaders, svals, caption, label, format);
}

bool printTable(std::ostream &o,
                const std::vector<std::string> &colHeaders,
                const std::vector<std::string> &lineHeaders,
                const std::vector<std::vector<std::string>> &vals,
                const char *caption,
                const char *label,
                const char *format)
{

    unsigned ncols;
    if (vals.size() > 0) {
        ncols = vals[0].size();
        if (lineHeaders.size()) ++ncols;
    } else {
        ncols = colHeaders.size();
    }

    if (colHeaders.size() > 0 && colHeaders.size() != ncols) {
        std::cerr << "Fatal error. " << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ << ": "
                  << "Header and lines have different numbers of columns.\n";
        std::cerr << "Number of column headers: " << colHeaders.size()
                  << "\tNumber of columns: " << ncols << "\n";

        abort();
    }

    for (unsigned i = 1; i < vals.size(); ++i) {
        if (vals[i].size() != vals[0].size()) {
            std::cerr << "Fatal error. " << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__
                      << ": "
                      << "All the rows do not contain the same number of columns.\n";
            abort();
        }
    }

    if (lineHeaders.size() > 0 && lineHeaders.size() != vals.size()) {
        std::cerr << "Fatal error. " << __FUNCTION__ << ":" << __FILE__ << ":" << __LINE__ << ": "
                  << "Invalid size of line header vector.\n";
        abort();
    }

    std::string colFormat;
    if (lineHeaders.size() > 0) {
        colFormat = "l|";
    }
    for (unsigned i = 1; i < ncols; ++i) colFormat += "c";

    o << "\\begin{table}\n"
         "\\begin{center}\n";
    if (caption && caption[0]) {
        o << "\\caption{" << caption << "}\n";
    }
    o << "\\resizebox{\\textwidth}{!}{%\n"
         "\\begin{tabular}{"
      << colFormat << "}\n";

    std::vector<unsigned> colWidth(ncols, 1);
    for (unsigned irow = 0; irow < vals.size(); ++irow) {
        int nHeader = 0;
        if (lineHeaders.size() > 0) {
            nHeader = 1;
            unsigned w = lineHeaders[irow].size();
            if (colWidth[0] < w) colWidth[0] = w;
        }
        for (unsigned icol = nHeader; icol < ncols; ++icol) {
            unsigned w = vals[irow][icol - nHeader].size();
            if (colWidth[icol] < w) colWidth[icol] = w;
        }
    }

    int lineLen = 0;
    int maxLineLen = 120;
    if (colHeaders.size()) {
        const char *sep = "";
        for (unsigned icol = 0; icol < ncols; ++icol) {
            int newLineLen = lineLen + colHeaders[icol].size() + strlen(sep);
            if (lineLen > 0 && lineLen > maxLineLen) {
                o << "\n";
                newLineLen = 0;
            }
            lineLen = newLineLen;
            o << sep << colHeaders[icol];
            sep = " & ";
        }
        o << "\\\\\n\\hline\n";
    }

    for (unsigned irow = 0; irow < vals.size(); ++irow) {
        const char *sep = "";
        //    if(lineHeaders.size() > 0){
        //      for(unsigned i = lineHeaders[0].size(); i < colWidth[0]; ++i) o << " ";
        //    }

        lineLen = 0;
        for (unsigned icol = 0; icol < ncols; ++icol) {
            int ival = icol;
            TString text;
            int newLineLen = lineLen + colWidth[icol];
            if (lineLen > 0 && newLineLen > maxLineLen) {
                o << "\n";
                newLineLen = 0;
            }
            lineLen = newLineLen;
            if (lineHeaders.size() > 0) {
                ival -= 1;
            }
            if (ival >= 0)
                text = vals[irow][ival];
            else
                text = lineHeaders[irow];
            // space padding:
            o << sep << text;
            for (unsigned i = text.Length(); i < colWidth[icol]; ++i) o << " ";
            sep = " & ";
        }
        o << "\\\\\n";
    }
    o << "\\end{tabular}}\n";
    if (label || strlen(label) > 0) o << "\\label{" << label << "}\n";
    o << "\\end{center}\n\\end{table}\n";

    return true;
}
