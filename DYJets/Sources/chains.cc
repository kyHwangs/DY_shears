#include "chains.h"

#include <stdexcept>
#include <string>

#include <TChain.h>
#include <TFile.h>

#include "logging.h"

namespace util
{

chains::chains(const std::vector<std::string> &files)
    : _events(std::make_shared<TChain>())
{

    for (const std::string &fullpath : files) {
        logging::debug << "Adding " << fullpath << " to the list of input files." << std::endl;
        std::string nanotreePath = fullpath + "/Events";
        

                _events->Add(nanotreePath.c_str());
}
}
} // namespace util
