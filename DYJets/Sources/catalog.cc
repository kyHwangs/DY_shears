#include "catalog.h"

#include <iostream>
#include <regex.h>

#include <boost/algorithm/string.hpp>

#include <TChain.h>
#include <TFile.h>

#include "functions.h"
#include "logging.h"

namespace logging = util::logging;

namespace data
{

catalog::catalog(const std::string &filename, const std::string &bonzaiDir, std::size_t maxFiles)
    : _event_chain(std::make_shared<TChain>()),
      _chains_initialized(false),
      _lumi(0),
      _xsec(0)
{
    std::string fullpath = filename;
    if (filename[0] == '/' || boost::contains(filename, "://")) { // Absolute path
        fullpath = filename;
    } else {
        fullpath = bonzaiDir + "/" + filename;
    }

    if (boost::starts_with(fullpath, "/store/")) {
        fullpath.insert(0, "root://eoscms.cern.ch//eos/cms");
    }

    regex_t xsecLine;
    int rc = regcomp(
        &xsecLine, "[#*][[:space:]]*sample xsec[[:space:]:=]\\+\\([[:digit:].eE+-]\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &xsecLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        logging::fatal << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer
                       << std::endl;
    }

    regex_t lumiLine;
    rc = regcomp(&lumiLine, "[#*][[:space:]]*lumi[[:space:]:=]\\+\\([[:digit:].eE+-]\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &lumiLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        logging::fatal << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer
                       << std::endl;
    }

    regex_t primEvtLine;
    rc = regcomp(
        &primEvtLine, "[#*][[:space:]]*primary events[[:space:]:=]\\+\\([[:digit:]]\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &primEvtLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        logging::fatal << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer
                       << std::endl;
    }

    regex_t dataTypeLine;
    rc = regcomp(&dataTypeLine, "[#*][[:space:]]*data type[[:space:]:=]\\+\\(.*\\+\\)", 0);
    if (rc) {
        char buffer[256];
        regerror(rc, &dataTypeLine, buffer, sizeof(buffer));
        buffer[sizeof(buffer) - 1] = 0;
        logging::fatal << "Bug found in " << __FILE__ << ":" << __LINE__ << ": " << buffer
                       << std::endl;
    }

    if (isRootFile(fullpath.c_str())) {
        logging::info << "Loading file: " << fullpath << endl;
        _files.push_back(fullpath);
    } else {
        int (*closeFunc)(FILE *);
        FILE *f = eosOpen(fullpath.c_str(), &closeFunc);
        if (!f) {
            logging::error << "Failed to open file: " << fullpath << std::endl;
        } else {
            logging::info << "Reading catalog: " << fullpath << std::endl;

            string line;
            char *buffer = 0;
            size_t buffer_size = 0;
            // if maxFiles = 0 only catalog header is read.
            while (!feof(f)) {
                ssize_t len = getline(&buffer, &buffer_size, f);
                if (len < 0) break;
                char *line = buffer;
                size_t n = len - 1;
                // trim white spaces:
                while (line[0] == ' ' || line[0] == '\t') {
                    ++line;
                    --n;
                }
                while (line[n] == ' ' || line[n] == '\t' || line[n] == '\r' || line[n] == '\n') {
                    line[n] = 0;
                    --n;
                }

                regmatch_t pmatch[2];
                if (!regexec(&xsecLine, line, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)) {
                    line[pmatch[1].rm_eo] = 0;
                    _xsec = strtod(line + pmatch[1].rm_so, 0);
                    if (_xsec == 0) {
                        logging::error << "Value of parameter 'sample xsec', "
                                       << line + pmatch[1].rm_so << " found in file " << fullpath
                                       << " is not valid." << std::endl;
                    }
                } else if (!regexec(
                               &lumiLine, line, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)) {
                    line[pmatch[1].rm_eo] = 0;
                    _lumi = strtod(line + pmatch[1].rm_so, 0);
                    if (_lumi == 0) {
                        logging::error << "Integrated luminosity parameter value, "
                                       << line + pmatch[1].rm_so << " found in file " << fullpath
                                       << " is not valid." << std::endl;
                    }
                } else if (!regexec(
                               &primEvtLine, line, sizeof(pmatch) / sizeof(pmatch[0]), pmatch, 0)) {
                    line[pmatch[1].rm_eo] = 0;
                    _primary_events = strtoll(line + pmatch[1].rm_so, 0, 10);
                    if (_primary_events == 0) {
                        logging::error << "Primary events parameter value, "
                                       << line + pmatch[1].rm_so << " found in file " << fullpath
                                       << " is not valid." << std::endl;
                    }
                } else if (!regexec(&dataTypeLine,
                                    line,
                                    sizeof(pmatch) / sizeof(pmatch[0]),
                                    pmatch,
                                    0)) {
                    line[pmatch[1].rm_eo] = 0;
                    std::string type = line + pmatch[1].rm_so;
                    if (type == "data") {
                        _isdata = true;
                    } else if (type == "mc") {
                        _isdata = false;
                    } else {
                        logging::error << "Data type parameter value, " << type << " found in file "
                                       << fullpath << " is not valid." << std::endl;
                    }
                }

                // skip empty lines,  comment lines and metadata lines:
                if (line[0] == 0 || line[0] == '#' || line[0] == '*') continue;

                if (maxFiles == 0 || (_event_chain == 0 && _event_chain == 0)) break;

                // keep content of first column only:
                char *p = line;
                while (*p != 0 && *p != ' ' && *p != '\t' && *p != '\r' && *p != '\n') ++p;
                *p = 0;

                if (_files.size() >= maxFiles) break;

                TString filePath = TString(line);

                if (filePath[0] != '/' && !filePath.Contains("://")) {
                    filePath.Insert(0, TString(bonzaiDir) + "/");
                }
                if (filePath.BeginsWith("/store/")) {
                    filePath.Insert(0, "root://eoscms.cern.ch//eos/cms");
                }

                _files.push_back(filePath.Data());
            } // next line
            logging::debug << "Found " << _files.size() << " file(s) in catalog." << std::endl;
            if (buffer) free(buffer);
            closeFunc(f);
        } // file opening succeeded
    }     // is root file
    regfree(&xsecLine);
    regfree(&lumiLine);
}

catalog::~catalog() {}

void catalog::initialize_chains()
{
    logging::warn << "Using the old, slower method to read data." << std::endl;
    logging::info << "Initializing chains (this may take a while)..." << std::endl;

    enum { False = 0, True, Unknown } isBonzai = Unknown;

    for (const std::string &fullpath : _files) {
        logging::debug << fullpath << std::endl;
        std::string nanotreePath = fullpath + "Events";
        if (_event_chain) _event_chain->Add(nanotreePath.c_str());
            TFile *f = TFile::Open(fullpath.c_str());
    }            
    _chains_initialized = true;
    logging::info << "Chains initialized." << std::endl;
}
} // namespace data
