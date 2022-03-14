#ifndef CONFIG_H
#define CONFIG_H

#include "TObject.h"
#include "TString.h"
#include "stdlib.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdio.h>
#include <string>
#include <utility>
#include <vector>

class Config
{
  public:
    Config(const char *filename = 0);

    void print(std::ostream &o = std::cout);

    bool read(const char *filename = 0);

    bool write(const char *filename);

    bool get(const char *key, bool &value);

    bool get(const char *key, int &value);

    bool get(const char *key, float &value);

    bool get(const char *key, double &value);

    bool get(const char *key, std::string &value);

    bool get(const char *key, std::vector<int> &value);

    bool get(const char *key, std::vector<float> &value);

    bool get(const char *key, std::vector<double> &value);

    bool get(const char *key, std::vector<std::string> &value);

    bool getB(const char *key, bool defaultValue = false);

    int getI(const char *key, int defaultValue = 0);

    int getL(const char *key, Long64_t defaultValue = 0);

    float getF(const char *key, float defaultValue = 0.);

    double getD(const char *key, double defaultValue = 0.);

    std::string getS(const char *key, std::string defaultValue = std::string(""));

    std::vector<bool> getVB(const char *key, std::vector<bool> defaultValue = std::vector<bool>());

    std::vector<int> getVI(const char *key, std::vector<int> defaultValue = std::vector<int>());

    std::vector<float> getVF(const char *key,
                             std::vector<float> defaultValue = std::vector<float>());

    std::vector<double> getVD(const char *key,
                              std::vector<double> defaultValue = std::vector<double>());

    std::vector<std::string>
    getVS(const char *key, std::vector<std::string> defaultValue = std::vector<std::string>());

    template <typename T> void set(const char *key, const T &value);

    void set(const char *key, const char *value)
    {
        std::string s(value);
        set(key, s);
    }

    void set(const char *key, const TString &value)
    {
        std::string s(value);
        set(key, s);
    }

  private:
    template <typename T> void convert(const T &x, std::string &s) const;

    template <typename T> void convert(const std::vector<T> &x, std::string &s) const;

    void convert(const std::string &str, float &x) const;

    void convert(const std::string &str, double &x) const;

    void convert(const std::string &str, bool &x) const;

    void convert(const std::string &str, int &x) const;

    void convert(const std::string &str, Long64_t &x) const;

    void convert(const std::string &str, std::string &x) const;

    template <typename T> void convert(const std::string &str, std::vector<T> &x) const;

    template <typename T> bool get(const char *key, T &value);

    template <typename T> T get2(const char *key, const T &defaultValue);

    std::string trim(std::string s, const char *chars = " \t") const;

    std::string tokenize(const std::string &s, const std::string &delim, int &pos) const;

    std::map<std::string, std::string> table_;
};

#endif // CONFIG_H not defined
