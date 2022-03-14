#include <map>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
#include <vector>
#include "stdlib.h"
#include "TObject.h"
#include <sstream>
#include "Config.h"

Config::Config(const char* filename){
    read(filename);
}

void Config::print(std::ostream& o){
    o << "#### vjets config file automatically generated on the ";
    o << __DATE__ << " at " << __TIME__ << std::endl;
    for(std::map<std::string, std::string>::const_iterator it = table_.begin();
            it != table_.end();
            ++it){
        o << it->first << " = " << it->second << std::endl;
    }
}

bool Config::read(const char* filename){
    if(filename == 0) filename = getenv("VJETS_CONFIG");
    if(filename == 0) filename = "vjets.cfg";

    std::ifstream f(filename);
    if(!f.good()){
        std::cerr << "Failed to read file " << filename << "!\n";
        return false;
    }

    std::cerr << "\n\tReading configuration from file " << filename << "\n";
    std::cerr << "--------------------------------------------------\n";

    std::string line;
    while(getline(f, line).good()){
        trim(line);
        if(line.size() == 0 || line[0] == '#') continue;
	int pos=0;
        std::string key = tokenize(line, " \t=", pos);
        std::string value = tokenize(line, "\n", pos);
	bool escaped = false;
	for(pos = 0; (unsigned) pos < value.size() && value[pos]!='#'; ++pos){
	  if(value[pos] == '\\' && !escaped){
	    value.erase(pos, 1); 
	    escaped = true;
	  } else{
	    escaped = false;
	  }
	}
	if(pos < (int)value.size()) value.erase(pos);
        value = trim(value, " \t=");
        std::cerr.width(16);
        std::cerr << std::left << key + ": " << value << std::endl;
        table_[key] = value;
    }
    std::cerr << "--------------------------------------------------\n" << std::endl;
    return true;
}

bool Config::write(const char* filename){
    std::ofstream f(filename);
    if(!f.good()){
        std::cerr << "Failed to read file " << filename << "!\n";
        return false;
    }
    print(f);
    return true;
}

bool Config::get(const char* key, bool& value){
    return get<bool>(key, value);
}

bool Config::get(const char* key, int& value){
    return get<int>(key, value);
}

bool Config::get(const char* key, float& value){
    return get<float>(key, value);
}

bool Config::get(const char* key, double& value){
    return get<double>(key, value);
}

bool Config::get(const char* key, std::string& value){
    return get<std::string>(key, value);
}

bool Config::get(const char* key, std::vector<int>& value){
    return get<std::vector<int> >(key, value);
}

bool Config::get(const char* key, std::vector<float>& value){
    return get<std::vector<float> >(key, value);
}

bool Config::get(const char* key, std::vector<double>& value){
    return get<std::vector<double> >(key, value);
}

bool Config::get(const char* key, std::vector<std::string>& value){
    return get<std::vector<std::string> >(key, value);
}

bool Config::getB(const char* key, bool defaultValue){
    return get2<bool>(key, defaultValue);
}

int Config::getI(const char* key, int defaultValue){
    return get2<int>(key, defaultValue);
}

int Config::getL(const char* key, Long64_t defaultValue){
    return get2<Long64_t>(key, defaultValue);
}

float Config::getF(const char* key, float defaultValue){
    return get2<float>(key, defaultValue);
}

double Config::getD(const char* key, double defaultValue){
    return get2<double>(key, defaultValue);
}

std::string Config::getS(const char* key, std::string defaultValue){
    return get2<std::string>(key, defaultValue);
}

std::vector<bool> Config::getVB(const char* key, std::vector<bool> defaultValue){
    return get2<std::vector<bool> >(key, defaultValue);
}

std::vector<int> Config::getVI(const char* key, std::vector<int> defaultValue){
    return get2<std::vector<int> >(key, defaultValue);
}

std::vector<float> Config::getVF(const char* key, std::vector<float> defaultValue){
    return get2<std::vector<float> >(key, defaultValue);
}

std::vector<double> Config::getVD(const char* key, std::vector<double> defaultValue){
    return get2<std::vector<double> >(key, defaultValue);
}

std::vector<std::string> Config::getVS(const char* key, std::vector<std::string> defaultValue){
    return get2<std::vector<std::string> >(key, defaultValue);
}


template<typename T>
void Config::convert(const T& x, std::string& s) const{
    std::stringstream ss;
    ss << x;
    s = ss.str();
}

template<typename T>
void Config::convert(const std::vector<T>& x, std::string& s) const{
    s = "";
    for(size_t i = 0; i < x.size(); ++i){
        if(i) s += ", ";
        std::string e;
        convert(x[i], e);
        s += e;
    }
}

void Config::convert(const std::string& str, float& x) const{
    x = strtod(str.c_str(), 0);
}

void Config::convert(const std::string& str, double& x) const{
    x = strtod(str.c_str(), 0);
}

void Config::convert(const std::string& str, bool& x) const{
  x = strtol(str.c_str(), 0, 0); //handles 0, 1. It will give also the corret value (0) for f(alse) and n(o)
  if(str.size()>0
     && (str[0]=='t' || str[0] =='T'
	 || str[0]=='y' || str[0] =='Y')) x = true;
}

void Config::convert(const std::string& str, int& x) const{
    x = strtol(str.c_str(), 0, 0);
}

void Config::convert(const std::string& str, Long64_t& x) const{
    x = strtoll(str.c_str(), 0, 0);
}

void Config::convert(const std::string& str, std::string&  x) const{  
    x = str;
}

template<typename T>
void Config::convert(const std::string& str, std::vector<T>& x) const{
    x.clear();
    int pos = 0;
    std::string s = str;
    trim(s, "{}");
    while(pos >= 0){
        std::string val = tokenize(s, ", ", pos);
        if(pos == 0) break;
        trim(val);
        T val_;
        convert(val, val_);
        // std::cerr << "Adding " << val_ << ", pos = " << pos << std::endl;
        x.push_back(val_);
    }
}

std::string Config::trim(std::string s, const char* chars) const{
    std::string::size_type pos0 = s.find_first_not_of(chars);
    if(pos0==std::string::npos){
        pos0=0;
    }
    std::string::size_type pos1 = s.find_last_not_of(chars) + 1;
    if(pos1==std::string::npos){
        pos1 = pos0;
    }
    return s.substr(pos0, pos1-pos0);
}

std::string Config::tokenize(const std::string& s, const std::string& delim,
			     int& pos) const{
    if(pos<0) return "";
    int pos0 = pos;
    int len = s.size();
    //eats delimeters at beginning of the string
    while(pos0<len && find(delim.begin(), delim.end(), s[pos0])!=delim.end()){
        ++pos0;
    }
    if(pos0==len) return "";
    pos = s.find_first_of(delim, pos0);
    return s.substr(pos0, (pos>0?pos:len)-pos0);
}


template<typename T>
bool Config::get(const char* key, T& value){
    std::map<std::string, std::string>::const_iterator it = table_.find(key);
    if(it!=table_.end()){
        convert(it->second, value);
        return true;
    } else{
        std::cerr << "Warning: parameter " << key << " was not found in configuration file!\n";
        std::string value_;
        convert(value, value_);
        std::cerr << "         will use " << value_ << " as default value or the command line option.\n"; 
        std::string key_(key);
        table_[key_] = value_;
        return false;
    }
}

template<typename T>
T Config::get2(const char* key, const T& defaultValue){
    T value = defaultValue;
    get(key, value);
    return value;
}

template<typename T>
void Config::set(const char* key, const T& value){
    std::string sVal;
    convert(value, sVal);    
    T prevValue;
    std::map<std::string, std::string>::iterator it = table_.find(key);
    if(it!=table_.end()){
      convert(it->second, prevValue);
      if(prevValue != value){
            std::cerr << "Info: the value of " << key << " from configuration has been overwritten "
                                                         "from "
                      << it->second << " to " << sVal << "\n";
            it->second = sVal;
        }
    } else {
        table_[key] = sVal;
    }
}

template
void Config::set<bool>(const char* key, const bool& value);
template
void Config::set<int>(const char* key, const int& value);
template
void Config::set<Long64_t>(const char* key, const Long64_t& value);
template
void Config::set<float>(const char* key, const float& value);
template
void Config::set<double>(const char* key, const double& value);
template
void Config::set<std::string>(const char* key, const std::string& value);
template
void Config::set<std::vector<int> >(const char* key, const std::vector<int>& value);
template
void Config::set<std::vector<float> >(const char* key, const std::vector<float>& value);
template
void Config::set<std::vector<double> >(const char* key, const std::vector<double>& value);
template
void Config::set<std::vector<std::string> >(const char* key, const std::vector<std::string>& value);

