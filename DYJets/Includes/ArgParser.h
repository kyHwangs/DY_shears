#ifndef _ARGPARSER_H_
#define _ARGPARSER_H_

#include <TString.h>
#include <iostream>

void getArg(TString fullArg, TString &arg);
void getArg(TString fullArg, int &arg);
void getArg(TString fullArg, long &arg);
void getArg(TString fullArg, double &arg);
void getArg(TString fullArg, bool &arg);

#endif
