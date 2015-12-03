#pragma once

#include <cassert>

#include "config.h"

#include "SayakaCommons.h"

namespace sayaka {

//
const int verbose = 1;
inline int IsVerbose() { return verbose; }



/*
 *
 */
void Initialize(int argc, char *argv[]);

CONFIG& GetConfig(const char *name=NULL);






} // namespace sayaka


