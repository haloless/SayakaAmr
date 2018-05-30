#pragma once

#include <cassert>

#include "config.h"

#include "SayakaCommons.h"

SAYAKA_NS_BEGIN;


//
int IsVerbose();



/*
 *
 */
void Initialize(int argc, char *argv[]);

CONFIG& GetConfig(const char *name=NULL);






SAYAKA_NS_END;

