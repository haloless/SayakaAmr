
#include <cstdio>
#include <cstdlib>

#include <vector>
#include <algorithm>

#include "log.h"

#include "Sayaka.h"



SAYAKA_NS_BEGIN;

const int verbose = 1;

int IsVerbose() { return verbose; }


// 
void Initialize(int argc, char *argv[]) {
	LOGPRINTF("%s\n", __FUNCTION__);
	return;
}


SAYAKA_NS_END;
