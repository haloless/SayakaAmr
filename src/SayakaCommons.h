#pragma once


#ifndef SAYAKA_SPACEDIM
#	error Must define SAYAKA_SPACEDIM
#endif

#if (SAYAKA_SPACEDIM==2)
#elif (SAYAKA_SPACEDIM==3)
#else
#	error Must have SAYAKA_SPACEDIM=2 or 3
#endif

#define SAYAKA_NS_BEGIN namespace sayaka {
#define SAYAKA_NS_END }

#include <cassert>

#include <tuple>

SAYAKA_NS_BEGIN;

// Spacial dimensions
const int MAX_DIM = 3;
//
const int NDIM = SAYAKA_SPACEDIM;
const int XDIM = 1;
const int YDIM = NDIM>=2 ? 1 : 0;
const int ZDIM = NDIM==3 ? 1 : 0;

//
const int MAX_FACE = 6;
const int NUM_FACE = NDIM*2;

// AMR tree child block
const int MAX_CHILD = 8;
const int NUM_CHILD = (XDIM+1)*(YDIM+1)*(ZDIM+1);


// 
inline void SelectDirIncr(int dir, int &ii, int &jj, int &kk) {
	assert(0<=dir && dir<NDIM);
	ii = dir==0 ? 1 : 0;
	jj = dir==1 ? 1 : 0;
	kk = dir==2 ? 1 : 0;
}
inline void SelectStaggerIncr(int dir, int &ii, int &jj, int &kk) {
	SelectDirIncr(dir, ii, jj, kk);
}

inline auto SelectStaggerIncr(int dir) {
	int ii, jj, kk;
	SelectStaggerIncr(dir, ii, jj, kk);
	return std::make_tuple(ii, jj, kk);
}

inline int SelectDirIndex(int dir, int i, int j, int k) {
	if (dir == 0) { 
		return i;
	} else if (dir == 1) {
		return j;
	} else {
		return k;
	}
}


SAYAKA_NS_END;



