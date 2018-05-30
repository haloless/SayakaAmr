#pragma once

#include <cassert>
#include <iosfwd>

#include "SayakaCommons.h"

SAYAKA_NS_BEGIN;

struct BlockFaceIndex
{
	enum {
		NxFace = 2,
		NyFace = NDIM>=2 ? 2 : 0,
		NzFace = NDIM==3 ? 2 : 0,

		NumFace = NDIM*2,
	};
	enum {
		FACE_XMINUS = 0,
		FACE_XPLUS, 
		FACE_YMINUS,
		FACE_YPLUS,
		FACE_ZMINUS,
		FACE_ZPLUS,
	};

	int m_index;

	BlockFaceIndex(int index=0)
		: m_index(index)
	{
		assert(0<=index && index<MAX_FACE);
	}
	BlockFaceIndex(int dir, int side)
		: m_index(dir*2+side)
	{
		assert(0<=dir && dir<MAX_DIM);
		assert(0<=side && side<2);
	}

	// conversion
	operator int() const { return m_index; }
	// allows ++
	operator int&() { return m_index; }

	int dir() const { return m_index/2; }
	int side() const { return m_index%2; }

	BlockFaceIndex opposite() const;
};

typedef BlockFaceIndex FaceIndex;

// IO
std::ostream& operator<<(std::ostream& os, const FaceIndex &iface);



SAYAKA_NS_END;


