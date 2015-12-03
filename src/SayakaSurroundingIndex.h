#pragma once

#include "vector3d.h"

#include "SayakaCommons.h"
#include "SayakaFaceIndex.h"

namespace sayaka
{


/**
 * block 3x3x3
 */
struct SurroundingIndex
{
	enum {
		NxSurr = XDIM*2+1,
		NySurr = YDIM*2+1,
		NzSurr = ZDIM*2+1,
		NumSurr = NxSurr*NySurr*NzSurr,

		PosLow = -1,
		PosMid = 0,
		PosHigh = 1,

		CenterIndex = NumSurr/2,
	};

	//
	int m_index;

	SurroundingIndex(int index)
		: m_index(index)
	{
		assert(0<=index && index<NumSurr);
	}
	SurroundingIndex(int ipos, int jpos, int kpos)
	{
		setPos(ipos, jpos, kpos);
	}




	int ipos() const {
		return (m_index%3) - 1;
	}
	int jpos() const {
		if (NDIM >= 2) {
			return (m_index/3)%3 - 1;
		} else {
			return 0;
		}
	}
	int kpos() const {
		if (NDIM == 3) {
			return (m_index/9) - 1;
		} else {
			return 0;
		}
	}
	Vector3i pos() const {
		Vector3i vpos; 
		vpos.setZero();
		vpos.x = ipos();
		if (NDIM >= 2) vpos.y = jpos();
		if (NDIM == 3) vpos.z = kpos();
		return vpos;
	}
	void setPos(int ipos, int jpos, int kpos) {
		assert(-1<=ipos && ipos<=1);
		assert(-1<=jpos && jpos<=1);
		assert(-1<=kpos && kpos<=1);

		m_index = (ipos+1);
		if (NDIM >= 2) {
			m_index += (jpos+1)*3;
		}
		if (NDIM == 3) {
			m_index += (kpos+1)*9;
		}
	}

	//
	bool isCenter() const {
		return m_index == CenterIndex;
	}
	bool isFaceSurr() const {
		int ijk = abs(ipos()) + abs(jpos()) + abs(kpos());
		return ijk == 1;
	}
	bool isFaceSurr(int dir) const {
		Vector3i vpos = pos();
		if (dir == 0) {
			return (vpos.x==-1 || vpos.x==1) && vpos.y==0 && vpos.z==0;
		} else if (dir == 1) {
			return vpos.x==0 && (vpos.y==-1 || vpos.y==1) && vpos.z==0;
		} else {
			return vpos.x==0 && vpos.y==0 && (vpos.z==-1 || vpos.z==1);
		}
	}
	bool isFaceSurr(int dir, int side) const {
		Vector3i vpos = pos();
		if (dir == 0) {
			return (side==0 ? vpos.x==-1 : vpos.x==1) && vpos.y==0 && vpos.z==0;
		} else if (dir == 1) {
			return vpos.x==0 && (side==0 ? vpos.y==-1 : vpos.y==1) && vpos.z==0;
		} else {
			return vpos.x==0 && vpos.y==0 && (side==0 ? vpos.z==-1 : vpos.z==1);
		}
	}

	// conversion
	operator int() const { return m_index; }
	// allows ++
	operator int&() { return m_index; }

	// surrounding connected to given face
	static inline SurroundingIndex FaceSurrounding(const FaceIndex &face)
	{
		int dir = face.dir();
		int side = face.side();

		int pos[MAX_DIM] = { 0 };
		pos[dir] = 2*side - 1;
		assert(pos[dir]==-1 || pos[dir]==1);

		SurroundingIndex face_surr(pos[0], pos[1], pos[2]);
		return face_surr;
	}
}; // struct_surroundingindex







} // namespace_sayaka



