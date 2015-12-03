#pragma once


#include "SayakaCommons.h"
#include "SayakaFaceIndex.h"
#include "SayakaBox.h"

namespace sayaka
{


struct TreeChildIndex
{
	enum {
		NxChild = 2,
		NyChild = 2,
		NzChild = NDIM==2 ? 1 : 2,

		NumChild = NxChild*NyChild*NzChild,
		NumChildOnFace = NumChild/2,
	};

	int m_index;

	TreeChildIndex(int index)
		: m_index(index)
	{
		assert(0<=index && index<MAX_CHILD);
	}
	TreeChildIndex(int iside, int jside, int kside)
		: m_index(iside + jside*2 + kside*4)
	{
		assert(0<=iside && iside<2);
		assert(0<=jside && jside<2);
		assert(0<=kside && kside<2);
	}

	// conversion
	operator int() const { return m_index; }
	// allows ++
	operator int&() { return m_index; }

	int iside() const { return m_index%2; }
	int jside() const { return (m_index/2)%2; }
	int kside() const { return m_index/4; }
	int side(int dir) const {
		int ind = m_index >> dir;
		return ind & 0x1;
	}
	void sides(int sides[MAX_DIM]) const {
		for (int dir=0; dir<MAX_DIM; dir++) {
			sides[dir] = side(dir);
		}
	}


	bool hasFaceBoundary(const FaceIndex &iface) const {
		int facedir = iface.dir();
		return side(facedir) == iface.side();
	}

	//
	TreeChildIndex neighborInDir(int dir) const {
		int ijk[MAX_DIM] = { 0 };
		this->sides(ijk);
		ijk[dir] = 1 - ijk[dir]; // swap another side
		return TreeChildIndex(ijk[0], ijk[1], ijk[2]);
	}
	
	// isSameBlock will be TRUE if the neighbor is in-block sibling.
	// Otherwise the neighbor is a child from adjacent block.
	TreeChildIndex neighbor(int iface, bool &isSameBlock) const;
	TreeChildIndex neighbor(int idir, int iside, bool &isSameBlock) const {
		assert(0<=idir && idir<MAX_DIM);
		assert(0<=iside && iside<2);

		FaceIndex iface(idir, iside);
		return neighbor(iface, isSameBlock);
	}

	//
	static inline RealBox ChildBoundBox(
		const RealBox &parentBox, 
		const TreeChildIndex &ichild) 
	{
		RealBox childBox = parentBox;
		for (int dir=0; dir<NDIM; dir++) {
			if (ichild.side(dir) == 0) {
				childBox.vhi(dir) = parentBox.center(dir);
			} else {
				childBox.vlo(dir) = parentBox.center(dir);
			}
		}
		return childBox;
	}
	static inline Vector3i ChildIndexOffset(
		const IndexBox &parentBox,
		const TreeChildIndex &ichild)
	{		
		Vector3i voffset; 
		voffset.setZero();

		for (int dir=0; dir<NDIM; dir++) {

			int side = ichild.side(dir);
			int nb = parentBox.size(dir);
			// currently only even block division is allowed
			assert(nb%2 == 0);

			voffset(dir) = side * (nb/2);
		}

		return voffset;
	}
	static inline IndexBox ChildSubIndexRange(
		const IndexBox &parentBox, const TreeChildIndex &ichild)
	{
		assert(parentBox.isCellBox());

		IndexBox subBox = parentBox;
		for (int dir=0; dir<NDIM; dir++) {
			// currently only even block division is allowed
			assert(parentBox.size(dir)%2 == 0);
			int mid = (parentBox.lo(dir)+parentBox.hi(dir)) / 2;
			if (ichild.side(dir) == 0) {
				subBox.vhi(dir) = mid;
			} else {
				subBox.vlo(dir) = mid+1;
			}
		}

		return subBox;
	}
}; // struct_childindex

// use this name
typedef TreeChildIndex ChildIndex;


} // namespace_sayaka



