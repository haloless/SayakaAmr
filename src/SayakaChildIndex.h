#pragma once


#include "SayakaCommons.h"
#include "SayakaFaceIndex.h"
#include "SayakaBox.h"

#include <utility>
#include <tuple>
#include <iosfwd>

SAYAKA_NS_BEGIN;

//
// Index of 4(2D) or 8(3D) child nodes.
//
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

public:

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

	TreeChildIndex(const TreeChildIndex &rhs) = default;

	// conversion to integer, allows direct indexing
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
	
	//
	// isSameBlock will be TRUE if the neighbor is in-block sibling.
	// Otherwise the neighbor is a child from adjacent block.
	//
	TreeChildIndex neighbor(int iface, bool &isSameBlock) const;

	TreeChildIndex neighbor(int idir, int iside, bool &isSameBlock) const {
		assert(0<=idir && idir<MAX_DIM);
		assert(0<=iside && iside<2);

		FaceIndex iface(idir, iside);
		return neighbor(iface, isSameBlock);
	}

	//auto neighbor(int iface) const {
	//	bool isSameBlock;
	//	TreeChildIndex ichild = neighbor(iface, isSameBlock);
	//	return std::make_tuple(ichild, isSameBlock);
	//}

	// BoundBox of child node
	static RealBox ChildBoundBox(
		const RealBox &parentBox, 
		const TreeChildIndex &ichild);

	// Offset of index of child in terms of parent box
	static Vector3i ChildIndexOffset(
		const IndexBox &parentBox,
		const TreeChildIndex &ichild);

	// Box of child 
	static IndexBox ChildSubIndexRange(
		const IndexBox &parentBox, const TreeChildIndex &ichild);


}; // struct_childindex

// use this alias instead
using ChildIndex = TreeChildIndex;


// IO
std::ostream& operator<<(std::ostream &os, const ChildIndex &ichild);



SAYAKA_NS_END;


