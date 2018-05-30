

#include "SayakaChildIndex.h"

#include <iostream>


SAYAKA_NS_BEGIN;


//
TreeChildIndex TreeChildIndex::neighbor(int iface, bool &isSameBlock) const {
	// index map of neighbor 
	const static int nbr[MAX_CHILD][MAX_FACE] = {
		//xl,xh,yl,yh,zl,zh
		{ 1, 1, 2, 2, 4, 4 }, // 0
		{ 0, 0, 3, 3, 5, 5 }, // 1
		{ 3, 3, 0, 0, 6, 6 }, // 2
		{ 2, 2, 1, 1, 7, 7 }, // 3
		{ 5, 5, 6, 6, 0, 0 }, // 4
		{ 4, 4, 7, 7, 1, 1 }, // 5
		{ 7, 7, 4, 4, 2, 2 }, // 6
		{ 6, 6, 5, 5, 3, 3 }, // 7
	};
	// boolean map of whether have the same parent
	const static int blk[MAX_CHILD][MAX_FACE] = {
		//xl,xh,yl,yh,zl,zh
		{ 0, 1, 0, 1, 0, 1 }, // 0
		{ 1, 0, 0, 1, 0, 1 }, // 1
		{ 0, 1, 1, 0, 0, 1 }, // 2
		{ 1, 0, 1, 0, 0, 1 }, // 3
		{ 0, 1, 0, 1, 1, 0 }, // 4
		{ 1, 0, 0, 1, 1, 0 }, // 5
		{ 0, 1, 1, 0, 1, 0 }, // 6
		{ 1, 0, 1, 0, 1, 0 }, // 7
	};
	
	assert(0<=iface && iface<MAX_FACE);

	isSameBlock = blk[m_index][iface];
	return nbr[m_index][iface];
}

//

RealBox TreeChildIndex::ChildBoundBox(
	const RealBox & parentBox, 
	const TreeChildIndex & ichild)
{
	RealBox childBox = parentBox;
	for (int dir = 0; dir<NDIM; dir++) {
		if (ichild.side(dir) == 0) {
			childBox.vhi(dir) = parentBox.center(dir);
		}
		else {
			childBox.vlo(dir) = parentBox.center(dir);
		}
	}
	return childBox;
}

// Offset of index of child box in parent box
Vector3i TreeChildIndex::ChildIndexOffset(
	const IndexBox & parentBox, 
	const TreeChildIndex & ichild)
{
	Vector3i voffset;
	voffset.setZero();

	for (int dir = 0; dir<NDIM; dir++) {

		int side = ichild.side(dir);
		int nb = parentBox.size(dir);
		// currently only even block division is allowed
		assert(nb % 2 == 0);

		voffset(dir) = side * (nb / 2);
	}

	return voffset;
}

// 
IndexBox TreeChildIndex::ChildSubIndexRange(
	const IndexBox & parentBox, 
	const TreeChildIndex & ichild)
{
	assert(parentBox.isCellBox());

	IndexBox subBox = parentBox;
	for (int dir = 0; dir<NDIM; dir++) {
		// currently only even block division is allowed
		assert(parentBox.size(dir) % 2 == 0);
		int mid = (parentBox.lo(dir) + parentBox.hi(dir)) / 2;
		if (ichild.side(dir) == 0) {
			subBox.vhi(dir) = mid;
		}
		else {
			subBox.vlo(dir) = mid + 1;
		}
	}

	return subBox;
}


//--------------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & os, const ChildIndex & ichild)
{
	os << static_cast<int>(ichild);
	return os;
}


SAYAKA_NS_END;

