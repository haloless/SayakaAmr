

#include "SayakaChildIndex.h"



namespace sayaka
{

//
ChildIndex ChildIndex::neighbor(int iface, bool &isSameBlock) const {
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
		{ 1, 0, 0, 1, 0, 1 },
		{ 0, 1, 1, 0, 0, 1 },
		{ 1, 0, 1, 0, 0, 1 },
		{ 0, 1, 0, 1, 1, 0 },
		{ 1, 0, 0, 1, 1, 0 },
		{ 0, 1, 1, 0, 1, 0 },
		{ 1, 0, 1, 0, 1, 0 },
	};
	
	assert(0<=iface && iface<MAX_FACE);

	isSameBlock = blk[m_index][iface];
	return nbr[m_index][iface];
}






} // namespace_sayaka



