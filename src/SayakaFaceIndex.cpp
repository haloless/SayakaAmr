

#include "SayakaFaceIndex.h"



SAYAKA_NS_BEGIN;

//
BlockFaceIndex BlockFaceIndex::opposite() const {
	const static int opp[] = {
		1, 0, 
		3, 2, 
		5, 4,
	};
	return opp[m_index]; 
}

//--------------------------------------------------------------------------------

std::ostream & operator<<(std::ostream & os, const FaceIndex & iface)
{
	os << static_cast<int>(iface);
	return os;
}


SAYAKA_NS_END;

