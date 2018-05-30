#include "SayakaBox.h"

#include <iostream>



SAYAKA_NS_BEGIN;


std::ostream & operator<<(std::ostream & os, const VariableLocation & varloc)
{
	//os << static_cast<int>(varloc);
	os << '<' << varloc[0] << ',' << varloc[1] << ',' << varloc[2] << '>';
	return os;
}

std::ostream & operator<<(std::ostream & os, const IndexBox & box)
{
	os << '{' << box.type() << ',' << box.lo() << ',' << box.hi() << '}';
	return os;
}

std::ostream & operator<<(std::ostream & os, const RealBox & rb)
{
	os << '{' << rb.lo() << ',' << rb.hi() << '}';
	return os;
}





//

IndexBox IndexBox::Intersection(const IndexBox & a, const IndexBox & b) {
	assert(a.type() == b.type());
	//int ilo = std::max(a.ilo(), b.ilo());
	//int jlo = std::max(a.jlo(), b.jlo());
	//int klo = std::max(a.klo(), b.klo());
	//int ihi = std::min(a.ihi(), b.ihi());
	//int jhi = std::min(a.jhi(), b.jhi());
	//int khi = std::min(a.khi(), b.khi());

	//Vector3i vlo = Vector3i::VecMake(ilo, jlo, klo);
	//Vector3i vhi = Vector3i::VecMake(ihi, jhi, khi);
	//IndexBox c(vlo, vhi, a.type());
	//return c;

	return IndexBox(
		a.vlo.maxVec(b.vlo), 
		a.vhi.minVec(b.vhi), 
		a.type());
}

IndexBox IndexBox::Extend(const IndexBox & a, int ngrow) {
	IndexBox b(a);
	b.extend(ngrow);
	return b;
}

IndexBox IndexBox::Stagger(const IndexBox & cellbox, int dir) {
	assert(cellbox.isCellBox());
	IndexBox facebox(cellbox);
	facebox.staggerInDir(dir);
	return facebox;
}

// Extract the one-layer slice inside the box most close to the face
IndexBox IndexBox::AdjacentFace(const IndexBox & box, int dir, int side) {
	assert(0 <= dir && dir<NDIM);
	assert(0 <= side && side <= 1);
	IndexBox bface = box;
	if (side == 0) { // low face, collapse high end
		bface.vhi(dir) = bface.vlo(dir);
	}
	else { // high face, collapse low end
		bface.vlo(dir) = bface.vhi(dir);
	}
	return bface;
}

// Extract the one-layer slice on the box face (i.e. face-centered)
IndexBox IndexBox::BoundaryFace(const IndexBox & box, int dir, int side) {
	assert(0 <= dir && dir<NDIM);
	assert(0 <= side && side <= 1);

	IndexBox bface = box;
	// enforce staggered in face direction
	bface.staggerInDir(dir);

	if (side == 0) { // low face, collapse high end
		bface.vhi(dir) = bface.vlo(dir);
	}
	else { // high face, collapse low end
		bface.vlo(dir) = bface.vhi(dir);
	}

	return bface;
}





SAYAKA_NS_END;

