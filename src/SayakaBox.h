#pragma once

#include <cassert>

#include <algorithm>
#include <iosfwd>

#include "vector3d.h"

#include "SayakaCommons.h"

/**
 * Box-based index
 */

SAYAKA_NS_BEGIN;

// 
enum VarLocType {
	VARLOC_CELL = 0,
	VARLOC_FACE_X = 1, // 0b0001
	VARLOC_FACE_Y = 2, // 0b0010
	VARLOC_FACE_Z = 4, // 0b0100
	//VARLOC_NODE = 7, // 0b0111
	VARLOC_NODE = (VARLOC_FACE_X*XDIM) | (VARLOC_FACE_Y*YDIM) | (VARLOC_FACE_Z*ZDIM),
};

//
// Index type of variable
// Specifies its position defined on a grid
//
struct VariableLocation
{
protected:
	int m_location;

public:
	constexpr VariableLocation() : m_location(VARLOC_CELL) {}
	constexpr VariableLocation(int loc) : m_location(loc) {}

	int test(int dir) const {
		assert(0 <= dir && dir < MAX_DIM);
		return (m_location & StaggerMask(dir)) >> dir;
	}

	int operator[](int dir) const { return test(dir); }

	//
	void setStaggered(int dir) {
		//assert(0<=dir && dir<NDIM);
		int mask = StaggerMask(dir);
		m_location |= mask;
	}
	int isStaggered(int dir) const {
		//assert(0<=dir && dir<NDIM);
		int mask = StaggerMask(dir);
		return (m_location & mask);
	}

	int isCellVar() const { 
		return m_location == VARLOC_CELL; 
	}
	
	int isFaceVar(int dir) const {
		assert(0<=dir && dir<NDIM);
		return m_location == StaggerMask(dir);
	}
	int isFaceVar() const {
		for (int dir=0; dir<NDIM; dir++) {
			if (isFaceVar(dir)) return 1;
		}
		return 0;
	}
	int isFacexVar() const { return isFaceVar(0); }
	int isFaceyVar() const { return isFaceVar(1); }
	int isFacezVar() const { return isFaceVar(2); }

	// Node location is slightly different
	// (1,1,0) in 2D, (1,1,1) in 3D
	int isNodeVar() const {
		return m_location == VARLOC_NODE;
	}

	//
	int getFaceVarDir() const {
		for (int dir=0; dir<NDIM; dir++) {
			if (isFaceVar(dir)) return dir;
		}
		return -1;
	}

	// conversion
	operator int() const { return m_location; }
	// allows ++
	//operator int&() { return m_location; }

	//
	static constexpr inline int StaggerMask(int dir) {
		return (0x1 << dir);
	}
};

// IO
std::ostream& operator<<(std::ostream &os, const VariableLocation &varloc);



//
// Box of indices
//
struct IndexBox
{
	Vector3i vlo;
	Vector3i vhi;
	VariableLocation varloc;

	//
	IndexBox() {
		vlo.setZero();
		vhi.setZero();
		varloc = VARLOC_CELL;
	}
	IndexBox(const Vector3i &inlo, const Vector3i &inhi) {
		vlo = inlo;
		vhi = inhi;
		varloc = VARLOC_CELL;
	}
	IndexBox(const Vector3i &inlo, const Vector3i &inhi, const VariableLocation &inloc) {
		vlo = inlo;
		vhi = inhi;
		varloc = inloc;
	}

	//
	const VariableLocation& type() const { return varloc; }

	const Vector3i& lo() const { return vlo; }
	const Vector3i& hi() const { return vhi; }

	const int& lo(int dir) const { return vlo(dir); }
	const int& hi(int dir) const { return vhi(dir); }

	const int& ilo() const { return vlo.x; }
	const int& jlo() const { return vlo.y; }
	const int& klo() const { return vlo.z; }
	const int& ihi() const { return vhi.x; }
	const int& jhi() const { return vhi.y; }
	const int& khi() const { return vhi.z; }

	const int& endpoint(int dir, int side) const {
		assert(0<=dir && dir<NDIM);
		assert(0<=side && side<=1);

		if (side == 0) {
			return lo(dir);
		} else {
			return hi(dir);
		}
	}

	//
	int size(int dir) const { 
		return vhi(dir)-vlo(dir)+1; 
	}
	Vector3i size() const {
		Vector3i sz = vhi - vlo;
		sz.x += 1;
		sz.y += 1;
		sz.z += 1;
		return sz;
	}

	int stride() const { return size(0)*size(1)*size(2); }

	int capacity(int ncomp=1) const { return stride() * ncomp; }

	// check if this is a valid box
	int isValid() const {
		return vlo.allLE(vhi);
		//return vlo.x<=vhi.x && vlo.y<=vhi.y && vlo.z<=vhi.z;
	}
	// check if covers target box
	int contains(const IndexBox &b) const {
		// they must be of the same index type
		assert(this->type() == b.type());

		return vlo.allLE(b.vlo) && vhi.allGE(b.vhi);
		//return vlo.x<=b.vlo.x && vlo.y<=b.vlo.y && vlo.z<=b.vlo.z &&
		//	vhi.x>=b.vhi.x && vhi.y>=b.vhi.y && vhi.z>=b.vhi.z;
	}
	//// check if larger in specified direction
	//int isLarger(const IndexBox &b, int dir) const {
	//	
	//}
	//// check if at least one side larger than target box
	//int isLarger(const IndexBox &b) const {
	//	assert(this->type() == b.type());

	//}
	// check if has intersection
	int intersects(const IndexBox &b) const {
		assert(this->type() == b.type());
		//
		IndexBox c = Intersection(*this, b);
		return c.isValid();
	}

	int containsIndex(int i, int j, int k) const {
		return vlo.x<=i && vlo.y<=j && vlo.z<=k &&
			i<=vhi.x && j<=vhi.y && k<=vhi.z;
	}

	//
	int offset(int i, int j, int k, int comp=0) const {
		return (i-ilo()) + (j-jlo())*size(0) + (k-klo())*size(0)*size(1) + 
			comp*size(0)*size(1)*size(2);
	}

	void index(const int offset, int &i, int &j, int &k) const {
		int comp = -1;
		this->index(offset, i, j, k, comp);
		assert(comp == 0);
	}
	void index(const int offset, int &i, int &j, int &k, int &comp) const {
		assert(offset >= 0);
		const int sz0 = size(0);
		const int sz1 = size(1);
		const int sz2 = size(2);
		i = (offset)%sz0 + ilo();
		j = (offset/sz0)%sz1 + jlo();
		k = (offset/(sz0*sz1))%sz2 + klo();
		comp = (offset/(sz0*sz1*sz2));
	}


	//
	IndexBox& extend(const Vector3i &elo, const Vector3i &ehi) {
		vlo -= elo;
		vhi += ehi;
		return *this;
	}
	IndexBox& extend(int ielo, int jelo, int kelo, int iehi, int jehi, int kehi) {
		Vector3i elo = Vector3i::VecMake(ielo, jelo, kelo);
		Vector3i ehi = Vector3i::VecMake(iehi, jehi, kehi);
		return extend(elo, ehi);
	}
	IndexBox& extend(int iext, int jext, int kext) {
		return extend(iext, jext, kext, iext, jext, kext);
	}
	IndexBox& extend(int ext) {
		int iext = ext;
		int jext = NDIM>=2 ? ext : 0;
		int kext = NDIM==3 ? ext : 0;
		return extend(iext, jext, kext);
	}

	IndexBox& extendInFace(int dir, int side, int ngrow) {
		assert(0<=dir && dir<NDIM);
		assert(0<=side && side<=1);

		if (side == 0) { // low end
			vlo(dir) -= ngrow;
		} else { // high end
			vhi(dir) += ngrow;
		}
		return *this;
	}
	IndexBox& extendInDir(int dir, int ngrow) {
		extendInFace(dir, 0, ngrow);
		extendInFace(dir, 1, ngrow);
		return *this;
	}
	
	IndexBox& extendLo(int dir, int ngrow) { return extendInFace(dir, 0, ngrow); }
	IndexBox& extendHi(int dir, int ngrow) { return extendInFace(dir, 1, ngrow); }

	IndexBox& shrinkLo(int dir, int ng) { return extendInFace(dir, 0, -ng); }
	IndexBox& shrinkHi(int dir, int ng) { return extendInFace(dir, 1, -ng); }

	//
	IndexBox& stagger(int istag, int jstag, int kstag) {
		assert(istag==0 || istag==1);
		assert(jstag==0 || jstag==1);
		assert(kstag==0 || kstag==1);

		// check if the current box is already staggered in that direction
		if (varloc.isStaggered(0)) {
			istag = 0;
		} else {
			if (istag) varloc.setStaggered(0);
		}
		if (varloc.isStaggered(1)) {
			jstag = 0;
		} else {
			if (jstag) varloc.setStaggered(1);
		}
		if (varloc.isStaggered(2)) {
			kstag = 0;
		} else {
			if (kstag) varloc.setStaggered(2);
		}

		extend(0, 0, 0, istag, jstag, kstag);
		return *this;
	}
	IndexBox& stagger(const int stag[MAX_DIM]) {
		return stagger(stag[0], stag[1], stag[2]);
	}
	IndexBox& staggerInDir(int dir) {
		assert(0<=dir && dir<NDIM);
		
		int stag[MAX_DIM] = { 0 };
		stag[dir] = 1;
		return stagger(stag);
	}
	IndexBox& staggerAll() {
		for (int dir=0; dir<NDIM; dir++) {
			staggerInDir(dir);
		}
		return *this;
	}

	//
	IndexBox& shift(const Vector3i &vshift) {
		vlo += vshift;
		vhi += vshift;
		return *this;
	}

	int isCellBox() const { 
		return varloc.isCellVar();
	}
	int isFaceBox(int dir) const {
		return varloc.isFaceVar(dir);
	}
	int isFaceBox() const {
		for (int dir=0; dir<NDIM; dir++) {
			if (isFaceBox(dir)) return 1;
		}
		return 0;
	}
	// check box staggered in all space directions
	// NOTE that in 2D, the staggered pattern of node box is (1,1,0)
	// but in 3D the node box is (1,1,1)
	int isNodeBox() const { return varloc.isNodeVar(); }
	int isStaggeredBox(int dir) const { return varloc.isStaggered(dir); }

	// intersection
	IndexBox operator&(const IndexBox &rhs) const {
		return Intersection(*this, rhs);
	}


public:

	// return intersection of two boxes
	static IndexBox Intersection(const IndexBox &a, const IndexBox &b);

	// return a new box by extending
	static IndexBox Extend(const IndexBox &a, int ngrow);

	//
	static IndexBox Stagger(const IndexBox &cellbox, int dir);
	
	// Extract the one-layer slice inside the box most close to the face
	static IndexBox AdjacentFace(const IndexBox &box, int dir, int side);

	// Extract the one-layer slice on the box face (i.e. face-centered)
	static IndexBox BoundaryFace(const IndexBox &box, int dir, int side);
};

inline bool operator==(const IndexBox &b0, const IndexBox &b1) {
	return b0.vlo==b1.vlo && b0.vhi==b1.vhi && b0.varloc==b1.varloc;
}
inline bool operator!=(const IndexBox &b0, const IndexBox &b1) {
	return !(b0==b1);
}

// IO
std::ostream& operator<<(std::ostream &os, const IndexBox &box);


// Macro used to completely tranverse a box
// Not very useful...
#define for_box_range_k(box,ind) for (int ind=box.klo(); ind<=box.khi(); ind++)
#define for_box_range_j(box,ind) for (int ind=box.jlo(); ind<=box.jhi(); ind++)
#define for_box_range_i(box,ind) for (int ind=box.ilo(); ind<=box.ihi(); ind++)

#define decl_box_range_k(box,ind) const int& ind##lo = box.klo(); \
	const int& ind##hi = box.khi();
#define decl_box_range_j(box,ind) const int& ind##lo = box.jlo(); \
	const int& ind##hi = box.jhi();
#define decl_box_range_i(box,ind) const int& ind##lo = box.ilo(); \
	const int& ind##hi = box.ihi();
#define decl_box_range(box,i,j,k) decl_box_range_i(box,i); \
	decl_box_range_j(box,j); \
	decl_box_range_k(box,k);

#define BEGIN_FOR_BOX_RANGE(box, ii,jj,kk) {\
const int ii##lo = box.ilo(); const int ii##hi = box.ihi(); \
const int jj##lo = box.jlo(); const int jj##hi = box.jhi(); \
const int kk##lo = box.klo(); const int kk##hi = box.khi(); \
for (int kk=kk##lo; kk<=kk##hi; ++kk) {\
for (int jj=jj##lo; jj<=jj##hi; ++jj) {\
for (int ii=ii##lo; ii<=ii##hi; ++ii) {\

#define END_FOR_BOX_RANGE }}}}



//
// Physical box region
//
struct RealBox
{
	Vector3d vlo;
	Vector3d vhi;

	//
	RealBox() 
	{
		vlo.setZero();
		vhi.setZero();
	}
	RealBox(const Vector3d& inlo, const Vector3d& inhi) 
	{
		vlo = inlo;
		vhi = inhi;
	}
	RealBox(const double xlo[], const double xhi[]) 
	{
		vlo.setZero();
		vhi.setZero();
		for (int dir = 0; dir < NDIM; dir++) {
			vlo(dir) = xlo[dir];
			vhi(dir) = xhi[dir];
		}
	}


	const Vector3d& lo() const { return vlo; }
	const Vector3d& hi() const { return vhi; }

	double center(int dir) const {
		return 0.5 * (vlo(dir) + vhi(dir));
	}
	double length(int dir) const {
		return vhi(dir)-vlo(dir);
	}
	Vector3d center() const {
		Vector3d vcen = 0.5 * (vlo+vhi);
		return vcen;
	}
	Vector3d length() const {
		Vector3d vlen = vhi - vlo;
		return vlen;
	}
}; // struct_realbox

// IO
std::ostream& operator<<(std::ostream &os, const RealBox &rb);



SAYAKA_NS_END;



