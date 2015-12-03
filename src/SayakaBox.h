#pragma once

#include <cassert>

#include <algorithm>

#include "vector3d.h"

#include "SayakaCommons.h"

/**
 * Box-based index
 */
namespace sayaka
{

// 
enum VarLocType {
	VARLOC_CELL = 0,
	VARLOC_FACE_X = 1, // 0b0001
	VARLOC_FACE_Y = 2, // 0b0010
	VARLOC_FACE_Z = 4, // 0b0100
	//VARLOC_NODE = 7, // 0b0111
	VARLOC_NODE = (VARLOC_FACE_X*XDIM) | (VARLOC_FACE_Y*YDIM) | (VARLOC_FACE_Z*ZDIM),
};

// Index type of variable
// Specifies its position defined on a grid
struct VariableLocation
{
	
	int m_location;

	VariableLocation() : m_location(VARLOC_CELL) {}
	VariableLocation(int loc) : m_location(loc) {}

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
		//int all_stag = 1;
		//for (int dir=0; dir<NDIM; dir++) {
		//	all_stag = all_stag && isStaggered(dir);
		//}
		//return all_stag;
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
	static inline int StaggerMask(int dir) {
		int mask = (0x1 << dir);
		return mask;
	}
};



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

	int stride() const {
		return size(0)*size(1)*size(2);
	}

	int capacity(int ncomp=1) const { 
		return stride() * ncomp; 
	}

	// check if this is a valid box
	int isValid() const {
		return vlo.x<=vhi.x && vlo.y<=vhi.y && vlo.z<=vhi.z;
	}
	// check if covers target box
	int contains(const IndexBox &b) const {
		// they must be of the same index type
		assert(this->type() == b.type());

		return vlo.x<=b.vlo.x && vlo.y<=b.vlo.y && vlo.z<=b.vlo.z &&
			vhi.x>=b.vhi.x && vhi.y>=b.vhi.y && vhi.z>=b.vhi.z;
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
	IndexBox& extendLow(int dir, int ngrow) {
		return extendInFace(dir, 0, ngrow);
	}
	IndexBox& extendHigh(int dir, int ngrow) {
		return extendInFace(dir, 1, ngrow);
	}


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


	//
	static inline IndexBox Intersection(const IndexBox &a, const IndexBox &b) {
		assert(a.type() == b.type());
		int ilo = std::max(a.ilo(), b.ilo());
		int jlo = std::max(a.jlo(), b.jlo());
		int klo = std::max(a.klo(), b.klo());
		int ihi = std::min(a.ihi(), b.ihi());
		int jhi = std::min(a.jhi(), b.jhi());
		int khi = std::min(a.khi(), b.khi());
		
		Vector3i vlo = Vector3i::VecMake(ilo, jlo, klo);
		Vector3i vhi = Vector3i::VecMake(ihi, jhi, khi);
		IndexBox c(vlo, vhi, a.type());
		return c;
	}
	static inline IndexBox Extend(const IndexBox &a, int ngrow) {
		IndexBox b(a);
		b.extend(ngrow);
		return b;
	}
	static inline IndexBox Stagger(const IndexBox &cellbox, int dir) {
		assert(cellbox.isCellBox());
		IndexBox facebox(cellbox);
		facebox.staggerInDir(dir);
		return facebox;
	}
	
	// Extract the one-layer slice inside the box most close to the face
	static inline IndexBox AdjacentFace(const IndexBox &box, int dir, int side) {
		assert(0<=dir && dir<NDIM);
		assert(0<=side && side<=1);
		IndexBox bface = box;
		if (side == 0) { // low face, collapse high end
			bface.vhi(dir) = bface.vlo(dir);
		} else { // high face, collapse low end
			bface.vlo(dir) = bface.vhi(dir);
		}
		return bface;
	}
	// Extract the one-layer slice on the box face (i.e. face-centered)
	static inline IndexBox BoundaryFace(const IndexBox &box, int dir, int side) {
		assert(0<=dir && dir<NDIM);
		assert(0<=side && side<=1);
		
		IndexBox bface = box;
		// enforce staggered in face direction
		bface.staggerInDir(dir);

		if (side == 0) { // low face, collapse high end
			bface.vhi(dir) = bface.vlo(dir);
		} else { // high face, collapse low end
			bface.vlo(dir) = bface.vhi(dir);
		}

		return bface;
	}
};

inline bool operator==(const IndexBox &b0, const IndexBox &b1) {
	return b0.vlo==b1.vlo && b0.vhi==b1.vhi && b0.varloc==b1.varloc;
}
inline bool operator!=(const IndexBox &b0, const IndexBox &b1) {
	return !(b0==b1);
}


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



} // namespace_sayaka




