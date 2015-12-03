#pragma once

#include "SayakaCommons.h"
#include "SayakaBox.h"


namespace sayaka
{


/**
 * Single-grid data used as building-blocks for AMR
 */
template<typename T>
struct GridData
{
	IndexBox m_box;
	int m_ncomp;
	T *m_data;

public:
	//typedef T data_t;

	GridData()
		: m_ncomp(0), m_data(NULL)
	{}

	GridData(const IndexBox &inbox, int ncomp=1)
		: m_box(inbox),
		m_ncomp(ncomp),
		m_data(NULL)
	{
		assert(ncomp>=1);

		// allocate data
		const int len = m_box.capacity(ncomp);
		m_data = new T[len];
	}

	// copy
	GridData(const GridData & src) 
		: m_box(src.box()),
		m_ncomp(src.ncomp()),
		m_data(NULL)
	{
		// allocate data
		const int len = m_box.capacity(m_ncomp);
		m_data = new T[len];

		setValue(src);
	}

	~GridData()
	{
		if (m_data) {
			delete[] m_data;
			m_data = NULL;
		}
	}

	//
	void reset(const IndexBox &box, int ncomp) {
		m_box = box;
		m_ncomp = ncomp;
		if (m_data) {
			delete[] m_data;
			m_data = NULL;
		}
		const int size = box.capacity(ncomp);
		m_data = new T[size];
	}

	//
	int stride() const { return m_box.stride(); }
	int capacity() const { return m_box.capacity(m_ncomp); }
	int ncomp() const { return m_ncomp; }
	const IndexBox& box() const { return m_box; }
	
	//
	const T* data(int comp=0) const { 
		return m_data + m_box.stride()*comp; 
	}
	T* data(int comp=0) { 
		return m_data + m_box.stride()*comp; 
	}


	// index-based 
	const T& operator()(int i, int j, int k, int comp) const {
		int offset = m_box.offset(i, j, k, comp);
		return m_data[offset];
	}
	T& operator()(int i, int j, int k, int comp) {
		int offset = m_box.offset(i, j, k, comp);
		return m_data[offset];
	}

	// offset-based
	const T& operator[](int offset) const { 
		assert(0<=offset && offset<this->capacity());
		return m_data[offset];
	}
	T& operator[](int offset) { 
		assert(0<=offset && offset<this->capacity());
		return m_data[offset];
	}

	GridData& operator=(const GridData &rhs) {
		assert(m_box == rhs.m_box);
		assert(m_ncomp == rhs.m_ncomp);
		this->setValue(rhs);
		return *this;
	}

	//
	void setValue(const T &val) {
		const int len = m_box.capacity(m_ncomp);
		for (int i=0; i<len; i++) {
			m_data[i] = val;
		}
	}
	void setValue(int scomp, int ncomp, const T &val) {
		assert(ncomp > 0);
		assert(0<=scomp && scomp+ncomp<=m_ncomp);
		const int len = m_box.stride() * ncomp;
		T *buf = this->data(scomp);
		for (int i=0; i<len; i++) {
			buf[i] = val;
		}
	}

	void setValue(int scomp, int ncomp, const IndexBox &box, const T &val) {
		const int ilo = box.ilo();
		const int jlo = box.jlo();
		const int klo = box.klo();
		const int ihi = box.ihi();
		const int jhi = box.jhi();
		const int khi = box.khi();

		for (int comp=scomp; comp<scomp+ncomp; comp++) {
			for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
					for (int i=ilo; i<=ihi; i++) {
						(*this)(i,j,k,comp) = val;
					}
				}
			}
		}
	}

	void setValue(const GridData &src) {
		assert(this->box() == src.box());
		assert(this->ncomp() == src.ncomp());

		const int len = this->capacity(); 
		for (int i=0; i<len; i++) {
			m_data[i] = src.m_data[i];
		}
	}

}; // struct_griddata<T>


/**
 * Separated operations from the GridData.
 * This is intended for those having no arithmatic.
 */
template<typename T>
class GridDataUtil
{
public:

	// copy
	void CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &dstbox,
		const IndexBox &srcbox) const;
	//
	void CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp, 
		const IndexBox &box) const;
	//
	void CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp) const;

	// 
	void AddEqual(GridData<T> &dst, 
		int comp,
		const IndexBox &box, 
		const T& val) const;
	void AddEqual(
		GridData<T> &dst, const GridData<T> &src,
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &box) const;
	void AddEqual(
		GridData<T> &dst, const GridData<T> &src, const T &coef,
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &box) const;
	//
	void MulEqual(GridData<T> &dst,
		const IndexBox &box, int comp,
		const T& val) const;

	//
	T DotProd(const GridData<T> &a, const GridData<T> &b,
		int acomp, int bcomp, const IndexBox &box);

	// linear algebra
	void Calc_axpby(
		GridData<T> &z,
		const T &a, const GridData<T> &x,
		const T &b, const GridData<T> &y,
		const IndexBox &zbox, 
		const IndexBox &xbox, 
		const IndexBox &ybox,
		int zcomp, int xcomp, int ycomp, 
		int ncomp) const;
	
	// 
	enum ReductionOp {
		REDUCE_SUM,
		REDUCE_SUM_ABS,
		REDUCE_SUM_POW2,
		REDUCE_MAX,
		REDUCE_MIN,
		REDUCE_MAX_ABS,
		REDUCE_MIN_ABS,
		//NUM_REDUCE_OP,
	};
	T ReduceValue(
		const GridData<T> &data,
		const IndexBox &box, int comp,
		ReductionOp op) const;
	//
	T ReduceSum(const GridData<T> &data,const IndexBox &box, int comp) const {
		return ReduceValue(data, box, comp, REDUCE_SUM);
	}
	T ReduceMaxAbs(const GridData<T> &data,const IndexBox &box, int comp) const {
		return ReduceValue(data, box, comp, REDUCE_MAX_ABS);
	}
};



//
typedef GridData<double> DoubleBlockData;
typedef GridData<int> IntBlockData;

//
extern GridDataUtil<double> DoubleGridDataUtil;
extern GridDataUtil<int> IntGridDataUtil;





} // namespace_sayaka


