

#include "SayakaBlockData.h"


namespace sayaka
{

	/*
	 * Member functions of GridDataUtil<T>
	 */

	// copy in-box value, slow?
	template<typename T>
	void GridDataUtil<T>::CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &dstbox,
		const IndexBox &srcbox) const
	{
		assert(dst.box().type() == src.box().type());
		// data capacity
		assert(dstcomp+ncomp <= dst.ncomp());
		assert(srccomp+ncomp <= src.ncomp());
		// box equivalence
		assert(dstbox.size() == srcbox.size());

		const int &ilo = dstbox.ilo();
		const int &ihi = dstbox.ihi();
		const int &jlo = dstbox.jlo();
		const int &jhi = dstbox.jhi();
		const int &klo = dstbox.klo();
		const int &khi = dstbox.khi();

		const int ioff = srcbox.ilo() - ilo;
		const int joff = srcbox.jlo() - jlo;
		const int koff = srcbox.klo() - klo;

		for (int icomp=0; icomp<ncomp; icomp++) {
			for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
					for (int i=ilo; i<=ihi; i++) {
						dst(i,j,k,dstcomp+icomp) = 
							src(i+ioff,j+joff,k+koff,srccomp+icomp);
					}
				}
			}
		}
	}
	
	template<typename T>
	void GridDataUtil<T>::CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp, 
		const IndexBox &box) const
	{
		assert(dst.box().type() == src.box().type());
		assert(dstcomp+ncomp <= dst.ncomp());
		assert(srccomp+ncomp <= src.ncomp());

		const int &ilo = box.ilo();
		const int &ihi = box.ihi();
		const int &jlo = box.jlo();
		const int &jhi = box.jhi();
		const int &klo = box.klo();
		const int &khi = box.khi();

		for (int icomp=0; icomp<ncomp; icomp++) {
			for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
					for (int i=ilo; i<=ihi; i++) {
						dst(i,j,k,dstcomp+icomp) = src(i,j,k,srccomp+icomp);
					}
				}
			}
		}
	}

	template<typename T>
	void GridDataUtil<T>::CopyData(
		GridData<T> &dst, const GridData<T> &src, 
		int dstcomp, int srccomp, int ncomp) const
	{
		assert(dst.box().type() == src.box().type());
		assert(dst.box() == src.box());
		assert(dstcomp+ncomp <= dst.ncomp());
		assert(srccomp+ncomp <= src.ncomp());

		const int len = dst.stride();
		for (int icomp=0; icomp<ncomp; icomp++) {
			T *dstdata = dst.data(dstcomp+icomp);
			const T* srcdata = src.data(srccomp+icomp);
			for (int i=0; i<len; i++) {
				dstdata[i] = srcdata[i];
			}
		}
	}

	template<typename T>
	void GridDataUtil<T>::AddEqual(
		GridData<T> &dst, int comp,
		const IndexBox &box, 
		const T& val) const
	{
		assert(dst.box().contains(box));
		assert(0<=comp && comp<dst.ncomp());

		decl_box_range(box, i,j,k);

		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					dst(i,j,k,comp) += val;
				}
			}
		}
	}
	template<typename T>
	void GridDataUtil<T>::AddEqual(
		GridData<T> &dst, const GridData<T> &src,
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &box) const
	{
		assert(dst.box().contains(box));
		assert(src.box().contains(box));
		assert(0<=dstcomp && dstcomp+ncomp<=dst.ncomp());
		assert(0<=srccomp && srccomp+ncomp<=src.ncomp());

		decl_box_range(box, i,j,k);

		for (int icomp=0; icomp<ncomp; icomp++) {
			for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
					for (int i=ilo; i<=ihi; i++) {
						dst(i,j,k,dstcomp+icomp) += src(i,j,k,srccomp+icomp);
					}
				}
			}
		}
	}
	template<typename T>
	void GridDataUtil<T>::AddEqual(
		GridData<T> &dst, const GridData<T> &src, const T &coef,
		int dstcomp, int srccomp, int ncomp,
		const IndexBox &box) const
	{
		assert(dst.box().contains(box));
		assert(src.box().contains(box));
		assert(0<=dstcomp && dstcomp+ncomp<=dst.ncomp());
		assert(0<=srccomp && srccomp+ncomp<=src.ncomp());

		decl_box_range(box, i,j,k);

		for (int icomp=0; icomp<ncomp; icomp++) {
			for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
					for (int i=ilo; i<=ihi; i++) {
						dst(i,j,k,dstcomp+icomp) += coef * src(i,j,k,srccomp+icomp);
					}
				}
			}
		}
	}
	template<typename T>
	void GridDataUtil<T>::MulEqual(
		GridData<T> &dst,
		const IndexBox &box, int comp,
		const T& val) const
	{
		assert(dst.box().contains(box));
		assert(0<=comp && comp<dst.ncomp());

		decl_box_range(box, i,j,k);

		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					dst(i,j,k,comp) *= val;
				}
			}
		}
	}

	// dot product
	template<typename T>
	T GridDataUtil<T>::DotProd(
		const GridData<T> &a, const GridData<T> &b,
		int acomp, int bcomp, const IndexBox &box)
	{
		assert(a.box().contains(box));
		assert(b.box().contains(box));
		assert(0<=acomp && acomp<a.ncomp());
		assert(0<=bcomp && bcomp<b.ncomp());

		decl_box_range(box, i,j,k);

		T dot_prod = (T) 0;

		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					dot_prod += a(i,j,k,acomp) * b(i,j,k,bcomp);
				}
			}
		}

		return dot_prod;
	}

	// linear algebra
	// z = a*x + b*y
	template<typename T>
	void GridDataUtil<T>::Calc_axpby(
		GridData<T> &z,
		const T& a, const GridData<T> &x,
		const T& b, const GridData<T> &y,
		const IndexBox &zbox, 
		const IndexBox &xbox, 
		const IndexBox &ybox,
		int zcomp, int xcomp, int ycomp, int ncomp) const
	{
		// check input box type
		assert(zbox.type() == z.box().type());
		assert(xbox.type() == x.box().type());
		assert(ybox.type() == y.box().type());
		// check data box type
		assert(z.box().type() == x.box().type());
		assert(z.box().type() == y.box().type());
		// check box size
		assert(zbox.size() == xbox.size());
		assert(zbox.size() == ybox.size());
		// check data capacity
		assert(0<=zcomp && zcomp+ncomp<=z.ncomp());
		assert(0<=xcomp && zcomp+ncomp<=x.ncomp());
		assert(0<=ycomp && zcomp+ncomp<=y.ncomp());


		// target box range
		decl_box_range(zbox, i,j,k);
		//
		const int ixoff = xbox.ilo() - ilo;
		const int jxoff = xbox.jlo() - jlo;
		const int kxoff = xbox.klo() - klo;
		//
		const int iyoff = ybox.ilo() - ilo;
		const int jyoff = ybox.jlo() - jlo;
		const int kyoff = ybox.klo() - klo;

		for (int icomp=0; icomp<ncomp; icomp++) {
			for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				double saxby = 
					a * x(i+ixoff,j+jxoff,k+kxoff,xcomp+icomp) +
					b * y(i+iyoff,j+jyoff,k+kyoff,ycomp+icomp);
				z(i,j,k,zcomp+icomp) = saxby;
			}
			}
			}
		}
	}


	template<typename T>
	T GridDataUtil<T>::ReduceValue(
		const GridData<T> &data,
		const IndexBox &box, int comp,
		ReductionOp op) const
	{
		assert(box.type() == data.box().type());
		assert(data.box().contains(box));
		assert(0<=comp && comp<data.ncomp());
		
		//
		decl_box_range(box, i,j,k);

		T result = (T) 0;
		if (op==REDUCE_MAX || op==REDUCE_MIN) {
			result = data(ilo,jlo,klo,comp);
		} else if (op==REDUCE_MAX_ABS || op==REDUCE_MIN_ABS) {
			result = abs(data(ilo,jlo,klo,comp));
		}

		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					const T &val = data(i,j,k,comp);
					if (op == REDUCE_SUM) {
						result += val;
					} else if (op == REDUCE_SUM_ABS) {
						result += abs(val);
					} else if (op == REDUCE_SUM_POW2) {
						result += val*val;
					} else if (op == REDUCE_MAX) {
						result = std::max(result, val);
					} else if (op == REDUCE_MIN) {
						result = std::min(result, val);
					} else if (op == REDUCE_MAX_ABS) {
						result = std::max(result, abs(val));
					} else if (op == REDUCE_MIN_ABS) {
						result = std::min(result, abs(val));
					}
				}
			}
		}

		return result;
	}



	/*
	 * enforce template explicit instantiation
	 */
	template class GridDataUtil<double>;
	template class GridDataUtil<int>;
	GridDataUtil<double> DoubleGridDataUtil;
	GridDataUtil<int> IntGridDataUtil;


} // namespace_sayaka




