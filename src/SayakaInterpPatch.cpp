
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <memory>

#include "log.h"

#include "SayakaInterpPatch.h"


namespace sayaka
{

void InterpPatch::prolongCrseToFineBox(
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, const Vector3i &fineoffset,
	int scomp, int ncomp) const
{
	assert(0<=scomp);
	assert(scomp+ncomp<=finedata.ncomp());
	assert(scomp+ncomp<=crsedata.ncomp());

	for (int comp=scomp; comp<scomp+ncomp; comp++) {
		const int prolong_method = prolong_flag[comp];

		if (var_loc.isCellVar()) {
			// cell-centered
			if (prolong_method == PROLONG_IGNORE) {
				continue;
			} else if (prolong_method == PROLONG_INJECTION) {
				prolongCellBlockData_Injection(
					ifine, finedata,
					icrse, crsedata,
					finebox, fineoffset,
					comp);
			} else if (prolong_method == PROLONG_CENTERED_LINEAR) {
				prolongCellBlockData_CenteredLinear(
					ifine, finedata,
					icrse, crsedata,
					finebox, fineoffset,
					comp);
			} else if (prolong_method == PROLONG_CENTERED_LINEAR_LIMITED) {
				prolongCellBlockData_CenteredLinearLimited(
					ifine, finedata,
					icrse, crsedata,
					finebox, fineoffset,
					comp);
			}
			else {
				LOGPRINTF("%s: unknown cell prolongation=%d\n", __FUNCTION__, prolong_method);
				exit(1);
			}
		} else if (var_loc.isFaceVar()) {
			// face-centered
			const int var_dir = var_loc.getFaceVarDir();

			if (prolong_method == PROLONG_IGNORE) {
				continue;
			} else if (prolong_method == PROLONG_INJECTION) {
				prolongFaceBlockData_Injection(var_dir,
					ifine, finedata,
					icrse, crsedata,
					finebox, fineoffset,
					comp);
			} else if (prolong_method==PROLONG_CENTERED_LINEAR ||
				prolong_method==PROLONG_CENTERED_LINEAR_LIMITED) 
			{
				int use_limiter = 
					prolong_method==PROLONG_CENTERED_LINEAR_LIMITED
					? 1 : 0;
				prolongFaceBlockData_CenteredLinear(var_dir,
					ifine, finedata,
					icrse, crsedata,
					finebox, fineoffset,
					comp, use_limiter);
			}
			else {
				LOGPRINTF("%s: unknown face prolongation=%d\n", __FUNCTION__, prolong_method);
				exit(1);
			}
		} 
		else {
			LOGPRINTF("%s: unsupported data index type=%d!\n", __FUNCTION__, static_cast<int>(var_loc));
			exit(1);
		}
	}
}


void InterpPatch::restrictFineToCrseBox(
	int icrse, DoubleBlockData &crsedata,
	int ifine, const DoubleBlockData &finedata,
	const IndexBox &crsebox,
	const Vector3i &fineoffset,
	int scomp, int ncomp) const
{
	assert(0<=scomp);
	assert(scomp+ncomp<=finedata.ncomp());
	assert(scomp+ncomp<=crsedata.ncomp());

	for (int comp=scomp; comp<scomp+ncomp; comp++) {
		const int restrict_method = restrict_flag[comp];

		if (var_loc.isCellVar()) {
			// cell-centered
			if (restrict_method == RESTRICT_IGNORE) {
				continue;
			} else if (restrict_method==RESTRICT_AVERAGE || 
				restrict_method==RESTRICT_SUM) 
			{
				int do_average = restrict_method==RESTRICT_AVERAGE ? 1 : 0;
				restrictCellBlockData_Average(
					icrse, crsedata,
					ifine, finedata,
					crsebox, fineoffset,
					comp, do_average);
			} else {
				LOGPRINTF("%s: unknown cell restriction=%d\n", __FUNCTION__, restrict_method);
				exit(1);
			}
		} else if (var_loc.isFaceVar()) {
			// face-centered
			const int var_dir = var_loc.getFaceVarDir();

			if (restrict_method == RESTRICT_IGNORE) {
				continue;
			} else if (restrict_method==RESTRICT_AVERAGE || 
				restrict_method==RESTRICT_SUM) 
			{
				int do_average = restrict_method==RESTRICT_AVERAGE ? 1 : 0;
				restrictFaceBlockData_Average(var_dir,
					icrse, crsedata,
					ifine, finedata,
					crsebox, fineoffset,
					comp, do_average);
			}
			else {
				LOGPRINTF("%s: unknown face restriction=%d\n", __FUNCTION__, restrict_method);
				exit(1);
			}
		}
		else {
			LOGPRINTF("%s: unsupported data index type=%d\n", __FUNCTION__, static_cast<int>(var_loc));
			exit(1);
		}
	}
}




void InterpPatch::prolongCellBlockData_Injection(
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, 
	const Vector3i &fineoffset,
	int comp) const
{
	//
	const int ioff = fineoffset.x;
	const int joff = fineoffset.y;
	const int koff = fineoffset.z;

	const int ilo = finebox.ilo();
	const int jlo = finebox.jlo();
	const int klo = finebox.klo();
	const int ihi = finebox.ihi();
	const int jhi = finebox.jhi();
	const int khi = finebox.khi();

	for (int k=klo; k<=khi; k++) {
		int kc = IndexMapping::DirFineToCrse<2>(k, koff);

		for (int j=jlo; j<=jhi; j++) {
			int jc = IndexMapping::DirFineToCrse<1>(j, joff);

			for (int i=ilo; i<=ihi; i++) {
				int ic = IndexMapping::DirFineToCrse<0>(i, ioff);

				finedata(i,j,k,comp) = crsedata(ic,jc,kc,comp);
			}
		}
	}
}

void InterpPatch::prolongCellBlockData_Bilinear(
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, 
	const Vector3i &fineoffset,
	int comp) const
{
	LOGPRINTF("%s: not implemented\n", __FUNCTION__);

	//
	const int ioff = fineoffset.x;
	const int joff = fineoffset.y;
	const int koff = fineoffset.z;

	const int ilo = finebox.ilo();
	const int jlo = finebox.jlo();
	const int klo = finebox.klo();
	const int ihi = finebox.ihi();
	const int jhi = finebox.jhi();
	const int khi = finebox.khi();

	for (int k=klo; k<=khi; k++) {
		int kc = IndexMapping::DirFineToCrse<2>(k, koff);

		for (int j=jlo; j<=jhi; j++) {
			int jc = IndexMapping::DirFineToCrse<1>(j, joff);

			for (int i=ilo; i<=ihi; i++) {
				int ic = IndexMapping::DirFineToCrse<0>(i, ioff);

				finedata(i,j,k,comp) = crsedata(ic,jc,kc,comp);
			}
		}
	}
}


void InterpPatch::prolongCellBlockData_CenteredLinear(
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, const Vector3i &fineoffset,
	int comp) const
{
	// NOTE that centered linear interpolation requires at least NGHOST>=2
	// This is because that with NGROW=1, 
	if (tree_data.numGrow() < 2) {
		LOGPRINTF("%s: NGROW=%d < 2\n", __FUNCTION__, tree_data.numGrow());
		exit(1);
	}

	//
	const int ioff = fineoffset.x;
	const int joff = fineoffset.y;
	const int koff = fineoffset.z;

	const int ilo = finebox.ilo();
	const int jlo = finebox.jlo();
	const int klo = finebox.klo();
	const int ihi = finebox.ihi();
	const int jhi = finebox.jhi();
	const int khi = finebox.khi();

	const int iclo = IndexMapping::DirFineToCrse<0>(ilo, ioff);
	const int ichi = IndexMapping::DirFineToCrse<0>(ihi, ioff);
	const int jclo = IndexMapping::DirFineToCrse<1>(jlo, joff);
	const int jchi = IndexMapping::DirFineToCrse<1>(jhi, joff);
	const int kclo = IndexMapping::DirFineToCrse<2>(klo, koff);
	const int kchi = IndexMapping::DirFineToCrse<2>(khi, koff);

	// check index out of bound
	assert(crsedata.box().ilo()<iclo && ichi<crsedata.box().ihi());
	if (NDIM >= 2) { assert(crsedata.box().jlo()<jclo && jchi<crsedata.box().jhi()); }
	if (NDIM == 3) { assert(crsedata.box().klo()<kclo && kchi<crsedata.box().khi()); }

	for (int kc=kclo; kc<=kchi; kc++) {
		for (int jc=jclo; jc<=jchi; jc++) {
			for (int ic=iclo; ic<=ichi; ic++) {
				// calculate gradient by centered difference
				// NOTE it is measured in coarse cell unit
				Vector3d crsegrad;
				crsegrad.setZero();

				for (int dir=0; dir<NDIM; dir++) {
					int ii, jj, kk;
					sayaka::SelectDirIncr(dir, ii, jj, kk);

					crsegrad(dir) = 0.5 * (
						crsedata(ic+ii,jc+jj,kc+kk,comp) - 
						crsedata(ic-ii,jc-jj,kc-kk,comp));
				}

				// interpolation to 8 (or 4) fine cells
				for (int kk=-ZDIM; kk<=ZDIM; kk+=2) {
					int k = IndexMapping::DirCrseToFine<2>(kc, (kk+1)/2, koff);
					double krel = 0.25 * (double) kk;
					if (k<klo || k>khi) continue; // this may happens for odd number ghost cells

					for (int jj=-YDIM; jj<=YDIM; jj+=2) {
						int j = IndexMapping::DirCrseToFine<1>(jc, (jj+1)/2, joff);
						double jrel = 0.25 * (double) jj;
						if (j<jlo || j>jhi) continue;

						for (int ii=-XDIM; ii<=XDIM; ii+=2) {
							int i = IndexMapping::DirCrseToFine<0>(ic, (ii+1)/2, ioff);
							double irel = 0.25 * (double) ii;
							if (i<ilo || i>ihi) continue;

							finedata(i,j,k,comp) = crsedata(ic,jc,kc,comp) + 
								irel*crsegrad.x + jrel*crsegrad.y + krel*crsegrad.z;
						}
					}
				}

			} // coarse i
		} // coarse j
	} // coarse k

}

void InterpPatch::prolongCellBlockData_CenteredLinearLimited(
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, const Vector3i &fineoffset,
	int comp) const
{
	// NOTE that centered linear interpolation requires at least NGHOST>=2
	// This is because that with NGROW=1, 
	if (tree_data.numGrow() < 2) {
		LOGPRINTF("%s: NGROW=%d < 2\n", __FUNCTION__, tree_data.numGrow());
		exit(1);
	}

	//
	const int ioff = fineoffset.x;
	const int joff = fineoffset.y;
	const int koff = fineoffset.z;

	const int ilo = finebox.ilo();
	const int jlo = finebox.jlo();
	const int klo = finebox.klo();
	const int ihi = finebox.ihi();
	const int jhi = finebox.jhi();
	const int khi = finebox.khi();

	const int iclo = IndexMapping::DirFineToCrse<0>(ilo, ioff);
	const int ichi = IndexMapping::DirFineToCrse<0>(ihi, ioff);
	const int jclo = IndexMapping::DirFineToCrse<1>(jlo, joff);
	const int jchi = IndexMapping::DirFineToCrse<1>(jhi, joff);
	const int kclo = IndexMapping::DirFineToCrse<2>(klo, koff);
	const int kchi = IndexMapping::DirFineToCrse<2>(khi, koff);

	// check index out of bound
	assert(crsedata.box().ilo()<iclo && ichi<crsedata.box().ihi());
	if (NDIM >= 2) { assert(crsedata.box().jlo()<jclo && jchi<crsedata.box().jhi()); }
	if (NDIM == 3) { assert(crsedata.box().klo()<kclo && kchi<crsedata.box().khi()); }

	for (int kc=kclo; kc<=kchi; kc++) {
		for (int jc=jclo; jc<=jchi; jc++) {
			for (int ic=iclo; ic<=ichi; ic++) {
				// calculate gradient by centered difference with limiter
				// NOTE it is measured in coarse cell unit
				Vector3d crsegrad;
				crsegrad.setZero();

				for (int dir=0; dir<NDIM; dir++) {
					int ii, jj, kk;
					sayaka::SelectDirIncr(dir, ii, jj, kk);

					double grad = limited_slope(
						crsedata(ic-ii,jc-jj,kc-kk,comp),
						crsedata(ic,jc,kc,comp),
						crsedata(ic+ii,jc+jj,kc+kk,comp));

					crsegrad(dir) = grad;
				}

				// interpolation to fine cells
				for (int kk=-ZDIM; kk<=ZDIM; kk+=2) {
					int k = IndexMapping::DirCrseToFine<2>(kc, (kk+1)/2, koff);
					double krel = 0.25 * (double) kk;
					if (k<klo || k>khi) continue; // this may happens for odd number ghost cells

					for (int jj=-YDIM; jj<=YDIM; jj+=2) {
						int j = IndexMapping::DirCrseToFine<1>(jc, (jj+1)/2, joff);
						double jrel = 0.25 * (double) jj;
						if (j<jlo || j>jhi) continue;

						for (int ii=-XDIM; ii<=XDIM; ii+=2) {
							int i = IndexMapping::DirCrseToFine<0>(ic, (ii+1)/2, ioff);
							double irel = 0.25 * (double) ii;
							if (i<ilo || i>ihi) continue;

							finedata(i,j,k,comp) = crsedata(ic,jc,kc,comp) + 
								irel*crsegrad.x + jrel*crsegrad.y + krel*crsegrad.z;
						}
					}
				}
			}
		}
	}
}






void InterpPatch::restrictCellBlockData_Average(
	int icrse, DoubleBlockData &crsedata,
	int ifine, const DoubleBlockData &finedata,
	const IndexBox &crsebox,
	const Vector3i &fineoffset,
	int comp, int do_average) const
{
	const int &ioff = fineoffset.x;
	const int &joff = fineoffset.y;
	const int &koff = fineoffset.z;

	const int &iclo = crsebox.ilo();
	const int &ichi = crsebox.ihi();
	const int &jclo = crsebox.jlo();
	const int &jchi = crsebox.jhi();
	const int &kclo = crsebox.klo();
	const int &kchi = crsebox.khi();

	const double wgt = 1.0 / ChildIndex::NumChild;

	for (int kc=kclo; kc<=kchi; kc++) {
		for (int jc=jclo; jc<=jchi; jc++) {
			for (int ic=iclo; ic<=ichi; ic++) {
				double sum = 0;

				// loop children
				for (int kk=0; kk<=ZDIM; kk++) {
					int kfine = IndexMapping::DirCrseToFine<2>(kc, kk, koff);
					for (int jj=0; jj<=YDIM; jj++) {
						int jfine = IndexMapping::DirCrseToFine<1>(jc, jj, joff);
						for (int ii=0; ii<=XDIM; ii++) {
							int ifine = IndexMapping::DirCrseToFine<0>(ic, ii, ioff);

							sum += finedata(ifine,jfine,kfine,comp);
						}
					}
				}

				if (do_average) {
					crsedata(ic,jc,kc,comp) = sum * wgt;
				} else {
					crsedata(ic,jc,kc,comp) = sum;
				}
			}
		}
	}
}


void InterpPatch::prolongFaceBlockData_Injection(int dir,
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, const Vector3i &fineoffset,
	int comp) const
{
	assert(0<=dir && dir<NDIM);
	assert(var_loc.isFaceVar(dir));
	assert(finebox.isFaceBox(dir));

	const int &ioff = fineoffset.x;
	const int &joff = fineoffset.y;
	const int &koff = fineoffset.z;

	int ii, jj, kk;
	sayaka::SelectStaggerIncr(dir, ii, jj, kk);

	const int &ilo = finebox.ilo();
	const int &jlo = finebox.jlo();
	const int &klo = finebox.klo();
	const int &ihi = finebox.ihi();
	const int &jhi = finebox.jhi();
	const int &khi = finebox.khi();

	// (I,j,k) is a face index
	// it can also be treated as the low face of cell(i,j,k)
	for (int k=klo; k<=khi; k++) {
		int kc = IndexMapping::DirFineToCrse<2>(k, koff);
		for (int j=jlo; j<=jhi; j++) {
			int jc = IndexMapping::DirFineToCrse<1>(j, joff);
			for (int i=ilo; i<=ihi; i++) {
				int ic = IndexMapping::DirFineToCrse<0>(i, ioff);
				
				// (ic,jc,kc) should be the cell index of coarse data
				// it is also the low face index of coarse data

				// ijk is the fine index in the face direction
				// if low side, then at low face of coarse cell
				// if high side, then at center face of coarse cell
				int ijk = sayaka::SelectDirIndex(dir, i, j, k);
				int fineside = IndexMapping::FineSide(ijk);
				assert(fineside==0 || fineside==1);

				if (fineside == 0) {
					// on the low face of coarse cell(ic,jc,kc)
					double vl = crsedata(ic,jc,kc,comp);
					finedata(i,j,k,comp) = vl;
				} else {
					// in the middle of coarse cell(ic,jc,kc)
					double vl = crsedata(ic,jc,kc,comp);
					double vh = crsedata(ic+ii,jc+jj,kc+kk,comp);
					finedata(i,j,k,comp) = 0.5 * (vl + vh);
				}
			}
		}
	}
}

void InterpPatch::prolongFaceBlockData_CenteredLinear(int facedir,
	int ifine, DoubleBlockData &finedata,
	int icrse, const DoubleBlockData &crsedata,
	const IndexBox &finebox, const Vector3i &fineoffset,
	int comp,
	int use_limiter) const
{
	assert(0<=facedir && facedir<NDIM);
	assert(var_loc.isFaceVar(facedir));
	assert(finebox.isFaceBox(facedir));

	if (tree_data.numGrow() < 2) {
		LOGPRINTF("%s: NGROW=%d < 2\n", __FUNCTION__, tree_data.numGrow());
		exit(1);
	}

	const int &ioff = fineoffset.x;
	const int &joff = fineoffset.y;
	const int &koff = fineoffset.z;

	const int ii = facedir==0 ? 1 : 0;
	const int jj = facedir==1 ? 1 : 0;
	const int kk = facedir==2 ? 1 : 0;
	const int ii1 = (1 - ii) * XDIM;
	const int jj1 = (1 - jj) * YDIM;
	const int kk1 = (1 - kk) * ZDIM;

	const int &ilo = finebox.ilo();
	const int &jlo = finebox.jlo();
	const int &klo = finebox.klo();
	const int &ihi = finebox.ihi();
	const int &jhi = finebox.jhi();
	const int &khi = finebox.khi();

	// coarse box range
	const int iclo = IndexMapping::DirFineToCrse<0>(ilo, ioff);
	const int ichi = IndexMapping::DirFineToCrse<0>(ihi, ioff);
	const int jclo = IndexMapping::DirFineToCrse<1>(jlo, joff);
	const int jchi = IndexMapping::DirFineToCrse<1>(jhi, joff);
	const int kclo = IndexMapping::DirFineToCrse<2>(klo, koff);
	const int kchi = IndexMapping::DirFineToCrse<2>(khi, koff);

	// 1st sweep, calculate coarse face gradient
	// and set fine data just on that face
	for (int kc=kclo; kc<=kchi; kc++) {
		for (int jc=jclo; jc<=jchi; jc++) {
			for (int ic=iclo; ic<=ichi; ic++) {
				// calculate face gradient
				// gradient in face direction is not calculated
				Vector3d crsegrad;
				crsegrad.setZero();

				for (int dir=0; dir<NDIM; dir++) {
					if (dir != facedir) {
						int iic, jjc, kkc;
						sayaka::SelectDirIncr(dir, iic, jjc, kkc);
						
						double grad = 0;
						if (use_limiter) {
							grad = limited_slope(
								crsedata(ic-iic,jc-jjc,kc-kkc,comp),
								crsedata(ic,jc,kc,comp),
								crsedata(ic+iic,jc+jjc,kc+kkc,comp));
						} else {
							grad = centered_slope(
								crsedata(ic-iic,jc-jjc,kc-kkc,comp),
								crsedata(ic,jc,kc,comp),
								crsedata(ic+iic,jc+jjc,kc+kkc,comp));
						}

						crsegrad(dir) = grad;
					}
				}

				// interpolate 4 (3D) or 2 (2D) fine face data
				for (int kind=-kk1; kind<=kk1; kind+=2) {
					int k = IndexMapping::DirCrseToFine<2>(kc, (kind+1)/2, koff);
					if (k<klo || k>khi) continue;
					double krel = 0.25 * (double) kind;
					
					for (int jind=-jj1; jind<=jj1; jind+=2) {
						int j = IndexMapping::DirCrseToFine<1>(jc, (jind+1)/2, joff);
						if (j<jlo || j>jhi) continue;
						double jrel = 0.25 * (double) jind;

						for (int iind=-ii1; iind<=ii1; iind+=2) {
							int i = IndexMapping::DirCrseToFine<0>(ic, (iind+1)/2, ioff);
							if (i<ilo || i>ihi) continue;
							double irel = 0.25 * (double) iind;

							finedata(i,j,k,comp) = crsedata(ic,jc,kc,comp) +
								irel*crsegrad.x + jrel*crsegrad.y + krel*crsegrad.z;

							assert(SelectDirIndex(facedir, i,j,k) ==
								2 * SelectDirIndex(facedir, ic-ioff,jc-joff,kc-koff));
						}
					}
				}
			}
		}
	}

	// 2nd sweep, set fine data in the middle of two coarse faces
	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				int ijk = sayaka::SelectDirIndex(facedir, i, j, k);
				int side = IndexMapping::FineSide(ijk);
				assert(side==0 || side==1);

				if (side == 1) { // middle face
					const double &vl = finedata(i-ii,j-jj,k-kk,comp);
					const double &vr = finedata(i+ii,j+jj,k+kk,comp);

					finedata(i,j,k,comp) = 0.5 * (vl + vr);
				}
			}
		}
	}
}


void InterpPatch::restrictFaceBlockData_Average(int facedir,
		int icrse, DoubleBlockData &crsedata,
		int ifine, const DoubleBlockData &finedata,
		const IndexBox &crsebox, 
		const Vector3i &fineoffset,
		int comp, int do_average) const
{
	assert(0<=facedir && facedir<NDIM);
	assert(var_loc.isFaceVar(facedir));
	assert(crsebox.isFaceBox(facedir));

	const int &ioff = fineoffset.x;
	const int &joff = fineoffset.y;
	const int &koff = fineoffset.z;

	const int ii = facedir==0 ? 1 : 0;
	const int jj = facedir==1 ? 1 : 0;
	const int kk = facedir==2 ? 1 : 0;
	const int ii1 = (1 - ii) * XDIM;
	const int jj1 = (1 - jj) * YDIM;
	const int kk1 = (1 - kk) * ZDIM;

	const int &ilo = crsebox.ilo();
	const int &ihi = crsebox.ihi();
	const int &jlo = crsebox.jlo();
	const int &jhi = crsebox.jhi();
	const int &klo = crsebox.klo();
	const int &khi = crsebox.khi();

	const double scale = 1.0 / ((XDIM+1)*(YDIM+1)*(ZDIM+1)/2);
	//const double scale = 1.0/ ChildIndex::NumChildOnFace;

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				double val = 0;

				// fine faces
				int count = 0;
				for (int kind=-kk1; kind<=kk1; kind+=2) {
					int kfine = IndexMapping::DirCrseToFine<2>(k, (kind+1)/2, koff);
					for (int jind=-jj1; jind<=jj1; jind+=2) {
						int jfine = IndexMapping::DirCrseToFine<1>(j, (jind+1)/2, joff);
						for (int iind=-ii1; iind<=ii1; iind+=2) {
							int ifine = IndexMapping::DirCrseToFine<0>(i, (iind+1)/2, ioff);
							
							val += finedata(ifine,jfine,kfine,comp);
							count += 1;

							assert(SelectDirIndex(facedir, ifine,jfine,kfine) ==
								2 * SelectDirIndex(facedir, i-ioff,j-joff,k-koff));
						}
					}
				}
				assert(count == ChildIndex::NumChildOnFace);

				if (do_average) {
					crsedata(i,j,k,comp) = val * scale;
				} else {
					crsedata(i,j,k,comp) = val;
				}
			}
		}
	}
}


} // namespace_sayaka


