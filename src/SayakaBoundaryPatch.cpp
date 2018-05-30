

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <memory>

#include "SayakaBoundaryPatch.h"


SAYAKA_NS_BEGIN;

BCPatch::BCPatch()
	: m_tree(nullptr)
	, m_varloc(VARLOC_CELL)
	, m_ncomp(0), m_ngrow(0)
{
	// nothing
}

BCPatch::BCPatch(
	AmrTree & tree_in, 
	TreeData & data_in, 
	VariableLocation varloc_in)
//: super_type(tree_in, data_in, varloc_in),
	: m_tree(&tree_in), m_varloc(varloc_in)
	, m_ncomp(data_in.numComp()), m_ngrow(data_in.numGrow())
	, bndryRegs(data_in.numComp())
{
	calcRanges();
}

void BCPatch::define(
	AmrTree & tree, 
	VariableLocation varloc, 
	int ncomp, int ngrow)
{
	m_tree = &tree;
	m_varloc = varloc;
	m_ncomp = ncomp;
	m_ngrow = ngrow;

	bndryRegs.resize(ncomp);

	//
	calcRanges();
}



// 
void BCPatch::calcRanges() {
	const IndexBox &validBox = m_tree->validBlockCellBox();
	const Vector3i vnum = validBox.size();

	//const int nlayer = tree_data.numGrow();
	const int nlayer = m_ngrow;
	const IndexBox layerBox = IndexBox(validBox).extend(nlayer);

	for (int kk=-ZDIM; kk<=ZDIM; kk++) {
		for (int jj=-YDIM; jj<=YDIM; jj++) {
			for (int ii=-XDIM; ii<=XDIM; ii++) {
				Vector3i voff;
				voff.x = ii * vnum.x;
				voff.y = jj * vnum.y;
				voff.z = kk * vnum.z;

				// boundary ranges of central block
				IndexBox bndryBox = validBox;
				bndryBox.shift(voff);
				bndryBox = IndexBox::Intersection(bndryBox, layerBox);
				assert(bndryBox.isValid());

				//
				bndryBoxes[ii+1][jj+1][kk+1] = bndryBox;
			}
		}
	}
}



void BCPatch::fillBlockBC(
	int physbc, 
	int isweep, int sweepdir,
	int xbndry, int ybndry, int zbndry,
	int iblock, DoubleBlockData &data, 
	int scomp, int ncomp) const
{
	if (m_varloc == VARLOC_CELL) {
		for (int comp=scomp; comp<scomp+ncomp; comp++) {

			// find BC value from BC register
			const BCRegister &bcreg = boundaryRegister(comp);
			int mathbc = bcreg.getBCMap(xbndry,ybndry,zbndry,physbc);
			double bcval = bcreg.getBCVal(xbndry, ybndry, zbndry, physbc);

			fillBasicCellBlockBC(mathbc, bcval, 
				isweep, sweepdir,
				xbndry, ybndry, zbndry,
				iblock, data, comp);
		}
	}
#if 0
	else if (m_varloc.isFaceVar()) {
		for (int comp = scomp; comp < scomp + ncomp; comp++) {
			// find BC value from BC register
			const BCRegister &bcreg = boundaryRegister(comp);
			int mathbc = bcreg.getBCMap(xbndry, ybndry, zbndry, physbc);
			double bcval = bcreg.getBCVal(xbndry, ybndry, zbndry, physbc);

			if (m_varloc.isFacexVar()) {
				fillBasicFacexBlockBC(mathbc, bcval,
					xbndry, ybndry, zbndry,
					iblock, data, comp);
			}
			else if (m_varloc.isFaceyVar()) {
				fillBasicFaceyBlockBC(mathbc, bcval,
					xbndry, ybndry, zbndry,
					iblock, data, comp);
			}
			else {
				fillBasicFacezBlockBC(mathbc, bcval,
					xbndry, ybndry, zbndry,
					iblock, data, comp);
			}
		}
	}
#endif
	
	else {
		LOGPRINTF("%s: cell-centered data only!\n", __FUNCTION__);
	}
}

void BCPatch::fillBasicCellBlockBC(
	int bctype, double bcval,
	int isweep, int sweepdir,
	int xbndry, int ybndry, int zbndry,
	int iblock, DoubleBlockData &blockdata, int comp) const
{
	assert(0<=comp && comp<blockdata.ncomp());
	assert(-XDIM<=xbndry && xbndry<=XDIM);
	assert(-YDIM<=ybndry && ybndry<=YDIM);
	assert(-ZDIM<=zbndry && zbndry<=ZDIM);
	assert(xbndry!=0 || ybndry!=0 || zbndry!=0);
	assert(0<=isweep && isweep<NDIM);
	assert(0<=sweepdir && sweepdir<NDIM);

	const SurroundingBlocks &surr = m_tree->surrBlocks[iblock];
	assert(surr.isExternalBoundarySurrounding(xbndry,ybndry,zbndry));

	if (m_varloc != VARLOC_CELL) {
		LOGPRINTF("%s: only for cell centered\n", __FUNCTION__);
		exit(1);
	}

	const IndexBox &validbox = m_tree->validBlockCellBox();
	const int nx = validbox.size(0);
	const int ny = validbox.size(1);
	const int nz = validbox.size(2);

	// mapping between ghost cell and corresponding real cell
	// (ig,jg,kg) ~ (isign*ig+irel,jsign*jg+jrel,ksign*kg+krel)
	//const int isign = xbndry==0 ? 1 : -1;
	//const int jsign = ybndry==0 ? 1 : -1;
	//const int ksign = zbndry==0 ? 1 : -1;
	//const int ioff = xbndry==0 ? 0 : xbndry==-1 ? -1 : 2*nx-1;
	//const int joff = ybndry==0 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
	//const int koff = zbndry==0 ? 0 : zbndry==-1 ? -1 : 2*nz-1;
	int isign = sweepdir==0 ? -1 : 1;
	int jsign = sweepdir==1 ? -1 : 1;
	int ksign = sweepdir==2 ? -1 : 1;
	int ioff = sweepdir!=0 ? 0 : xbndry==-1 ? -1 : 2*nx-1;
	int joff = sweepdir!=1 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
	int koff = sweepdir!=2 ? 0 : zbndry==-1 ? -1 : 2*nz-1;

	if (surr.bcref(xbndry,ybndry,zbndry) >= 0) {
		SurroundingIndex bcref = surr.bcref(xbndry,ybndry,zbndry);
		int xref = bcref.ipos();
		int yref = bcref.jpos();
		int zref = bcref.kpos();

		int bcdir = -1;
		if (xbndry != xref) bcdir = 0;
		if (ybndry != yref) bcdir = 1;
		if (zbndry != zref) bcdir = 2;
		assert(0<=bcdir && bcdir<NDIM);

		isign = bcdir==0 ? -1 : 1;
		jsign = bcdir==1 ? -1 : 1;
		ksign = bcdir==2 ? -1 : 1;
		ioff = bcdir!=0 ? 0 : xbndry==-1 ? -1 : 2*nx-1;
		joff = bcdir!=1 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
		koff = bcdir!=2 ? 0 : zbndry==-1 ? -1 : 2*nz-1;
	}

	const IndexBox &bndrybox = bndryBox(xbndry, ybndry, zbndry);
	const int ilo = bndrybox.ilo();
	const int jlo = bndrybox.jlo();
	const int klo = bndrybox.klo();
	const int ihi = bndrybox.ihi();
	const int jhi = bndrybox.jhi();
	const int khi = bndrybox.khi();

	for (int k=klo; k<=khi; k++) {
		int kref = ksign*k + koff;
		for (int j=jlo; j<=jhi; j++) {
			int jref = jsign*j + joff;
			for (int i=ilo; i<=ihi; i++) {
				int iref = isign*i + ioff;

				const double &vref = blockdata(iref,jref,kref,comp);

				if (bctype == BCType::Neumann) {
					// homogeneous Neumann only
					blockdata(i,j,k,comp) = vref;
				} else if (bctype == BCType::Dirichlet) {
					blockdata(i,j,k,comp) = 2.0*bcval - vref;
				} else if (bctype == BCType::SimpleFill) {
					blockdata(i,j,k,comp) = bcval;
				}
			}
		}
	}
} // boundaryconditionpatch_fillbasiccellblockbc


void BCPatch::fillBasicFacexBlockBC(
	int bctype, double bcval,
	int xbndry, int ybndry, int zbndry,
	int iblock, DoubleBlockData &blockdata, int comp) const
{
	assert(0<=comp && comp<blockdata.ncomp());
	assert(-XDIM<=xbndry && xbndry<=XDIM);
	assert(-YDIM<=ybndry && ybndry<=YDIM);
	assert(-ZDIM<=zbndry && zbndry<=ZDIM);
	assert(xbndry!=0 || ybndry!=0 || zbndry!=0);

	if (m_varloc!=VARLOC_FACE_X) {
		LOGPRINTF("%s: only for x-face data\n", __FUNCTION__);
		exit(1);
	}

	const IndexBox &validbox = m_tree->validBlockCellBox();
	const int nx = validbox.size(0);
	const int ny = validbox.size(1);
	const int nz = validbox.size(2);

	// mapping between ghost cell and corresponding real cell
	// (ig,jg,kg) ~ (isign*ig+irel,jsign*jg+jrel,ksign*kg+krel)
	//const int isign = xbndry==0 ? 1 : -1;
	const int jsign = ybndry==0 ? 1 : -1;
	const int ksign = zbndry==0 ? 1 : -1;

	//const int ioff = xbndry==0 ? 0 : xbndry==-1 ? 0 : 2*nx;
	const int joff = ybndry==0 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
	const int koff = zbndry==0 ? 0 : zbndry==-1 ? -1 : 2*nz-1;

	const IndexBox &bndrybox = bndryBox(xbndry, ybndry, zbndry);
	const int ilo = xbndry==0 ? bndrybox.ilo()+1 : bndrybox.ilo();
	const int jlo = bndrybox.jlo();
	const int klo = bndrybox.klo();
	const int ihi = xbndry==0 ? bndrybox.ihi() : bndrybox.ihi()+1;
	const int jhi = bndrybox.jhi();
	const int khi = bndrybox.khi();

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				if (xbndry == 0) {
					int iref = i;
					int jref = jsign*j + joff;
					int kref = ksign*k + koff;
					const double &vref = blockdata(iref,jref,kref,comp);
					if (bctype == BCType::Neumann) {
						blockdata(i,j,k,comp) = vref;
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = 2.0*bcval - vref;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				} else {
					if (bctype == BCType::Neumann) {
						int iref = xbndry==-1 ? 0 : nx;
						int jref = jsign*j + joff;
						int kref = ksign*k + koff;
						blockdata(i,j,k,comp) = blockdata(iref,jref,kref,comp);
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = bcval;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				}
			}
		}
	}
}

void BCPatch::fillBasicFaceyBlockBC(
	int bctype, double bcval,
	int xbndry, int ybndry, int zbndry,
	int iblock, DoubleBlockData &blockdata, int comp) const
{
	assert(0<=comp && comp<blockdata.ncomp());
	assert(-XDIM<=xbndry && xbndry<=XDIM);
	assert(-YDIM<=ybndry && ybndry<=YDIM);
	assert(-ZDIM<=zbndry && zbndry<=ZDIM);
	assert(xbndry!=0 || ybndry!=0 || zbndry!=0);

	if (m_varloc!=VARLOC_FACE_Y) {
		LOGPRINTF("%s: only for y-face data\n", __FUNCTION__);
		exit(1);
	}
	static_assert(NDIM >= 2, "NDIM invalid");

	const IndexBox &validbox = m_tree->validBlockCellBox();
	const int nx = validbox.size(0);
	const int ny = validbox.size(1);
	const int nz = validbox.size(2);

	// mapping between ghost cell and corresponding real cell
	// (ig,jg,kg) ~ (isign*ig+irel,jsign*jg+jrel,ksign*kg+krel)
	const int isign = xbndry==0 ? 1 : -1;
	//const int jsign = ybndry==0 ? 1 : -1;
	const int ksign = zbndry==0 ? 1 : -1;

	const int ioff = xbndry==0 ? 0 : xbndry==-1 ? -1 : 2*nx-1;
	//const int joff = ybndry==0 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
	const int koff = zbndry==0 ? 0 : zbndry==-1 ? -1 : 2*nz-1;

	const IndexBox &bndrybox = bndryBox(xbndry, ybndry, zbndry);
	const int ilo = bndrybox.ilo();
	const int jlo = ybndry==0 ? bndrybox.jlo()+1 : bndrybox.jlo();
	const int klo = bndrybox.klo();
	const int ihi = bndrybox.ihi();
	const int jhi = ybndry==0 ? bndrybox.jhi() : bndrybox.jhi()+1;
	const int khi = bndrybox.khi();

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				if (ybndry == 0) {
					int iref = isign*i + ioff;
					int jref = j;
					int kref = ksign*k + koff;
					const double &vref = blockdata(iref,jref,kref,comp);
					if (bctype == BCType::Neumann) {
						blockdata(i,j,k,comp) = vref;
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = 2.0*bcval - vref;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				} else {
					if (bctype == BCType::Neumann) {
						int iref = isign*i + ioff;
						int jref = ybndry==-1 ? 0 : ny;
						int kref = ksign*k + koff;
						blockdata(i,j,k,comp) = blockdata(iref,jref,kref,comp);
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = bcval;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				}
			}
		}
	}
}

void BCPatch::fillBasicFacezBlockBC(
	int bctype, double bcval,
	int xbndry, int ybndry, int zbndry,
	int iblock, DoubleBlockData &blockdata, int comp) const
{
	assert(0<=comp && comp<blockdata.ncomp());
	assert(-XDIM<=xbndry && xbndry<=XDIM);
	assert(-YDIM<=ybndry && ybndry<=YDIM);
	assert(-ZDIM<=zbndry && zbndry<=ZDIM);
	assert(xbndry!=0 || ybndry!=0 || zbndry!=0);

	if (m_varloc!=VARLOC_FACE_Z) {
		LOGPRINTF("%s: only for z-face data\n", __FUNCTION__);
		exit(1);
	}
	static_assert(NDIM == 3, "NDIM invalid");

	const IndexBox &validbox = m_tree->validBlockCellBox();
	const int nx = validbox.size(0);
	const int ny = validbox.size(1);
	const int nz = validbox.size(2);

	// mapping between ghost cell and corresponding real cell
	// (ig,jg,kg) ~ (isign*ig+irel,jsign*jg+jrel,ksign*kg+krel)
	const int isign = xbndry==0 ? 1 : -1;
	const int jsign = ybndry==0 ? 1 : -1;
	//const int ksign = zbndry==0 ? 1 : -1;

	const int ioff = xbndry==0 ? 0 : xbndry==-1 ? -1 : 2*nx-1;
	const int joff = ybndry==0 ? 0 : ybndry==-1 ? -1 : 2*ny-1;
	//const int koff = zbndry==0 ? 0 : zbndry==-1 ? -1 : 2*nz-1;

	const IndexBox &bndrybox = bndryBox(xbndry, ybndry, zbndry);
	const int ilo = bndrybox.ilo();
	const int jlo = bndrybox.jlo();
	const int klo = zbndry==0 ? bndrybox.klo()+1 : bndrybox.klo();
	const int ihi = bndrybox.ihi();
	const int jhi = bndrybox.jhi();
	const int khi = zbndry==0 ? bndrybox.khi() : bndrybox.khi()+1;

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				if (zbndry == 0) {
					int iref = isign*i + ioff;
					int jref = jsign*j + joff;
					int kref = k;
					const double &vref = blockdata(iref,jref,kref,comp);
					if (bctype == BCType::Neumann) {
						blockdata(i,j,k,comp) = vref;
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = 2.0*bcval - vref;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				} else {
					if (bctype == BCType::Neumann) {
						int iref = isign*i + ioff;
						int jref = jsign*j + joff;
						int kref = zbndry==-1 ? 0 : nz;
						blockdata(i,j,k,comp) = blockdata(iref,jref,kref,comp);
					} else if (bctype == BCType::Dirichlet) {
						blockdata(i,j,k,comp) = bcval;
					} else if (bctype == BCType::SimpleFill) {
						blockdata(i,j,k,comp) = bcval;
					}
				}
			}
		}
	}
}




SAYAKA_NS_END;




