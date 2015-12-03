

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <memory>

#include "log.h"

#include "SayakaVarLocInterpolator.h"


namespace sayaka
{

// cell -> face(dir)
void VarLocInterpolator::interp_block_cell_to_face(
	int dir,
	DoubleBlockData &face_data,
	const DoubleBlockData &cell_data, 
	const IndexBox &dstbox,
	int dstcomp, int srccomp) const
{
	assert(0<=dir && dir<NDIM);
	// box must be the correct type
	assert(face_data.box().isFaceBox(dir));
	assert(cell_data.box().isCellBox());
	assert(dstbox.isFaceBox(dir));
	// must have enough size of the cell box
	assert(cell_data.box().contains(dstbox));

	int ii, jj, kk;
	sayaka::SelectStaggerIncr(dir, ii, jj, kk);

	const int &ilo = dstbox.ilo();
	const int &ihi = dstbox.ihi();
	const int &jlo = dstbox.jlo();
	const int &jhi = dstbox.jhi();
	const int &klo = dstbox.klo();
	const int &khi = dstbox.khi();

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				// NOTE the staggered data
				// cell(i-1,j) | face(I,j) | cell(i,j)
				const double &vl = cell_data(i-ii,j-jj,k-kk,srccomp);
				const double &vh = cell_data(i,j,k,srccomp);
				face_data(i,j,k,dstcomp) = 0.5 * (vl + vh);
			}
		}
	}
}

// face(dir) -> cell
void VarLocInterpolator::interp_block_face_to_cell(
	int dir,
	DoubleBlockData &cell_data,
	const DoubleBlockData &face_data,
	const IndexBox &dstbox,
	int dstcomp, int srccomp) const
{
	assert(0<=dir && dir<NDIM);
	// check box type
	assert(cell_data.box().isCellBox());
	assert(face_data.box().isFaceBox(dir));
	assert(dstbox.isCellBox());
	//

	int ii, jj, kk;
	sayaka::SelectStaggerIncr(dir, ii, jj, kk);

	const int &ilo = dstbox.ilo();
	const int &ihi = dstbox.ihi();
	const int &jlo = dstbox.jlo();
	const int &jhi = dstbox.jhi();
	const int &klo = dstbox.klo();
	const int &khi = dstbox.khi();

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				// NOTE the staggered data
				// face(I,j) | cell(i,j) | face(I+1,j)
				const double &vl = face_data(i,j,k,srccomp);
				const double &vr = face_data(i+ii,j+jj,k+kk,srccomp);
				cell_data(i,j,k,dstcomp) = 0.5 * (vl + vr);
			}
		}
	}
}

// cell -> node
void VarLocInterpolator::interp_block_cell_to_node(
	DoubleBlockData &node_data,
	const DoubleBlockData &cell_data,
	const IndexBox &dstbox,
	int dstcomp, int srccomp) const
{
	// box type
	assert(node_data.box().isNodeBox());
	assert(cell_data.box().isCellBox());
	assert(dstbox.isNodeBox());

	const int &ilo = dstbox.ilo();
	const int &ihi = dstbox.ihi();
	const int &jlo = dstbox.jlo();
	const int &jhi = dstbox.jhi();
	const int &klo = dstbox.klo();
	const int &khi = dstbox.khi();

	const double scale = 1.0 / ((XDIM+1)*(YDIM+1)*(ZDIM+1));

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				double val = 0;
				// loop surrounding cells
				for (int kc=k-ZDIM; kc<=k; kc++) {
					for (int jc=j-YDIM; jc<=j; jc++) {
						for (int ic=i-XDIM; ic<=i; ic++) {
							val += cell_data(ic,jc,kc,srccomp);
						}
					}
				}
				node_data(i,j,k,dstcomp) = val * scale;
			}
		}
	}
}

// node -> cell
void VarLocInterpolator::interp_block_node_to_cell(
	DoubleBlockData &cell_data,
	const DoubleBlockData &node_data,
	const IndexBox &dstbox,
	int dstcomp, int srccomp) const
{
	// box type
	assert(cell_data.box().isCellBox());
	assert(node_data.box().isNodeBox());
	assert(dstbox.isCellBox());

	const int &ilo = dstbox.ilo();
	const int &ihi = dstbox.ihi();
	const int &jlo = dstbox.jlo();
	const int &jhi = dstbox.jhi();
	const int &klo = dstbox.klo();
	const int &khi = dstbox.khi();

	const double scale = 1.0 / ((XDIM+1)*(YDIM+1)*(ZDIM+1));

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				double val = 0;
				// loop surrounding nodes
				for (int kn=k; kn<=k+ZDIM; kn++) {
					for (int jn=j; jn<=j+YDIM; jn++) {
						for (int in=i; in<=i+XDIM; in++) {
							val += node_data(in,jn,kn,srccomp);
						}
					}
				}
				cell_data(i,j,k,dstcomp) = val * scale;
			}
		}
	}
}

// face(dir) -> node
void VarLocInterpolator::interp_block_face_to_node(int dir,
	DoubleBlockData &node_data,
	const DoubleBlockData &face_data,
	const IndexBox &dstbox,
	int dstcomp, int srccomp) const
{
	assert(0<=dir && dir<NDIM);
	// box type
	assert(node_data.box().isNodeBox());
	assert(face_data.box().isFaceBox(dir));
	assert(dstbox.isNodeBox());

	const int ii = dir==0 ? 0 : 1;
	const int jj = dir==1 ? 0 : 1;
	const int kk = dir==2 ? 0 : 1;

	const int &ilo = dstbox.ilo();
	const int &ihi = dstbox.ihi();
	const int &jlo = dstbox.jlo();
	const int &jhi = dstbox.jhi();
	const int &klo = dstbox.klo();
	const int &khi = dstbox.khi();

	const double scale = 1.0 / ((XDIM+1)*(YDIM+1)*(ZDIM+1)/2);

	for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
			for (int i=ilo; i<=ihi; i++) {
				double val = 0;
				// loop surrounding faces
				for (int kface=k-kk*ZDIM; kface<=k; kface++) {
					for (int jface=j-jj*YDIM; jface<=j; jface++) {
						for (int iface=i-ii*XDIM; iface<=i; iface++) {
							val += face_data(iface,jface,kface,srccomp);
						}
					}
				}
				node_data(i,j,k,dstcomp) = val * scale;
			}
		}
	}
}

//
void VarLocInterpolator::interp_cell_to_face(int dir,
	TreeData &face_data, const TreeData &cell_data,
	int dstcomp, int srccomp, int ncomp, int ngrow) const
{
	assert(m_tree);
	assert(0<=dir && dir<NDIM);
	//
	assert(face_data.isFaceData(dir));
	assert(cell_data.isCellData());
	//
	assert(0<=dstcomp && dstcomp+ncomp<=face_data.numComp());
	assert(0<=srccomp && srccomp+ncomp<=cell_data.numComp());
	//
	assert(face_data.numGrow() >= ngrow);
	assert(cell_data.numGrow() > ngrow);

	//
	IndexBox dstbox = face_data.validBox();
	dstbox.extend(ngrow);

	for (int i=0; i<m_tree->numBlocks; i++) {
		for (int icomp=0; icomp<ncomp; icomp++) {
			interp_block_cell_to_face(dir,
				face_data[i], cell_data[i],
				dstbox, dstcomp+icomp, srccomp+icomp);
		}
	}
}

//
void VarLocInterpolator::interp_face_to_cell(int dir,
	TreeData &cell_data, const TreeData &face_data,
	int dstcomp, int srccomp, int ncomp, int ngrow) const
{
	assert(m_tree);
	assert(0<=dir && dir<NDIM);
	//
	assert(cell_data.isCellData());
	assert(face_data.isFaceData(dir));
	//
	assert(0<=dstcomp && dstcomp+ncomp<=cell_data.numComp());
	assert(0<=srccomp && srccomp+ncomp<=face_data.numComp());
	//
	assert(cell_data.numGrow() >= ngrow);
	assert(face_data.numGrow() > ngrow);

	// target box
	IndexBox dstbox = cell_data.validBox();
	dstbox.extend(ngrow);

	for (int i=0; i<m_tree->numBlocks; i++) {
		for (int icomp=0; icomp<ncomp; icomp++) {
			interp_block_face_to_cell(dir,
				cell_data[i], face_data[i],
				dstbox, dstcomp+icomp, srccomp+icomp);
		}
	}
}

//
void VarLocInterpolator::interp_cell_to_node(
	TreeData &node_data, const TreeData &cell_data,
	int dstcomp, int srccomp, int ncomp, int ngrow) const
{
	assert(m_tree);
	//
	assert(node_data.isNodeData());
	assert(cell_data.isCellData());
	//
	assert(0<=dstcomp && dstcomp+ncomp<=node_data.numComp());
	assert(0<=srccomp && srccomp+ncomp<=cell_data.numComp());
	//
	assert(node_data.numGrow() >= ngrow);
	assert(cell_data.numGrow() > ngrow);

	// target box
	IndexBox dstbox = node_data.validBox();
	dstbox.extend(ngrow);

	for (int i=0; i<m_tree->numBlocks; i++) {
		for (int icomp=0; icomp<ncomp; icomp++) {
			interp_block_cell_to_node(
				node_data[i], cell_data[i],
				dstbox, dstcomp+icomp, srccomp+icomp);
		}
	}
}

} // namespace_sayaka


