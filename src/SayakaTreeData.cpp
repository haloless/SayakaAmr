
#include "SayakaTreeData.h"



namespace sayaka
{

//
TreeData* TreeData::CreateCellData(const AmrTree &tree, int ncomp, int ngrow) {
	const IndexBox &validbox = tree.validBlockCellBox();
	
	IndexBox grownbox = validbox;
	grownbox.extend(ngrow);

	TreeData *pdata = new TreeData(const_cast<AmrTree&>(tree), validbox, grownbox, ngrow, ncomp);
	return pdata;
}
//
TreeData* TreeData::CreateFaceData(const AmrTree &tree, int dir, int ncomp, int ngrow) {
	const IndexBox &cellbox = tree.validBlockCellBox();

	IndexBox facebox = cellbox;
	facebox.staggerInDir(dir);

	IndexBox grownbox = facebox;
	grownbox.extend(ngrow);

	TreeData *pdata = new TreeData(const_cast<AmrTree&>(tree), facebox, grownbox, ngrow, ncomp);
	return pdata;
}
//
TreeData* TreeData::CreateNodeData(const AmrTree &tree, int ncomp, int ngrow) {
	const IndexBox &cellbox = tree.validBlockCellBox();

	IndexBox nodebox = cellbox;
	nodebox.staggerAll();

	IndexBox grownbox = nodebox;
	grownbox.extend(ngrow);

	TreeData *pdata = new TreeData(const_cast<AmrTree&>(tree), nodebox, grownbox, ngrow, ncomp);
	return pdata;
}
//
std::vector<TreeData*> TreeData::CreateFaceDataPArray(const AmrTree &tree, int ncomp, int ngrow) {
	std::vector<TreeData*> parray(NDIM, static_cast<TreeData*>(NULL));
	for (int dir=0; dir<NDIM; dir++) {
		TreeData *pdata = CreateFaceData(tree, dir, ncomp, ngrow);
		parray[dir] = pdata;
	}
	return parray;
}
//
void TreeData::ReleaseDataPArray(std::vector<TreeData*> &parray, bool clear_array) {
	for (int i=0; i<parray.size(); i++) {
		TreeData *pdata = parray[i];
		if (pdata) {
			delete pdata;
			parray[i] = NULL;
		}
	}
	if (clear_array) {
		parray.clear();
	}
}


void TreeData::Copy(TreeData &dst, const TreeData& src,
	int dcomp, int scomp, int ncomp, int ngrow)
{
	assert(dst.validBox() == src.validBox());
	assert(0<=dcomp && dcomp+ncomp<=dst.numComp());
	assert(0<=scomp && scomp+ncomp<=src.numComp());
	assert(0<=ngrow && ngrow<=dst.numGrow());
	assert(0<=ngrow && ngrow<=src.numGrow());

	const AmrTree &tree = dst.tree();
	const IndexBox &vbox = dst.validBox();
	const IndexBox &box = IndexBox::Extend(vbox, ngrow);

	for (int i=0; i<tree.numBlocks; i++) {
		DoubleGridDataUtil.CopyData(
			dst[i], src[i], 
			dcomp, scomp, ncomp,
			box);
	}
}
void TreeData::SetValue(TreeData &dst, double val,
	int dcomp, int ncomp, int ngrow)
{
	assert(0<=dcomp && dcomp+ncomp<=dst.numComp());
	assert(0<=ngrow && ngrow<=dst.numGrow());

	const AmrTree &tree = dst.tree();
	const IndexBox &vbox = dst.validBox();
	const IndexBox &box = IndexBox::Extend(vbox, ngrow);

	for (int i=0; i<tree.numBlocks; i++) {
		dst[i].setValue(dcomp, ncomp, box, val);
	}
}


//
void TreeData::Copy(TreeData &dst, const TreeData &src, 
	int dcomp, int scomp, int ncomp, int ngrow,
	const int iblocks[], const int ibegin, const int iend)
{
	assert(dst.validBox() == src.validBox());
	assert(0<=dcomp && dcomp+ncomp<=dst.numComp());
	assert(0<=scomp && scomp+ncomp<=src.numComp());
	assert(0<=ngrow && ngrow<=dst.numGrow());
	assert(0<=ngrow && ngrow<=src.numGrow());

	const AmrTree &tree = dst.tree();
	const IndexBox &vbox = dst.validBox();
	const IndexBox &box = IndexBox::Extend(vbox, ngrow);

	for (int igrid=ibegin; igrid<iend; igrid++) {
		const int iblock = iblocks[igrid];
		assert(iblock >= 0);

		DoubleGridDataUtil.CopyData(
			dst[iblock], src[iblock], 
			dcomp, scomp, ncomp,
			box);
	}
}
//
void TreeData::SetValue(TreeData &dst, double val,
	int dcomp, int ncomp, int ngrow,
	const int iblocks[], const int ibegin, const int iend)
{
	assert(0<=dcomp && dcomp+ncomp<=dst.numComp());
	assert(0<=ngrow && ngrow<=dst.numGrow());

	const AmrTree &tree = dst.tree();
	const IndexBox &vbox = dst.validBox();
	const IndexBox &box = IndexBox::Extend(vbox, ngrow);

	for (int igrid=ibegin; igrid<iend; igrid++) {
		const int iblock = iblocks[igrid];
		assert(iblock >= 0);

		dst[iblock].setValue(dcomp, ncomp, box, val);
	}	
}
//
void TreeData::AddEqual(TreeData &dst, const TreeData &src, double coef,
	int dcomp, int scomp, int ncomp, int ngrow,
	const int iblocks[], const int ibegin, const int iend)
{
	assert(dst.validBox() == src.validBox());
	assert(0<=dcomp && dcomp+ncomp<=dst.numComp());
	assert(0<=scomp && scomp+ncomp<=src.numComp());
	assert(0<=ngrow && ngrow<=dst.numGrow());
	assert(0<=ngrow && ngrow<=src.numGrow());

	const AmrTree &tree = dst.tree();
	const IndexBox &vbox = dst.validBox();
	const IndexBox &box = IndexBox::Extend(vbox, ngrow);

	for (int igrid=ibegin; igrid<iend; igrid++) {
		const int iblock = iblocks[igrid];
		assert(iblock >= 0);

		DoubleGridDataUtil.AddEqual(
			dst[iblock], src[iblock], coef, 
			dcomp, scomp, ncomp,
			box);
	}
}


int TreeData::writeTreeLeafBlockData(const char *filename, int step, double time) const {
	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}

	FILE *fp = fopen(filename, "w");
	if (!fp) {
		LOGPRINTF("%s: failed to open %s\n", __FUNCTION__, filename);
		return 1;
	}

	const int ncomp = numComp();

	//
	fprintf(fp, "TITLE=\"Amr Data\"\n");
	//
	fprintf(fp, "FILETYPE=SOLUTION\n");
	//
	fprintf(fp, "VARIABLES=");
	for (int comp=0; comp<ncomp; comp++) {
		fprintf(fp, "\"Comp%d\" ", comp);
	}
	fprintf(fp, "\n");

	const IndexBox &validBox = this->validBox();
	const Vector3i vnb = this->tree().validBlockCellBox().size();

	for (int i=0; i<blockNum(); i++) {
		const AmrTreeNode &block = (*m_tree)[i];
		const DoubleBlockData &blockdata = (*this)[i];

		if (m_tree->blocks[i].isLeaf()) {
			//double xlo = block.boundBox.vlo.x;
			//double ylo = block.boundBox.vlo.y;
			//double zlo = block.boundBox.vlo.z;
			//double dx = block.blockLength.x / vnb.x;
			//double dy = block.blockLength.y / vnb.y;
			//double dz = block.blockLength.z / vnb.z;

			fprintf(fp, "ZONE\n");
			fprintf(fp, "I=%d, J=%d\n", vnb.x+1, vnb.y+1);
			fprintf(fp, "DATAPACKING=BLOCK\n");
			fprintf(fp, "VARLOCATION=([1");
			for (int comp=1; comp<ncomp; comp++) {
				fprintf(fp, ",%d", comp+1);
			}
			if (validBox.isCellBox()) {
				fprintf(fp, "]=CELLCENTERED)\n");
			} else if (validBox.isNodeBox()) {
				fprintf(fp, "]=NODAL)\n");
			}
			fprintf(fp, "STRANDID=%d, SOLUTIONTIME=%e\n", step, time);

			for (int comp=0; comp<ncomp; comp++) {
				for (int kk=validBox.klo(); kk<=validBox.khi(); kk++) {
					for (int jj=validBox.jlo(); jj<=validBox.jhi(); jj++) {
						for (int ii=validBox.ilo(); ii<=validBox.ihi(); ii++) {
							fprintf(fp, "%e ", blockdata(ii,jj,kk,comp));
						}
					}
				}
				fprintf(fp, "\n");
			}
		}
	}

	fclose(fp);

	LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}


int TreeData::writeBlockData(int iblock, 
	const char *filename, int step, double time) const 
{
	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}

	FILE *fp = fopen(filename, "w");
	if (!fp) {
		LOGPRINTF("%s: failed to open %s\n", __FUNCTION__, filename);
		return 1;
	}

	const int ncomp = numComp();
	const int ngrow = numGrow();

	//
	fprintf(fp, "TITLE=\"Amr Block Data\"\n");
	//
	fprintf(fp, "FILETYPE=FULL\n");
	//
	fprintf(fp, "VARIABLES=\"X\" \"Y\" ");
	for (int comp=0; comp<ncomp; comp++) {
		fprintf(fp, "\"Comp%d\" ", comp);
	}
	fprintf(fp, "\n");

	const IndexBox &validBox = this->validBox();
	const IndexBox &grownBox = this->indexBox();
	
	const int ilo = grownBox.ilo();
	const int jlo = grownBox.jlo();
	const int klo = grownBox.klo();
	const int ihi = grownBox.ihi();
	const int jhi = grownBox.jhi();
	const int khi = grownBox.khi();

	const AmrTreeNode &block = (*m_tree)[iblock];
	const DoubleBlockData &blockdata = (*this)[iblock];

	double xlo = block.boundBox.vlo.x;
	double ylo = block.boundBox.vlo.y;
	double zlo = block.boundBox.vlo.z;
	const Vector3i vnb = validBox.size();
	double dx = block.blockLength.x / vnb.x;
	double dy = block.blockLength.y / vnb.y;
	double dz = block.blockLength.z / vnb.z;

	fprintf(fp, "ZONE T=\"%d\"\n", iblock);
	fprintf(fp, "I=%d, J=%d\n", ihi-ilo+2, jhi-jlo+2);
	fprintf(fp, "DATAPACKING=BLOCK\n");
	fprintf(fp, "VARLOCATION=([%d", NDIM+1);
	for (int comp=1; comp<ncomp; comp++) {
		fprintf(fp, ",%d", comp+1+NDIM);
	}
	fprintf(fp, "]=CELLCENTERED)\n");
	fprintf(fp, "STRANDID=%d, SOLUTIONTIME=%e\n", step, time);

	// x
	for (int kk=klo; kk<=khi+ZDIM; kk++) {
		for (int jj=jlo; jj<=jhi+YDIM; jj++) {
			for (int ii=ilo; ii<=ihi+XDIM; ii++) {
				double xx = xlo + dx * ii;
				double yy = ylo + dy * jj;
				double zz = zlo + dz * kk;
				fprintf(fp, "%lf ", xx);
			}
		}
	}
	fprintf(fp, "\n");
	// y
	for (int kk=klo; kk<=khi+ZDIM; kk++) {
		for (int jj=jlo; jj<=jhi+YDIM; jj++) {
			for (int ii=ilo; ii<=ihi+XDIM; ii++) {
				double xx = xlo + dx * ii;
				double yy = ylo + dy * jj;
				double zz = zlo + dz * kk;
				fprintf(fp, "%lf ", yy);
			}
		}
	}
	fprintf(fp, "\n");

	// data
	for (int comp=0; comp<ncomp; comp++) {
		for (int kk=klo; kk<=khi; kk++) {
			for (int jj=jlo; jj<=jhi; jj++) {
				for (int ii=ilo; ii<=ihi; ii++) {
					fprintf(fp, "%lf ", blockdata(ii,jj,kk,comp));
				}
			}
		}
		fprintf(fp, "\n");
	}


	fclose(fp);

	// LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}

int TreeData::writeTreeLevelBlockData(int ilevel, 
	const char *filename, int step, double time) const
{
	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}

	FILE *fp = fopen(filename, "w");
	if (!fp) {
		LOGPRINTF("%s: failed to open %s\n", __FUNCTION__, filename);
		return 1;
	}

	// header
	fprintf(fp, "TITLE=\"Amr Level %d\"\n", ilevel);
	//
	fprintf(fp, "FILETYPE=FULL\n");
	//
	fprintf(fp, "VARIABLES=");
	for (int dir=0; dir<NDIM; dir++) {
		const char* axis[MAX_DIM] = { "X", "Y", "Z" };
		fprintf(fp, "\"%s\" ", axis[dir]);
	}
	for (int comp=0; comp<this->numComp(); comp++) {
		fprintf(fp, "\"Comp%d\" ", comp);
	}
	fprintf(fp, "\n");

	for (int iblock=0; iblock<this->blockNum(); iblock++) {
		if (this->block(iblock).getLevel() == ilevel) {
			int writeGrid = 1;
			writeBlockZone(fp, iblock, step, time, writeGrid);
		}
	}


	fclose(fp);

	LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}

void TreeData::writeBlockZone(FILE *fp, 
	int iblock, int step, double time, 
	int writeGrid) const
{
	assert(fp);

	if (NDIM != 2) {
		LOGPRINTF("%s: 2D only\n", __FUNCTION__);
		exit(1);
	}

	const IndexBox &cellBox = this->tree().validBlockCellBox();
	const Vector3i ncell = cellBox.size();
	const IndexBox &validBox = this->validBox();

	const AmrTreeNode &block = (*m_tree)[iblock];
	const DoubleBlockData &blockdata = (*this)[iblock];

	//double xlo = block.boundBox.vlo.x;
	//double ylo = block.boundBox.vlo.y;
	//double zlo = block.boundBox.vlo.z;
	//double dx = block.blockLength.x / vnb.x;
	//double dy = block.blockLength.y / vnb.y;
	//double dz = block.blockLength.z / vnb.z;

	const int ncomp = this->numComp();

	fprintf(fp, "ZONE\n");
	fprintf(fp, "I=%d, J=%d\n", ncell.x+1, ncell.y+1);
	fprintf(fp, "DATAPACKING=BLOCK\n");

	fprintf(fp, "VARLOCATION=([");
	{
		const int offset = writeGrid ? NDIM+1 : 1;
		for (int comp=0; comp<ncomp; comp++) {
			if (comp > 0) fprintf(fp, ",");
			fprintf(fp, "%d", comp+offset);
		}
		if (validBox.isCellBox()) {
			fprintf(fp, "]=CELLCENTERED)\n");
		} else if (validBox.isNodeBox()) {
			fprintf(fp, "]=NODAL)\n");
		}
	}
	fprintf(fp, "STRANDID=%d, SOLUTIONTIME=%e\n", step, time);

	if (writeGrid) {
		const Vector3d &xlo = block.boundBox.lo();
		const Vector3d dx = m_tree->getBlockCellSize(iblock);

		for (int dir=0; dir<NDIM; dir++) {
			for (int k=cellBox.klo(); k<=cellBox.khi()+ZDIM; k++) {
				for (int j=cellBox.jlo(); j<=cellBox.jhi()+YDIM; j++) {
					for (int i=cellBox.ilo(); i<=cellBox.ihi()+XDIM; i++) {
						double pos[] = {
							xlo.x + dx.x * i, 
							xlo.y + dx.y * j,
							xlo.z + dx.z * k,
						};
						fprintf(fp, "%lf ", pos[dir]);
					}
				}
			}
			fprintf(fp, "\n");
		}
	}

	for (int comp=0; comp<ncomp; comp++) {
		for (int kk=validBox.klo(); kk<=validBox.khi(); kk++) {
			for (int jj=validBox.jlo(); jj<=validBox.jhi(); jj++) {
				for (int ii=validBox.ilo(); ii<=validBox.ihi(); ii++) {
					fprintf(fp, "%lf ", blockdata(ii,jj,kk,comp));
				}
			}
		}
		fprintf(fp, "\n");
	}
}



} // namespace_sayaka

