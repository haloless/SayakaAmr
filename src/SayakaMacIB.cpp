
#include "SayakaMacSolver.h"


namespace sayaka 
{

static inline bool IsSolid(const double &ls) {
	return (ls <= 0);
}

void MacSolver::setUnitIBFrac() {
	setUnitIBVolumeFrac(*m_ib_volfrac, 0);
	setUnitIBAreaFrac(m_ib_areafrac, 0);
}
void MacSolver::setUnitIBVolumeFrac(TreeData &volfrac, int volcomp) const {
	assert(volfrac.isCellData());
	assert(0<=volcomp && volcomp<volfrac.numComp());
	TreeData::SetValue(volfrac, 1.0, volcomp, 1, 0);
}
void MacSolver::setUnitIBAreaFrac(std::vector<TreeData*> &areafrac, int areacomp) const {
	assert(areafrac.size() == NDIM);
	for (int dir=0; dir<NDIM; dir++) {
		assert(areafrac[dir]);
		setUnitIBAreaFrac(dir, *areafrac[dir], areacomp);
	}
}
void MacSolver::setUnitIBAreaFrac(int dir, TreeData &areafrac, int areacomp) const {
	assert(areafrac.isFaceData(dir));
	assert(0<=areacomp && areacomp<areafrac.numComp());
	TreeData::SetValue(areafrac, 1.0, areacomp, 1, 0);
}

void MacSolver::setIBFracByLS(const TreeData &ls, int lscomp) {
	calcIBFracByLS(ls, *m_ib_volfrac, m_ib_areafrac, lscomp, 0, 0);
}

void MacSolver::calcIBFracByLS(const TreeData &lsdata, 
	TreeData &volfrac, std::vector<TreeData*> &areafrac,
	int lscomp, int volcomp, int areacomp) const
{
	assert(lsdata.isCellData());
	// ghost cell is needed 
	if (lsdata.numGrow() < 1) {
		LOGPRINTF("%s: LS ngrow must >= 1\n", __FUNCTION__);
		exit(1);
	}

	const AmrTree &tree = getTree();
	int finest_level = tree.currentFinestLevel();
	
	const MGLevelTower &tower = getTower(finest_level);
	assert(!tower.isEmptyLevel());

	// calculate on leaf blocks
	for (int igrid=0; igrid<tower.numLevelBlock; igrid++) {
		int iblock = tower[igrid];
		assert(iblock>=0 && tree[iblock].isLeaf());

		// volume
		calcBlockIBVolumeFrac(iblock, lsdata, volfrac, lscomp, volcomp);
		// area
		for (int dir=0; dir<NDIM; dir++) {
			assert(areafrac[dir]);
			calcBlockIBAreaFrac(iblock, dir, lsdata, *areafrac[dir], lscomp, areacomp);
		}
	}

	// enforce crse/fine consistency on leaf level
	// needed for face area fraction
	// since it is the 'area fraction' (not 'area'), so average is needed
	for (int ilevel=finest_level; ilevel>=2; ilevel--) {
		for (int dir=0; dir<NDIM; dir++) {
			tower.syncSubLevelCrseFineFlux(ilevel, ilevel-1, 
				dir, *areafrac[dir], areacomp, 
				MGLevelTower::FLUX_SYNC_AVG);
		}
	}

	// Now fractions are ready on the leaf level
	// begin average down to parent levels covered by the finest level
	// NOTE for those fractions, AVG is used, not SUM
	{
		FillPatch fill(const_cast<AmrTree&>(tree), volfrac);
		fill.setRestriction(volcomp, InterpPatch::RESTRICT_AVERAGE);
		for (int ilevel=finest_level-1; ilevel>=1; ilevel--) {
			fill.restrictLevel(ilevel, volcomp, 1);
		}
	}
	for (int dir=0; dir<NDIM; dir++) {
		FillPatch fill(const_cast<AmrTree&>(tree), *areafrac[dir]);
		fill.setRestriction(areacomp, InterpPatch::RESTRICT_AVERAGE);
		for (int ilevel=finest_level-1; ilevel>=1; ilevel--) {
			fill.restrictLevel(ilevel, areacomp, 1);
		}
	}
}

void MacSolver::calcBlockIBVolumeFrac(int iblock,
	const TreeData &lsdata, TreeData &volfrac,
	int lscomp, int volcomp) const
{
	assert(lsdata.isCellData());
	assert(volfrac.isCellData());

	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	const DoubleBlockData &ls = lsdata[iblock];
	DoubleBlockData &vfrac = volfrac[iblock];

	for (int k=klo; k<=khi; k++) {
	for (int j=jlo; j<=jhi; j++) {
	for (int i=ilo; i<=ihi; i++) {
		if (IsSolid(ls(i,j,k,lscomp))) {
			vfrac(i,j,k,volcomp) = 0.0;
		} else {
			vfrac(i,j,k,volcomp) = 1.0;
		}
	}
	}
	}
}
void MacSolver::calcBlockIBAreaFrac(int iblock, int dir,
	const TreeData &lsdata, TreeData &areafrac,
	int lscomp, int areacomp) const
{
	assert(lsdata.numGrow() >= 1);
	assert(lsdata.isCellData());
	assert(areafrac.isFaceData(dir));

	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	int ii, jj, kk;
	SelectStaggerIncr(dir, ii,jj,kk);

	const DoubleBlockData &ls = lsdata[iblock];
	DoubleBlockData &sfrac = areafrac[iblock];

	for (int k=klo; k<=khi+kk; k++) {
	for (int j=jlo; j<=jhi+jj; j++) {
	for (int i=ilo; i<=ihi+ii; i++) {
		if (IsSolid(ls(i-ii,j-jj,k-kk,lscomp)) || IsSolid(ls(i,j,k,lscomp))) {
			sfrac(i,j,k,areacomp) = 0.0;
		} else {
			sfrac(i,j,k,areacomp) = 1.0;
		}
	}
	}
	}
}


void MacSolver::mg_restict_wgt(int crse_ilevel, TreeData &phidata, int phicomp)
{
	const AmrTree &tree = getTree();
	if (crse_ilevel == tree.maxAmrLevel) return;
	assert(1<=crse_ilevel && crse_ilevel<tree.maxAmrLevel);

	const AmrTreeLevel &crse_lvl = tree.getTreeLevel(crse_ilevel);
	if (crse_lvl.isEmptyLevel()) return;

	const int fine_ilevel = crse_ilevel + 1;
	const AmrTreeLevel &fine_lvl = tree.getTreeLevel(fine_ilevel);
	if (fine_lvl.isEmptyLevel()) return;

	//
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	//
	TreeData &tmpdata = getDataBuf();
	const int tmpcomp = COMP_BUF_TMP;
	
	//
	const TreeData &ibvfrac = (*m_ib_volfrac);
	const int vfraccomp = 0;

	// copy data on fine level to tmp buffer
	TreeData::Copy(tmpdata, phidata, tmpcomp, phicomp, 1, 0, fine_lvl.levelBlock);

	// on fine level, apply IB volume fraction
	for (int igrid=0; igrid<fine_lvl.numLevelBlock; igrid++) {
		const int iblock = fine_lvl[igrid];

		DoubleBlockData &phi = phidata[iblock];
		const DoubleBlockData &vfrac = ibvfrac[iblock];

		for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
		for (int i=ilo; i<=ihi; i++) {
			phi(i,j,k,phicomp) *= vfrac(i,j,k,vfraccomp);
		}
		}
		}
	}

	// average down
	FillPatch phifill(const_cast<AmrTree&>(tree), phidata);
	phifill.setRestriction(phicomp, InterpPatch::RESTRICT_AVERAGE);
	phifill.restrictLevel(crse_ilevel, phicomp, 1);

	// coarse level, correct by volume fraction
	for (int igrid=0; igrid<crse_lvl.numLevelBlock; igrid++) {
		const int iblock = crse_lvl[igrid];
		if (tree[iblock].isLeaf()) continue; // only non-leaf blocks

		DoubleBlockData &phi = phidata[iblock];
		const DoubleBlockData &vfrac = ibvfrac[iblock];

		for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
		for (int i=ilo; i<=ihi; i++) {
			const double &vf = vfrac(i,j,k,vfraccomp);
			assert(0.0<=vf && vf<=1.0);

			if (vf > 0) {
				phi(i,j,k,phicomp) /= vf;
			} else {
				phi(i,j,k,phicomp) = 0;
			}
		}
		}
		}
	}

	// restore fine data
	TreeData::Copy(phidata, tmpdata, phicomp, tmpcomp, 1, 0, fine_lvl.levelBlock);
}

void MacSolver::mg_prolong_wgt(int fine_ilevel, TreeData &phidata, int phicomp)
{
	const AmrTree &tree = getTree();
	if (fine_ilevel == 1) return;
	assert(1<fine_ilevel && fine_ilevel<=tree.maxAmrLevel);

	const AmrTreeLevel &fine_lvl = tree.getTreeLevel(fine_ilevel);
	if (fine_lvl.isEmptyLevel()) return;

	//const int crse_ilevel = fine_ilevel - 1;
	//const AmrTreeLevel &crse_lvl = tree.getTreeLevel(crse_ilevel);
	//assert(!crse_lvl.isEmptyLevel());

	//
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);
	
	//
	const TreeData &ibvfrac = (*m_ib_volfrac);
	const int vfraccomp = 0;

	// interpolation from coarse
	FillPatch phifill(const_cast<AmrTree&>(tree), phidata);
	phifill.setProlongation(phicomp, InterpPatch::PROLONG_INJECTION);
	phifill.prolongLevel(fine_ilevel, phicomp, 1);

	// fine level, correct by volume fraction
	for (int igrid=0; igrid<fine_lvl.numLevelBlock; igrid++) {
		const int iblock = fine_lvl[igrid];
		assert(tree[iblock].parent >= 0);

		DoubleBlockData &phi = phidata[iblock];
		const DoubleBlockData &vfrac = ibvfrac[iblock];

		for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
		for (int i=ilo; i<=ihi; i++) {
			const double &vf = vfrac(i,j,k,vfraccomp);
			assert(0.0<=vf && vf<=1.0);

			if (vf > 0) {
				// do nothing
			} else {
				phi(i,j,k,phicomp) = 0;
			}
		}
		}
		}
	}
}





} // namespace_sayaka;

