
#include "SayakaMacSolver.h"


namespace sayaka 
{

void MacSolver::setDefaultCoefs() {
	setAlpha(0.0);
	setBeta(1.0);
	setUniformCoefA(0.0);
	setUniformCoefB(1.0);
}

void MacSolver::setUniformCoefA(double a) {
	assert(m_acoef);
	m_acoef->setValue(a);
}
void MacSolver::setUniformCoefB(double b) {
	assert(m_bcoef.size() == NDIM);

	for (int dir=0; dir<NDIM; dir++) {
		assert(m_bcoef[dir]);
		m_bcoef[dir]->setValue(b);
	}
}

void MacSolver::setCoefB(const std::vector<TreeData*> &bdata, int bcomp)
{
	assert(bdata.size() == NDIM);
	for (int dir=0; dir<NDIM; dir++) {
		assert(bdata[dir]);
		setCoefB(dir, *bdata[dir], bcomp);
	}
}
void MacSolver::setCoefB(int dir, const TreeData &bdata, int srccomp)
{
	assert(0<=dir && dir<NDIM);
	assert(bdata.isFaceData(dir));
	assert(0<=srccomp && srccomp<bdata.numComp());

	const AmrTree &tree = getTree();
	int finest_level = tree.currentFinestLevel();
	
	const MGLevelTower &tower = getTower(finest_level);
	assert(!tower.isEmptyLevel());

	assert(m_bcoef.size()==NDIM && m_bcoef[dir]!=NULL);
	TreeData &bcoef = *m_bcoef[dir];
	const int dstcomp = 0;

	// copy on leaf blocks
	TreeData::Copy(bcoef, bdata, dstcomp, srccomp, 1, 0, tower.levelBlock);

	// enforce crse/fine consistency on leaf level
	// since it is the 'coefficient', so average is needed
	// NOTE it assumes coefficients are continuous even at IB surface!!!
	for (int ilevel=finest_level; ilevel>=2; ilevel--) {
		tower.syncSubLevelCrseFineFlux(ilevel, ilevel-1, 
			dir, bcoef, dstcomp, 
			MGLevelTower::FLUX_SYNC_AVG);
	}

	// Now coefficients are ready on the leaf level
	// begin average down to parent levels covered by the finest level
	// NOTE for those coefficients, AVG is used, not SUM

	{
		FillPatch fill(const_cast<AmrTree&>(tree), bcoef);
		fill.setRestriction(dstcomp, InterpPatch::RESTRICT_AVERAGE);
		for (int ilevel=finest_level-1; ilevel>=1; ilevel--) {
			fill.restrictLevel(ilevel, dstcomp, 1);
		}
	}

}


}

