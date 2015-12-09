

#include "SayakaMacSolver.h"


namespace sayaka 
{

typedef MGLevelTower MGLevelRegister;

static const int NumFace = FaceIndex::NumFace;

void MacSolver::define(const AmrTree &tree) {
	const IndexBox &validbox = tree.validBlockCellBox();
	const IndexBox &grownbox = IndexBox::Extend(validbox, NUM_GROW_BUF);

	m_buf = new TreeData(const_cast<AmrTree&>(tree), 
		validbox, grownbox, NUM_GROW_BUF, NUM_COMP_BUF);
	m_buf->setValue(0.0);

	m_fillbuf = new FillPatch(const_cast<AmrTree&>(tree), *m_buf);

	// residual restriction
	m_fillbuf->setProlongation(COMP_BUF_RESID, InterpPatch::PROLONG_INJECTION);
	//m_fillbuf->setRestriction(COMP_BUF_RESID, InterpPatch::RESTRICT_SUM);
	m_fillbuf->setRestriction(COMP_BUF_RESID, InterpPatch::RESTRICT_AVERAGE);

	// correction prolongation
	m_fillbuf->setProlongation(COMP_BUF_CORR, InterpPatch::PROLONG_INJECTION);
	//m_fillbuf->setRestriction(COMP_BUF_CORR, InterpPatch::RESTRICT_SUM);
	m_fillbuf->setRestriction(COMP_BUF_CORR, InterpPatch::RESTRICT_AVERAGE);

	// flux buffer to treat coarse/fine boundary
	m_crse_fine_save.clear();
	m_crse_fine_save = TreeData::CreateFaceDataPArray(tree, 1, 0);
	for (int dir=0; dir<NDIM; dir++) {
		assert(m_crse_fine_save[dir] != NULL);
		m_crse_fine_save[dir]->setValue(0.0);
	}

	// IB fraction
	m_ib_volfrac = TreeData::CreateCellData(tree, 1, 0);
	m_ib_areafrac = TreeData::CreateFaceDataPArray(tree, 1, 0);

	// A & B coefficients
	m_acoef = TreeData::CreateCellData(tree, 1, 0);
	m_bcoef = TreeData::CreateFaceDataPArray(tree, 1, 0);

	//
	m_tower.resize(tree.maxRefineLevel() + 1);
	for (int ilevel=1; ilevel<=tree.maxRefineLevel(); ilevel++) {
		m_tower[ilevel].define(tree, ilevel);
	}

	// BC
	m_bc_rec.resize(tree.maxBlockNum);
	for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		m_bc_rec[iblock].setInteriorBC();
	}
}

void MacSolver::registerMGLevels() {
	const AmrTree &tree = getTree();

	for (int ilevel=tree.maxRefineLevel(); ilevel>=1; ilevel--) {
		registerMGLevel(ilevel);
	}
}

void MacSolver::registerMGLevel(int mg_level) {
	const AmrTree &tree = getTree();
	assert(1<=mg_level && mg_level<=tree.maxRefineLevel());
	
	MGLevelRegister &levelReg = getMGLevelReg(mg_level);
	assert(levelReg.level == mg_level);

	levelReg.update();
}

/*****************************************************************
 * MAC operator
 *****************************************************************/

void MacSolver::mg_apply_level(int mg_level, 
	TreeData &out, TreeData &in, int outcomp, int incomp)
{
	assert(1<=mg_level && mg_level<=getTree().maxRefineLevel());

	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &tower = getTower(mg_level);
	assert(!tower.isEmptyLevel());

	if (useBndryFill()) { // use ghost cell
		// fill boundary first for input data
		const int nlayer = 1;
		mg_fillbndry_level(mg_level, in, incomp, 1, nlayer);
		
		// flux buffer for coarse/fine 
		const int savecomp = 0;
		std::vector<TreeData*> &savedata = getFluxSaveDataBuf();

		// fine -> coarse
		for (int ilevel=mg_level; ilevel>=1; ilevel--) {
			for (int igrid=tower.begin(ilevel); igrid<tower.end(ilevel); igrid++) {
				const int iblock = tower[igrid];
				assert(iblock>=0 && tree[iblock].getLevel()==ilevel);

				mg_apply_block_usefill(
					mg_level, iblock,
					out, in, outcomp, incomp,
					savedata, savecomp);
			}

			// after applying this level, enforce crse/fine block face consistency
			// NOTE applying MAC operator uses area weighting, so SUM is needed, not AVG
			if (ilevel > 1) {
				for (int dir=0; dir<NDIM; dir++) {
					tower.syncSubLevelCrseFineFlux(
						ilevel, ilevel-1, 
						dir, *savedata[dir], savecomp,
						MGLevelTower::FLUX_SYNC_SUM);
				}
			}
		}
	} else { // do not use ghost cell
		for (int igrid=0; igrid<tower.numLevelBlock; igrid++) {
			const int iblock = tower[igrid];
			assert(tree[iblock].getLevel() <= mg_level);

			mg_apply_block_nofill(
				mg_level, iblock,
				out, in, outcomp, incomp);
		}
	}
}


void MacSolver::mg_smooth_level(int mg_level, 
	TreeData &sol, const TreeData &rhs, int solcomp, int rhscomp,
	int redblack_order)
{
	assert(1<=mg_level && mg_level<=getTree().maxRefineLevel());

	if (useBndryFill()) {
		// need flux buffer to save the crse/fine interaction
		std::vector<TreeData*> &cfsave = getFluxSaveDataBuf();
		const int savecomp = 0;
		//
		mg_smooth_level_usefill(mg_level, 
			sol, rhs, solcomp, rhscomp, 
			cfsave, savecomp, 
			redblack_order);
	} else {
		mg_smooth_level_nofill(mg_level, 
			sol, rhs, solcomp, rhscomp);
	}
}

void MacSolver::mg_flux_level(int mg_level,
	int dir, TreeData &flux, TreeData &phi,
	int fluxcomp, int phicomp)
{
	// check flux 
	assert(0<=dir && dir<NDIM);
	assert(flux.isFaceData(dir));

	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &mgLevelReg = getTower(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	// set flux to zero on mg_level
	// this is necessary because flux on coarse/fine face will be accumulated
	mg_zero_level(mg_level, flux, fluxcomp, 1);

	// calculate flux on mg_level
	if (useBndryFill()) { // flux calculation using ghost cell
		// fill boundary
		const int nlayer = 1;
		mg_fillbndry_level(mg_level, phi, phicomp, 1, nlayer);

		// go downwards fine -> coarse
		// as we need fine-coarse flux ready for the coarse level
		for (int ilevel=mg_level; ilevel>=1; ilevel--) {
			// loop blocks on this level
			for (int igrid=mgLevelReg.begin(ilevel); igrid<mgLevelReg.end(ilevel); igrid++) {
				const int iblock = mgLevelReg[igrid];
				assert(tree[iblock].getLevel() == ilevel);

				mg_flux_block_usefill(
					mg_level, iblock, 
					dir, flux, phi, 
					fluxcomp, phicomp);
			}

			// after calculating flux on the current level
			// restrict flux to coarse level and enforce them at coarse/fine interfaces
			// NOTE flux without area weighting needs averaging
			if (ilevel > 1) {
				mgLevelReg.syncSubLevelCrseFineFlux(
					ilevel, ilevel-1,
					dir, flux, fluxcomp, 
					MGLevelTower::FLUX_SYNC_AVG);
			}
		}
	} else { // flux calculation no ghost
		if (1) {
			for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
				const int iblock = mgLevelReg[igrid];
				assert(tree[iblock].getLevel() <= mg_level);

				mg_flux_block_nofill(
					mg_level, iblock, 
					dir, flux, phi, 
					fluxcomp, phicomp);
			}
		} else {
			// go downwards fine -> coarse
			// as we need fine-coarse flux ready for the coarse level
			for (int ilevel=mg_level; ilevel>=1; ilevel--) {
				// range is close/open, i.e. [begin, end)
				const int begin_level = mgLevelReg.begin(ilevel);
				const int end_level = mgLevelReg.end(ilevel);

				for (int igrid=begin_level; igrid<end_level; igrid++) {
					const int iblock = mgLevelReg[igrid];
					assert(tree[iblock].getLevel() == ilevel);

					mg_flux_block_nofill2(mg_level, iblock, 
						dir, flux, phi, 
						fluxcomp, phicomp);
				}
			}
		}
	}
}


void MacSolver::mg_smooth_level_nofill(int mg_level, 
	TreeData &phinew, const TreeData &rhs,
	int outcomp, int rhscomp)
{
	assert(1<=mg_level && mg_level<=getTree().maxRefineLevel());

	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelRegister &mgLevelReg = getMGLevelReg(mg_level);
	assert(!mgLevelReg.isEmptyLevel());

	for (int igrid=0; igrid<mgLevelReg.numLevelBlock; igrid++) {
		const int iblock = mgLevelReg[igrid];
		assert(tree[iblock].getLevel() <= mg_level);

		mg_smooth_block_nofill(
			mg_level, iblock, 
			phinew, rhs,
			outcomp, rhscomp);
	}
}


static void apply_cell_inblock(
	const DoubleBlockData &buf,
	int comp,
	int i, int j, int k, 
	int ilo, int jlo, int klo,
	int ihi, int jhi, int khi,
	const Vector3d &dx, 
	const Vector3d &ds,
	const double dv,
	const double alpha, const double beta,
	double &diag, double &offdiag,
	int onface[FaceIndex::NumFace])
{
	assert(onface);
	// set do not touches face
	std::fill(onface, onface+FaceIndex::NumFace, 0);

	// center contribution
	diag += alpha * dv;

	if (i > ilo) {
		double coef = beta * ds(0) / dx(0);
		offdiag += coef * buf(i-1,j,k,comp);
		diag += coef;
	} else {
		onface[FaceIndex::FACE_XMINUS] = 1;
	} 
	if (i < ihi) {
		double coef = beta * ds(0) / dx(0);
		offdiag += coef * buf(i+1,j,k,comp);
		diag += coef;
	} else {
		onface[FaceIndex::FACE_XPLUS] = 1;
	}
	if (NDIM >= 2) {
		if (j > jlo) {
			double coef = beta * ds(1) / dx(1);
			offdiag += coef * buf(i,j-1,k,comp);
			diag += coef;
		} else {
			onface[FaceIndex::FACE_YMINUS] = 1;
		}
		if (j < jhi) {
			double coef = beta * ds(1) / dx(1);
			offdiag += coef * buf(i,j+1,k,comp);
			diag += coef;
		} else {
			onface[FaceIndex::FACE_YPLUS] = 1;
		}
	}
	if (NDIM == 3) {
		if (k > klo) {
			double coef = beta * ds(2) / dx(2);
			offdiag += coef * buf(i,j,k-1,comp);
			diag += coef;
		} else {
			onface[FaceIndex::FACE_ZMINUS] = 1;
		}
		if (k < khi) {
			double coef = beta *  ds(2) / dx(2);
			offdiag += coef * buf(i,j,k+1,comp);
			diag += coef;
		} else {
			onface[FaceIndex::FACE_ZPLUS] = 1;
		}
	}
}

static void apply_cell_onface(
	const DoubleBlockData &buf, 
	const TreeData &data,
	int comp,
	const FaceIndex &face, 
	const BlockFaceRegister &faceReg,
	int i, int j, int k, 
	int ilo, int jlo, int klo,
	int ihi, int jhi, int khi,
	const Vector3d &dx, 
	const Vector3d &ds,
	const double dv,
	const double alpha, const double beta,
	double &diag, double &offdiag)
{
	const int dir = face.dir();
	const int side = face.side();

	if (faceReg.face_type == FACE_FINE_BC) { // physical BC
		LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
		exit(1);
	} else if (faceReg.face_type == FACE_FINE_FINE) { // in-level fine-fine face
		// mapped index in neighbor block at same level
		int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
		int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
		int kmap = dir!=2 ? k : (side==0 ? khi : klo);

		int neigh = faceReg.face_neigh;
		double val = data[neigh](imap,jmap,kmap,comp);
		double coef = beta * ds(dir) / dx(dir);

		offdiag += coef * val;
		diag += coef;
	} else if (faceReg.face_type == FACE_FINE_CRSE) {
		// mapped index in neighbor block at same level
		int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
		int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
		int kmap = dir!=2 ? k : (side==0 ? khi : klo);

		// the face neighbor holds the coarse neighbor
		// but the mapped index is fine
		// so coarsen the mapped index
		const Vector3i &voff = faceReg.fine_crse_offset;
		int ic = IndexMapping::DirFineToCrse<0>(imap, voff.x);
		int jc = IndexMapping::DirFineToCrse<1>(jmap, voff.y);
		int kc = IndexMapping::DirFineToCrse<2>(kmap, voff.z);

		int neigh_crse = faceReg.face_neigh;
		double val = data[neigh_crse](ic,jc,kc,comp);
		double coef = beta * ds(dir) / (dx(dir)*1.5);

		offdiag += coef * val;
		diag += coef;
	} else if (faceReg.face_type == FACE_CRSE_FINE) {
		// mapped index in neighbor block at same level
		int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
		int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
		int kmap = dir!=2 ? k : (side==0 ? khi : klo);

		// locate mapped coarse index in neighbor's children
		int iloc = -1;
		for (int ind=0; ind<ChildIndex::NumChildOnFace; ind++) {
			if (faceReg.crse_fine_subbox[ind].containsIndex(imap,jmap,kmap)) {
				iloc = ind; break;
			}
		} 
		assert(iloc != -1);

		const Vector3i &crse_fine_offset = faceReg.crse_fine_offset[iloc];
		assert(imap >= crse_fine_offset.x);
		assert(jmap >= crse_fine_offset.y);
		assert(kmap >= crse_fine_offset.z);

		const int &crse_fine_neigh = faceReg.crse_fine_neigh[iloc];
		assert(crse_fine_neigh >= 0);

		// the face neighbor holds the coarse neighbor
		// the mapped index is also coarse
		// we need to refine the mapped index
		for (int kind=0; kind<=ZDIM; kind++) {
			for (int jind=0; jind<=YDIM; jind++) {
				for (int iind=0; iind<=XDIM; iind++) {
					if (SelectDirIndex(dir,iind,jind,kind) == 1-side) {
						int ifine = IndexMapping::DirCrseToFine<0>(imap, iind,
							crse_fine_offset.x);
						int jfine = IndexMapping::DirCrseToFine<1>(jmap, jind,
							crse_fine_offset.y);
						int kfine = IndexMapping::DirCrseToFine<2>(kmap, kind,
							crse_fine_offset.z);

						double val = data[crse_fine_neigh](ifine,jfine,kfine,comp);
						double dsfine = ds(dir) / ChildIndex::NumChildOnFace;
						double coef = beta * dsfine / (dx(dir)*0.75);
						offdiag += coef * val;
						diag += coef;
					}
				}
			}
		}
	} else {
		LOGPRINTF("%s: should never be here\n", __FUNCTION__);
		exit(1);
	}
}

void MacSolver::mg_smooth_block_nofill(
	int mg_level, int iblock, 
	TreeData &phi, const TreeData &rhsold,
	int dstcomp, int rhscomp) 
{
	const AmrTree &tree = getTree();
	
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	const Vector3d ds = tree.getBlockCellArea(iblock);
	const double dv = tree.getBlockCellVolume(iblock);
	// coef
	const double alpha = getAlpha();
	const double beta = getBeta();

	// data buffer
	DoubleBlockData &buf = phi[iblock];
	//const DoubleBlockData &indata = phiold[iblock];
	const DoubleBlockData &rhsdata = rhsold[iblock];

	const int solcomp = dstcomp;

	// MG level
	const MGLevelRegister &level = getMGLevelReg(mg_level);
	assert(!level.isEmptyLevel());

	// total correction
	//double corr = 0;

	for (int rb=0; rb<=1; rb++) { // use GSRB
		for (int k=klo; k<=khi; k++) {
		for (int j=jlo; j<=jhi; j++) {
		for (int i=ilo; i<=ihi; i++) {
			if ((i+j+k)%2 != rb) continue;

			double diag = 0;
			double offdiag = 0;

			int onface[FaceIndex::NumFace] = { 0 };

			// treat in-block contribution
			apply_cell_inblock(
				buf, solcomp,
				i, j, k, 
				ilo, jlo, klo, 
				ihi, jhi, khi,
				dx, ds, dv,
				alpha, beta,
				diag, offdiag,
				onface);

			for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
				if (onface[face]) {
					const BlockFaceRegister &faceReg = level.getFaceReg(iblock, face);

					apply_cell_onface(
						buf, phi, solcomp,
						face, faceReg,
						i, j, k,
						ilo, jlo, klo,
						ihi, jhi, khi,
						dx, ds, dv,
						alpha, beta,
						diag, offdiag);
				}
			}

			// GS relax
			double rhs = rhsdata(i,j,k,rhscomp);
			double uold = buf(i,j,k,solcomp);
			double unew = (rhs*dv + offdiag) / diag;
			buf(i,j,k,solcomp) = unew;

			//corr += abs(unew - uold);
		}
		}
		}
	} // end RB cycle

	//return corr;
}

void MacSolver::mg_apply_block_nofill(
	int mg_level, int iblock,
	TreeData &outdata, const TreeData &indata,
	int dstcomp, int srccomp)
{
	const AmrTree &tree = getTree();
	
	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	const Vector3d ds = tree.getBlockCellArea(iblock);
	const double dv = tree.getBlockCellVolume(iblock);
	// coef
	const double alpha = getAlpha();
	const double beta = getBeta();

	// data buffer
	DoubleBlockData &outbuf = outdata[iblock];
	const DoubleBlockData &inbuf = indata[iblock];

	// MG level
	const MGLevelRegister &level = getMGLevelReg(mg_level);
	assert(!level.isEmptyLevel());

	// range
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	for (int k=klo; k<=khi; k++) {
	for (int j=jlo; j<=jhi; j++) {
	for (int i=ilo; i<=ihi; i++) {
		double diag = 0;
		double offdiag = 0;

		int onface[FaceIndex::NumFace] = { 0 };

		// in-block contribution
		apply_cell_inblock(
			inbuf, srccomp, 
			i, j, k, 
			ilo, jlo, klo,
			ihi, jhi, khi,
			dx, ds, dv,
			alpha, beta,
			diag, offdiag,
			onface);

		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			if (onface[face]) {
				const BlockFaceRegister &faceReg = level.getFaceReg(iblock, face);
				// block face contribution
				apply_cell_onface(
					inbuf, indata, srccomp, 
					face, faceReg,
					i, j, k, 
					ilo, jlo, klo, 
					ihi, jhi, khi,
					dx, ds, dv,
					alpha, beta,
					diag, offdiag);
			}
		}

		// apply operator
		// NOTE the sign of off-diagonal part
		double uin = inbuf(i,j,k,srccomp);
		double val = diag*uin - offdiag;
		// TODO current PPE is in 'conservative form'
		val /= dv;

		outbuf(i,j,k,dstcomp) = val;
	}
	}
	}
}




void MacSolver::mg_flux_block_nofill(
	int mg_level, int iblock,
	int dir, TreeData &fluxdata,
	const TreeData &phidata,
	int fluxcomp, int phicomp)
{
	const AmrTree &tree = getTree();
	const MGLevelRegister &mgLevelReg = getMGLevelReg(mg_level);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	//const Vector3d ds = tree.getBlockCellArea(iblock);
	//const double dv = tree.getBlockCellVolume(iblock);
	const double dh = dx(dir);
	const double dh_crse_fine = dh * 0.75;
	const double dh_fine_crse = dh * 1.5;
	//const double dsdir = ds(dir);

	// staggered index
	//const int ii = dir==0 ? 1 : 0;
	//const int jj = dir==1 ? 1 : 0;
	//const int kk = dir==2 ? 1 : 0;
	int ii, jj, kk;
	SelectStaggerIncr(dir, ii,jj,kk);

	//
	DoubleBlockData &flux = fluxdata[iblock];
	const DoubleBlockData &phi = phidata[iblock];

	// 
	//const IndexBox &fluxbox = fluxdata.validBox();
	//assert(fluxbox.isFaceBox(dir));
	
	//decl_box_range(fluxbox, i,j,k);

	const IndexBox &cellbox = tree.validBlockCellBox();
	decl_box_range(cellbox, i,j,k);

	const int face_end[2] = {
		SelectDirIndex(dir, ilo,jlo,klo),
		SelectDirIndex(dir, ihi+ii,jhi+jj,khi+kk)
	};

	for (int k=klo; k<=khi+kk; k++) {
	for (int j=jlo; j<=jhi+jj; j++) {
	for (int i=ilo; i<=ihi+ii; i++) {
		// the staggered component of index
		const int ijk = SelectDirIndex(dir, i,j,k);

		int side = -1;
		if (ijk == face_end[0]) { // hits low end
			side = 0;
		} else if (ijk == face_end[1]) { // hits high end
			side = 1;
		}

		if (side >= 0) {
			//const double fluxsign = side==0 ? -1.0 : 1.0;

			// current face
			FaceIndex face(dir, side);
			const BlockFaceRegister &faceReg = mgLevelReg.getFaceReg(iblock, face);

			if (faceReg.face_type == FACE_FINE_BC) { // physical BC
				LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
				exit(1);
			} else if (faceReg.face_type == FACE_FINE_FINE) { // in-level fine-fine face
				// (i,j,k) is face index
				// it can also be treated as cell index owing the low face(i,j,k)
				// mapped cell index in neighbor block at same level
				int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
				int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
				int kmap = dir!=2 ? k : (side==0 ? khi : klo);

				int neigh = faceReg.face_neigh;
				double val = phidata[neigh](imap,jmap,kmap,phicomp);
				assert(cellbox.containsIndex(imap,jmap,kmap));

				double phil, phir;
				if (side == 0) {
					phil = val; 
					phir = phi(i,j,k,phicomp);
				} else {
					phil = phi(i-ii,j-jj,k-kk,phicomp);
					phir = val;
				}
				flux(i,j,k,fluxcomp) = -(phir-phil) / dh;
			} else if (faceReg.face_type == FACE_FINE_CRSE) {
				// (i,j,k) is face index
				// it can also be treated as cell index owing the low face(i,j,k)
				// mapped cell index in neighbor block at same level
				int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
				int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
				int kmap = dir!=2 ? k : (side==0 ? khi : klo);

				// the face neighbor holds the coarse neighbor
				// but the mapped index is fine
				// so coarsen the mapped index
				const Vector3i &voff = faceReg.fine_crse_offset;
				int ic = IndexMapping::DirFineToCrse<0>(imap, voff.x);
				int jc = IndexMapping::DirFineToCrse<1>(jmap, voff.y);
				int kc = IndexMapping::DirFineToCrse<2>(kmap, voff.z);

				int neigh_crse = faceReg.face_neigh;
				double val = phidata[neigh_crse](ic,jc,kc,phicomp);

				double phil, phir;
				if (side == 0) {
					phil = val; 
					phir = phi(i,j,k,phicomp);
				} else {
					phil = phi(i-ii,j-jj,k-kk,phicomp);
					phir = val;
				}

				flux(i,j,k,fluxcomp) = -(phir-phil) / dh_fine_crse;
			} else if (faceReg.face_type == FACE_CRSE_FINE) {
				// mapped index in neighbor block at same level
				int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
				int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
				int kmap = dir!=2 ? k : (side==0 ? khi : klo);

				// locate mapped coarse index in neighbor's children
				int iloc = -1;
				for (int ind=0; ind<ChildIndex::NumChildOnFace; ind++) {
					if (faceReg.crse_fine_subbox[ind].containsIndex(imap,jmap,kmap)) {
						iloc = ind; break;
					}
				} 
				assert(iloc != -1);

				const Vector3i &crse_fine_offset = faceReg.crse_fine_offset[iloc];
				assert(imap >= crse_fine_offset.x);
				assert(jmap >= crse_fine_offset.y);
				assert(kmap >= crse_fine_offset.z);

				const int &crse_fine_neigh = faceReg.crse_fine_neigh[iloc];
				assert(crse_fine_neigh >= 0);

				//
				flux(i,j,k,fluxcomp) = 0;

				// the face neighbor holds the coarse neighbor
				// the mapped index is also coarse
				// we need to refine the mapped index
				for (int kind=0; kind<=ZDIM; kind++) {
				for (int jind=0; jind<=YDIM; jind++) {
				for (int iind=0; iind<=XDIM; iind++) {
					if (SelectDirIndex(dir,iind,jind,kind) == 1-side) {
						int ifine = IndexMapping::DirCrseToFine<0>(imap, iind,
							crse_fine_offset.x);
						int jfine = IndexMapping::DirCrseToFine<1>(jmap, jind,
							crse_fine_offset.y);
						int kfine = IndexMapping::DirCrseToFine<2>(kmap, kind,
							crse_fine_offset.z);

						double val = phidata[crse_fine_neigh](ifine,jfine,kfine,phicomp);
				
						double phil, phir;
						if (side == 0) {
							phil = val; 
							phir = phi(i,j,k,phicomp);
						} else {
							phil = phi(i-ii,j-jj,k-kk,phicomp);
							phir = val;
						}
						flux(i,j,k,fluxcomp) += -(phir-phil) / dh_crse_fine;
					}
				}
				}
				}

				flux(i,j,k,fluxcomp) /= ChildIndex::NumChildOnFace;
			} else {
				LOGPRINTF("%s: should never be here\n", __FUNCTION__);
				exit(1);
			}
		} else { // in-block
			// flux = -d(phi)/dx
			flux(i,j,k,fluxcomp) = 
				-(phi(i,j,k,phicomp) - phi(i-ii,j-jj,k-kk,phicomp)) / dh;
		}
	}
	}
	}
}
//
void MacSolver::mg_flux_block_nofill2(
	int mg_level, int iblock,
	int dir, TreeData &fluxdata,
	const TreeData &phidata,
	int fluxcomp, int phicomp)
{
	assert(fluxdata.isFaceData(dir));

	const AmrTree &tree = getTree();
	const MGLevelRegister &mgLevelReg = getMGLevelReg(mg_level);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	//const Vector3d ds = tree.getBlockCellArea(iblock);
	//const double dv = tree.getBlockCellVolume(iblock);
	const double dh = dx(dir);
	const double dh_crse_fine = dh * 0.75;
	const double dh_fine_crse = dh * 1.5;
	//const double dsdir = ds(dir);

	// staggered index
	int ii, jj, kk;
	SelectStaggerIncr(dir, ii,jj,kk);

	//
	DoubleBlockData &flux = fluxdata[iblock];
	const DoubleBlockData &phi = phidata[iblock];

	const IndexBox &cellbox = tree.validBlockCellBox();
	decl_box_range(cellbox, i,j,k);

	const int face_end[2] = {
		SelectDirIndex(dir, ilo,jlo,klo),
		SelectDirIndex(dir, ihi+ii,jhi+jj,khi+kk)
	};

	for (int k=klo; k<=khi+kk; k++) {
	for (int j=jlo; j<=jhi+jj; j++) {
	for (int i=ilo; i<=ihi+ii; i++) {
		// the staggered component of index
		const int ijk = SelectDirIndex(dir, i,j,k);

		int side = -1;
		if (ijk == face_end[0]) { // hits low end
			side = 0;
		} else if (ijk == face_end[1]) { // hits high end
			side = 1;
		}

		if (side >= 0) {
			// current face
			FaceIndex face(dir, side);
			const BlockFaceRegister &faceReg = mgLevelReg.getFaceReg(iblock, face);

			if (faceReg.face_type == FACE_FINE_BC) { // physical BC
				LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
				exit(1);
			} else if (faceReg.face_type == FACE_FINE_FINE) { // in-level fine-fine face
				// (i,j,k) is face index
				// it can also be treated as cell index owing the low face(i,j,k)
				// mapped cell index in neighbor block at same level
				int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
				int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
				int kmap = dir!=2 ? k : (side==0 ? khi : klo);

				int neigh = faceReg.face_neigh;
				double val = phidata[neigh](imap,jmap,kmap,phicomp);
				assert(cellbox.containsIndex(imap,jmap,kmap));

				double phil, phir;
				if (side == 0) {
					phil = val; 
					phir = phi(i,j,k,phicomp);
				} else {
					phil = phi(i-ii,j-jj,k-kk,phicomp);
					phir = val;
				}
				flux(i,j,k,fluxcomp) = -(phir-phil) / dh;
			} else if (faceReg.face_type == FACE_FINE_CRSE) {
				assert(tree[iblock].getLevel() > tree.minAmrLevel);
				// (i,j,k) is face index
				// it can also be treated as cell index owing the low face(i,j,k)
				// mapped cell index in neighbor block at same level
				int imap = dir!=0 ? i : (side==0 ? ihi : ilo);
				int jmap = dir!=1 ? j : (side==0 ? jhi : jlo);
				int kmap = dir!=2 ? k : (side==0 ? khi : klo);

				// the face neighbor holds the coarse neighbor
				// but the mapped index is fine
				// so coarsen the mapped index
				const Vector3i &voff = faceReg.fine_crse_offset;
				int ic = IndexMapping::DirFineToCrse<0>(imap, voff.x);
				int jc = IndexMapping::DirFineToCrse<1>(jmap, voff.y);
				int kc = IndexMapping::DirFineToCrse<2>(kmap, voff.z);

				int neigh_crse = faceReg.face_neigh;
				double val = phidata[neigh_crse](ic,jc,kc,phicomp);

				double phil, phir;
				if (side == 0) {
					phil = val; 
					phir = phi(i,j,k,phicomp);
				} else {
					phil = phi(i-ii,j-jj,k-kk,phicomp);
					phir = val;
				}

				// the fine-coarse flux
				flux(i,j,k,fluxcomp) = -(phir-phil) / dh_fine_crse;

				// accumulate in its coarse neighbor
				const double fine_crse_area_ratio = 1.0 / ChildIndex::NumChildOnFace;
				if (side == 0) {
					fluxdata[neigh_crse](ic+ii,jc+jj,kc+kk,fluxcomp) += 
						fine_crse_area_ratio * flux(i,j,k,fluxcomp);
				} else {
					fluxdata[neigh_crse](ic,jc,kc,fluxcomp) += 
						fine_crse_area_ratio * flux(i,j,k,fluxcomp);
				}
			} else if (faceReg.face_type == FACE_CRSE_FINE) {
				assert(tree[iblock].getLevel()<mg_level && 
					tree[iblock].getLevel()>=tree.minAmrLevel);

				// nothing need to be done here
				// coarse-fine flux is ready as accumulation 
				// of neighbor fine-coarse neighbor fluxes
			} else {
				LOGPRINTF("%s: should never be here\n", __FUNCTION__);
				exit(1);
			}
		} else { // in-block
			// flux = -d(phi)/dx
			flux(i,j,k,fluxcomp) = 
				-(phi(i,j,k,phicomp) - phi(i-ii,j-jj,k-kk,phicomp)) / dh;
		}
	}
	}
	}
}
//
void MacSolver::mg_flux_block_usefill(
	int mg_level, int iblock,
	int dir, TreeData &fluxdata,
	const TreeData &phidata,
	int fluxcomp, int phicomp)
{
	assert(fluxdata.isFaceData(dir));
	// need at least ngrow=1
	assert(phidata.numGrow() >= 1);

	const AmrTree &tree = getTree();
	const MGLevelRegister &mgLevelReg = getMGLevelReg(mg_level);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	//const Vector3d ds = tree.getBlockCellArea(iblock);
	//const double dv = tree.getBlockCellVolume(iblock);
	const double dh = dx(dir);
	const double dh_crse_fine = dh * 0.75;
	const double dh_fine_crse = dh * 1.5;
	//const double dsdir = ds(dir);

	// staggered index
	int ii, jj, kk;
	SelectStaggerIncr(dir, ii,jj,kk);

	//
	DoubleBlockData &flux = fluxdata[iblock];
	const DoubleBlockData &phi = phidata[iblock];

	const IndexBox &cellbox = tree.validBlockCellBox();
	decl_box_range(cellbox, i,j,k);

	const int face_end[2] = {
		SelectDirIndex(dir, ilo,jlo,klo),
		SelectDirIndex(dir, ihi+ii,jhi+jj,khi+kk)
	};

	for (int k=klo; k<=khi+kk; k++) {
	for (int j=jlo; j<=jhi+jj; j++) {
	for (int i=ilo; i<=ihi+ii; i++) {
		// the staggered component of index
		const int ijk = SelectDirIndex(dir, i,j,k);

		int side = -1;
		if (ijk == face_end[0]) { // hits low end
			side = 0;
		} else if (ijk == face_end[1]) { // hits high end
			side = 1;
		}

		if (side >= 0) {
			// current face
			FaceIndex face(dir, side);
			const BlockFaceRegister &faceReg = mgLevelReg.getFaceReg(iblock, face);

			if (faceReg.face_type == FACE_FINE_BC) { // physical BC
				LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
				exit(1);
			} else if (faceReg.face_type == FACE_FINE_FINE) { // in-level fine-fine face
				// flux calculated directly
				flux(i,j,k,fluxcomp) = 
					-(phi(i,j,k,phicomp) - phi(i-ii,j-jj,k-kk,phicomp)) / dh;
			} else if (faceReg.face_type == FACE_FINE_CRSE) {
				assert(tree[iblock].getLevel() > tree.minAmrLevel);

				// coarse value is filled in ghost cells
				// remember to use the fine_crse_dh
				flux(i,j,k,fluxcomp) = 
					-(phi(i,j,k,phicomp) - phi(i-ii,j-jj,k-kk,phicomp)) / dh_fine_crse;
			} else if (faceReg.face_type == FACE_CRSE_FINE) {
				assert(tree[iblock].getLevel()<mg_level && 
					tree[iblock].getLevel()>=tree.minAmrLevel);

				// nothing need to be done here
				// coarse-fine flux is ready as accumulation 
				// of neighbor fine-coarse neighbor fluxes
			} else {
				LOGPRINTF("%s: should never be here\n", __FUNCTION__);
				exit(1);
			}
		} else { // in-block
			// flux = -d(phi)/dx
			flux(i,j,k,fluxcomp) = 
				-(phi(i,j,k,phicomp) - phi(i-ii,j-jj,k-kk,phicomp)) / dh;
		}
	}
	}
	}
}


void MacSolver::mg_smooth_level_usefill(int mg_level,
	TreeData &phi, const TreeData &rhs, int phicomp, int rhscomp,
	std::vector<TreeData*> &fluxsave, int savecomp,
	int redblack_order)
{
	assert(1<=mg_level && mg_level<=getTree().maxRefineLevel());
	assert(fluxsave.size() == NDIM);

	if (phi.numGrow() < 1) {
		LOGPRINTF("%s: require ngrow>=1\n", __FUNCTION__);
		exit(1);
	}

	const AmrTree &tree = getTree();
	if (tree.getTreeLevel(mg_level).isEmptyLevel()) return;

	const MGLevelTower &mg_tower = getTower(mg_level);
	assert(!mg_tower.isEmptyLevel());

	// 
	FillPatch fillphi(const_cast<AmrTree&>(tree), phi);
	fillphi.setProlongation(phicomp, InterpPatch::PROLONG_INJECTION);
	//
	const int ncomp = 1;
	const int nlayer = 1;
	const int nface = FaceIndex::NumFace;

	// fill boundary for whole tower
	// this get all boundary (including crse/fine) ready
	mg_fillbndry_level(mg_level, phi, phicomp, ncomp, nlayer);

	for (int dir=0; dir<NDIM; dir++) {
		fluxsave[dir]->setValue(0.0);
	}

	for (int ilevel=mg_level; ilevel>=1; ilevel--) {
		const int begin_level = mg_tower.begin(ilevel);
		const int end_level = mg_tower.end(ilevel);
		if (begin_level>=end_level) {
			// come down to empty level
			assert(ilevel<mg_level);
			break; 
		}

		if (ilevel > 1) {
			// coarser level exists
			// set crse-fine save data to zero
			//mg_zero_tree_level(mg_level, ilevel-1, crse_fine_save, savecomp, nface);
		}

		// Red-Black
		for (int rb=0; rb<2; rb++) {
			const int rb_phase = (rb+redblack_order) % 2;

			// fill in-level boundary only
			const int onlySameLevel = 1;
			const int fillPhysBC = 0;
			const int fillCorner = 0;
			fillphi.fillLevelBoundary(ilevel, phicomp, ncomp, nlayer,
				onlySameLevel, fillPhysBC, fillCorner);

			// loop blocks on this level
			for (int igrid=begin_level; igrid<end_level; igrid++) {
				const int iblock = mg_tower[igrid];

				mg_smooth_block_rbcycle_usefill(
					mg_level, iblock,
					phi, rhs, phicomp, rhscomp,
					fluxsave, savecomp, 
					rb_phase);
			} // end loop blocks on this level
		} // end Red-Black

		if (ilevel > 1) {
			// coarser level exists
			// enforce fine(level)->coarse(level-1) contribution consistency
			for (int dir=0; dir<NDIM; dir++) {
				mg_tower.syncSubLevelCrseFineFlux(
					ilevel, ilevel-1, 
					dir, *fluxsave[dir], savecomp, 
					MGLevelTower::FLUX_SYNC_SUM);
			}
		}
	}
}


static inline void cell_eval_coef(	
	//const DoubleBlockData &phi, int comp,
	int i, int j, int k, 
	int ilo, int jlo, int klo,
	int ihi, int jhi, int khi,
	const int face_range[MAX_FACE],
	const Vector3d &dx, const Vector3d &ds, const double dv,
	const double alpha, const double beta,
	const BlockFaceRegister* block_face[NumFace],
	const BlockBCRecord &block_bc,
	double &diagcoef, double offcoef[NumFace],
	int bcflag[NumFace],
	int fcflag[NumFace],
	int cfflag[NumFace])
{
	for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
		int dir = face.dir();
		int side = face.side();

		const int ijk = SelectDirIndex(dir, i,j,k);
		if (ijk == face_range[face]) {
			// cell hits face
			const int &face_type = block_face[face]->face_type;

			if (face_type == FACE_FINE_FINE) {
				offcoef[face] = beta * ds(dir) / dx(dir);
			} else if (face_type == FACE_FINE_BC) {
				int bc_type = block_bc.type(face);
				assert(bc_type >= 0);
				if (bc_type == BCType_Neumann) {
					offcoef[face] = 0.0;
					bcflag[face] = 1;
				} else {
					assert(bc_type==BCType_Dirichlet || bc_type==BCType_SimpleFill);
					offcoef[face] = beta * ds(dir) / dx(dir);
					bcflag[face] = 2;
				}
				// BC should have been processed
				//LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
				//exit(1);
			} else if (face_type == FACE_FINE_CRSE) {
				offcoef[face] = beta * ds(dir) / (dx(dir)*1.5);
				fcflag[face] = 1;
			} else if (face_type == FACE_CRSE_FINE) {
				offcoef[face] = beta * ds(dir) / (dx(dir)*0.75);
				cfflag[face] = 1;
			}
		} else {
			// internal cell
			offcoef[face] = beta * ds(dir) / dx(dir);
		}
	}
}

// stencil in each face direction
static const int ioff[MAX_FACE] = {
	-XDIM, +XDIM, 0, 0, 0, 0,
};
static const int joff[MAX_FACE] = {
	0, 0, -YDIM, +YDIM, 0, 0,
};
static const int koff[MAX_FACE] = {
	0, 0, 0, 0, -ZDIM, +ZDIM,
};
// stagger index for each face
static const int istag[MAX_FACE] = {
	0, XDIM, 0, 0, 0, 0,
};
static const int jstag[MAX_FACE] = {
	0, 0, 0, YDIM, 0, 0,
};
static const int kstag[MAX_FACE] = {
	0, 0, 0, 0, 0, ZDIM,
};


void MacSolver::mg_smooth_block_rbcycle_usefill(
	int mg_level, int iblock, 
	TreeData &phidata, const TreeData &rhsdata, int dstcomp, int rhscomp, 
	std::vector<TreeData*> &fluxsave, int savecomp,
	int rb_phase) 
{
	// algorithm requires ghost cells
	assert(phidata.numGrow() >= 1);

	const AmrTree &tree = getTree();
	
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	const Vector3d ds = tree.getBlockCellArea(iblock);
	const double dv = tree.getBlockCellVolume(iblock);
	// coef
	const double alpha = getAlpha();
	const double beta = getBeta();

	// data buffer
	DoubleBlockData &phi = phidata[iblock];
	const DoubleBlockData &rhs = rhsdata[iblock];
	const DoubleBlockData &vfrac = (*m_ib_volfrac)[iblock];

	// MG level
	const MGLevelTower &level = getTower(mg_level);
	assert(!level.isEmptyLevel());
	// block face
	const BlockFaceRegister* blockface[FaceIndex::NumFace];
	level.getBlockFaceRegs(iblock, blockface);
	// block BC
	const BlockBCRecord &blockbc = m_bc_rec[iblock];

	// cell range in each direction
	const int face_range[MAX_FACE] = {
		ilo, ihi, jlo, jhi, klo, khi,
	};

	for (int k=klo; k<=khi; k++) {
	for (int j=jlo; j<=jhi; j++) {
	for (int i=ilo; i<=ihi; i++) {
		if ((i+j+k)%2 != rb_phase) continue;

		double diagcoef = alpha * dv;
		double offcoef[FaceIndex::NumFace] = { 0 };
		// coarse/fine flag
		int fcflag[FaceIndex::NumFace] = { 0 };
		int cfflag[FaceIndex::NumFace] = { 0 };
		// BC flag
		int bcflag[FaceIndex::NumFace] = { 0 };

		cell_eval_coef(
			i, j, k, 
			ilo, jlo, klo, 
			ihi, jhi, khi,
			face_range, 
			dx, ds, dv, 
			alpha, beta, 
			blockface, blockbc, 
			diagcoef, offcoef, 
			bcflag, fcflag, cfflag);

		// TODO move this 
		// correction by IB fraction
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			const int facedir = face.dir();

			const DoubleBlockData &bcoef = (*m_bcoef[facedir])[iblock];
			offcoef[face] *= bcoef(i+istag[face],j+jstag[face],k+kstag[face],0);

			const DoubleBlockData &sfrac = (*m_ib_areafrac[facedir])[iblock];
			offcoef[face] *= sfrac(i+istag[face],j+jstag[face],k+kstag[face],0);
		}
		// currently we assume this...
		assert(diagcoef == 0);

		double offdiag = 0;

		//
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			double offval = 0;
			if (cfflag[face] == 0) {
				if (bcflag[face] == 1) { // Neumann
					offval = 0;
					offcoef[face] = 0;
				} else if (bcflag[face] == 2) { // Dirichlet
					offval = 0;
					offcoef[face] *= 2.0;
				} else {
					assert(bcflag[face] == 0);
					offval = offcoef[face] * phi(i+ioff[face],j+joff[face],k+koff[face],dstcomp);
				} 
			} else {
				//offval = cfsave(i,j,k,savecomp+face);
				const DoubleBlockData &cfflux = (*fluxsave[face.dir()])[iblock];
				offval = cfflux(i+istag[face],j+jstag[face],k+kstag[face],savecomp);
			}

			//
			diagcoef += offcoef[face];
			offdiag += offval;
		}

		// GS relax
		const double rhsval = rhs(i,j,k,rhscomp);
		const double vf = vfrac(i,j,k,0);
		double unew = 0;
		if (vf > 0) {
			assert(diagcoef != 0);
			unew = (rhsval*dv*vf + offdiag) / diagcoef;
		} else {
			assert(diagcoef == 0);
			unew = 0;
		}
		phi(i,j,k,dstcomp) = unew;

		// for fine-crse face, save fine contribution
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			if (fcflag[face]) {
				DoubleBlockData &cfflux = (*fluxsave[face.dir()])[iblock];

				double fc_contrib = offcoef[face] * unew;

				cfflux(i+istag[face],j+jstag[face],k+kstag[face],savecomp) = fc_contrib;
			}
		}
	}
	}
	}
}

void MacSolver::mg_apply_block_usefill(int mg_level, int iblock,
	TreeData &outdata, const TreeData &phidata, int outcomp, int phicomp,
	std::vector<TreeData*> &savedata, int savecomp)
{
	// cell-centered
	assert(outdata.isCellData());
	assert(phidata.isCellData());
	assert(0<=outcomp && outcomp<outdata.numComp());
	assert(0<=phicomp && phicomp<phidata.numComp());
	// algorithm requires ghost cells
	assert(phidata.numGrow() >= 1);
	// 
	assert(savedata.size() == NDIM);

	const AmrTree &tree = getTree();
	
	// index range
	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	const Vector3d ds = tree.getBlockCellArea(iblock);
	const double dv = tree.getBlockCellVolume(iblock);
	// coef
	const double alpha = getAlpha();
	const double beta = getBeta();

	// data buffer
	const DoubleBlockData &phi = phidata[iblock];
	DoubleBlockData &out = outdata[iblock];
	const DoubleBlockData &vfrac = (*m_ib_volfrac)[iblock];

	// MG level
	const MGLevelTower &level = getTower(mg_level);
	assert(!level.isEmptyLevel());
	// block face
	const BlockFaceRegister* blockface[FaceIndex::NumFace];
	level.getBlockFaceRegs(iblock, blockface);
	// block BC
	const BlockBCRecord &blockbc = m_bc_rec[iblock];

	// cell range in each direction
	const int face_range[MAX_FACE] = {
		ilo, ihi, jlo, jhi, klo, khi,
	};

	for (int k=klo; k<=khi; k++) {
	for (int j=jlo; j<=jhi; j++) {
	for (int i=ilo; i<=ihi; i++) {
		// coefficient to be updated
		double diagcoef = alpha * dv;
		double offcoef[NumFace] = { 0 };
		// coarse/fine flag
		int fcflag[NumFace] = { 0 };
		int cfflag[NumFace] = { 0 };
		// BC flag
		int bcflag[NumFace] = { 0 };

		cell_eval_coef(
			i, j, k,
			ilo, jlo, klo,
			ihi, jhi, khi,
			face_range,
			dx, ds, dv,
			alpha, beta,
			blockface, blockbc,
			diagcoef, offcoef,
			bcflag, fcflag, cfflag);

		// TODO move this 
		// correction by IB fraction
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			const int facedir = face.dir();

			const DoubleBlockData &bcoef = (*m_bcoef[facedir])[iblock];
			offcoef[face] *= bcoef(i+istag[face],j+jstag[face],k+kstag[face],0);

			const DoubleBlockData &sfrac = (*m_ib_areafrac[facedir])[iblock];
			offcoef[face] *= sfrac(i+istag[face],j+jstag[face],k+kstag[face],0);
		}
		// currently we assume this...
		assert(diagcoef == 0);

		double offdiag = 0;
		//
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			double offval = 0;
			if (cfflag[face] == 0) {
				if (bcflag[face] == 1) { // Neumann
					offval = 0;
					offcoef[face] = 0;
				} else if (bcflag[face] == 2) { // Dirichlet
					offval = -offcoef[face] * phi(i,j,k,phicomp);
					//offcoef[face] *= 2.0;
				} else {
					assert(bcflag[face] == 0);
					offval = offcoef[face] * phi(i+ioff[face],j+joff[face],k+koff[face],phicomp);
				}
			} else {
				const DoubleBlockData &cfsave = (*savedata[face.dir()])[iblock];
				offval = cfsave(i+istag[face],j+jstag[face],k+kstag[face],savecomp);
			}

			//
			diagcoef += offcoef[face];
			offdiag += offval;
		}

		// apply operator
		const double &phic = phi(i,j,k,phicomp);
		const double Lphi = (phic * diagcoef - offdiag) / dv;

		const double &vf = vfrac(i,j,k,0);
		if (vf > 0) {
			out(i,j,k,outcomp) = Lphi / vf;
		} else {
			assert(diagcoef == 0);
			out(i,j,k,outcomp) = 0;
		}

		// for fine-crse face, save fine contribution
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			if (fcflag[face]) {
				DoubleBlockData &cfsave = (*savedata[face.dir()])[iblock];

				double fc_contrib = offcoef[face] * phic;

				cfsave(i+istag[face],j+jstag[face],k+kstag[face],savecomp) = fc_contrib;
			}
		}
	}
	}
	}
}


void MacSolver::setBC(const BCRegister &bcreg) {
	const AmrTree &tree = getTree();

#pragma omp parallel for
	for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		m_bc_rec[iblock].setBlock(tree[iblock], bcreg);
	}
}
void MacSolver::fixBC(
	const TreeData &phi, TreeData &rhs, 
	int phicomp, int rhscomp) const 
{
	assert(phi.isCellData());
	assert(rhs.isCellData());
	assert(0<=phicomp && phicomp<phi.numComp());
	assert(0<=rhscomp && rhscomp<rhs.numComp());

	// PHI is assumed to hold BC values
	if (phi.numGrow() < 1) {
		LOGPRINTF("%s: BC value NGROW<1\n", __FUNCTION__);
		exit(1);
	}

	const AmrTree &tree = getTree();
	const int finest_level = tree.currentFinestLevel();
	const MGLevelTower &tower = getTower(finest_level);
	assert(!tower.isEmptyLevel());

	// fix RHS on leaf blocks
#pragma omp parallel for
	for (int igrid=0; igrid<tower.numLevelBlock; igrid++) {
		const int iblock = tower[igrid];
	//for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		assert(iblock>=0/* && tree[iblock].isLeaf()*/);

		fixBlockBC(finest_level, iblock,
			phi, rhs, phicomp, rhscomp);
	}
}
void MacSolver::fixBlockBC(int mg_level, int iblock, 
	const TreeData &phidata, TreeData &rhsdata, 
	int phicomp, int rhscomp) const 
{
	const AmrTree &tree = getTree();
	const AmrTreeNode &block = tree[iblock];

	const IndexBox &validbox = tree.validBlockCellBox();
	decl_box_range(validbox, i,j,k);

	// metric
	const Vector3d dx = tree.getBlockCellSize(iblock);
	const Vector3d ds = tree.getBlockCellArea(iblock);
	const double dv = tree.getBlockCellVolume(iblock);
	// coef
	const double alpha = getAlpha();
	const double beta = getBeta();

	//
	const MGLevelTower &tower = getTower(mg_level);
	assert(!tower.isEmptyLevel());
	// block face
	const BlockFaceRegister* block_face[FaceIndex::NumFace];
	tower.getBlockFaceRegs(iblock, block_face);

	// BC
	const BlockBCRecord &block_bcrec = m_bc_rec[iblock];

	const DoubleBlockData &phi = phidata[iblock];
	DoubleBlockData &rhs = rhsdata[iblock];

	const DoubleBlockData &vfrac = (*m_ib_volfrac)[iblock];

	// loop block faces
	for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
		if (block_face[face]->face_type != FACE_FINE_BC) continue;


		int bctype = block_bcrec.type(face);
		//double bcval = block_bcrec.value(face);
		// BCREC must hold math BC >= 0
		assert(bctype >= 0);

		int dir = face.dir();
		int side = face.side();

		const DoubleBlockData &bcoef = (*m_bcoef[dir])[iblock];
		const DoubleBlockData &sfrac = (*m_ib_areafrac[dir])[iblock];

		int imin = dir==0 ? (side==0 ? ilo : ihi) : ilo;
		int imax = dir==0 ? (side==0 ? ilo : ihi) : ihi;
		int jmin = dir==1 ? (side==0 ? jlo : jhi) : jlo;
		int jmax = dir==1 ? (side==0 ? jlo : jhi) : jhi;
		int kmin = dir==2 ? (side==0 ? klo : khi) : klo;
		int kmax = dir==2 ? (side==0 ? klo : khi) : khi;

		for (int k=kmin; k<=kmax; k++) {
		for (int j=jmin; j<=jmax; j++) {
		for (int i=imin; i<=imax; i++) {
			// value in valid cell
			double phi0 = phi(i,j,k,phicomp);
			// value in ghost cell
			double phi1 = phi(i+ioff[face],j+joff[face],k+koff[face],phicomp);

			double coef = beta * ds(dir) / dx(dir);
			coef *= bcoef(i+istag[face],j+jstag[face],k+kstag[face],0);
			coef *= sfrac(i+istag[face],j+jstag[face],k+kstag[face],0);

			double corr = 0;
			if (bctype == BCType_Neumann) {
				// assume homogeneous Neumann
				// do nothing
			} else if (bctype == BCType_Dirichlet) {
				double phibc = 0.5 * (phi0 + phi1);
				corr = phibc * 2.0 * coef;
			} else if (bctype == BCType_SimpleFill) {
				double phibc = phi1;
				corr = phibc * 2.0 * coef;
			} else {
				LOGPRINTF("%s: block=%d face=%d unknown BC=%d\n", __FUNCTION__,
					iblock, (int) face, bctype);
				exit(1);
			}

			if (vfrac(i,j,k,0) > 0) {
				corr /= (dv*vfrac(i,j,k,0));
				rhs(i,j,k,rhscomp) += corr;
			} 
		}
		}
		}
		
	}
}

} // namespace_sayaka


