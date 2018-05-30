
#include "SayakaTower.h"
#include "SayakaFillPatch.h"

SAYAKA_NS_BEGIN;


void MGLevelTower::define(const AmrTree &tree, int mg_level) {
	assert(m_tree == NULL);
	// bind to tree
	m_tree = &tree;

	const int maxlevel = tree.maxAmrLevel;
	assert(1<=mg_level && mg_level<=maxlevel);

	super_type::define(mg_level, maxlevel);

	// pre-alloc memory
	const int maxblock = tree.maxBlockNum;

	for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
		this->face_reg[iface].clear();
		this->face_reg[iface].resize(maxblock);
	}
}


void MGLevelTower::update(bool change_to_current_finest_level) {
	const AmrTree &tree = getTree();

	if (change_to_current_finest_level) {
		// change the tower to be the current finest AMR level
		int finest_level = tree.currentFinestLevel();
		assert(tree.minAmrLevel<=finest_level && finest_level<=tree.maxAmrLevel);
		
		setMGLevel(finest_level);
	}

	const int mg_level = getMGLevel();
	assert(1<=mg_level && mg_level<=tree.maxRefineLevel());
	
	// clear current MG level
	clearMGLevel();

	if (tree.getTreeLevel(mg_level).isEmptyLevel()) {
		// if the tree is empty at this level
		// then the tower is also empty
		return;
	} 

	for (int ilevel=mg_level; ilevel>=1; ilevel--) {
		const AmrTreeLevel &level = tree.getTreeLevel(ilevel);

		const IndexBox &validbox = tree.validBlockCellBox();

		// begin index of current level in MG structure
		setBegin(ilevel, numLevelBlock);

		for (int igrid=0; igrid<level.numLevelBlock; igrid++) {
			const int iblock = level.levelBlock[igrid];
			
			registerBlock(mg_level, ilevel, iblock);
		} // end loop blocks on tree level

		// end index of current level in MG structure
		setEnd(ilevel, numLevelBlock);
	} // end loop levels <= mg_level
}


bool MGLevelTower::registerBlock(int mg_level,
	int ilevel, int iblock)
{
	assert(iblock >= 0);
	assert(mg_level == this->getMGLevel());
	assert(1<=ilevel && ilevel<=mg_level);

	const AmrTree &tree = getTree();
	const IndexBox &validbox = tree.validBlockCellBox();

	const AmrTreeNode &block = tree[iblock];
	assert(block.getLevel() == ilevel);

	if (block.isLeaf() || ilevel==mg_level) { 
		// block belongs to the current MG level in two cases: 
		// (a) leaf node with (level<=mg_level)
		// and (b) non-leaf node with (level==mg_level)

		// add to current MG level
		addBlock(iblock);

		// register block faces
		// 1st sweep, decide face types
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			BlockFaceRegister &faceReg = getFaceReg(iblock, face);

			const int ineigh = block.neighbor[face];
			if (ineigh <= NeighborType::PhysBndry) {
				// hit physical boundary
				faceReg.face_type = FACE_FINE_BC;
				faceReg.face_neigh = ineigh;
			} else if (NeighborType::PhysBndry<ineigh && ineigh<0) {
				// should be fine-coarse face
				assert(ineigh == -1);
				assert(ilevel > tree.minAmrLevel);

				// set face type
				faceReg.face_type = FACE_FINE_CRSE;

				// set neighbor to the exact coarse neighbor
				assert(block.parent >= 0);
				int ineigh_crse = tree[block.parent].neighbor[face];
				assert(ineigh_crse >= 0);
				assert(tree[ineigh_crse].isLeaf());
				faceReg.face_neigh = ineigh_crse;
			} else { // ineigh>=0, is valid neighbor
				if (ilevel<mg_level && tree[ineigh].isNotLeaf()) {
					// coarse-fine face
					faceReg.face_type = FACE_CRSE_FINE;
					// still save the coarse neighbor
					faceReg.face_neigh = ineigh;
				} else {
					// fine-fine face
					faceReg.face_type = FACE_FINE_FINE;
					faceReg.face_neigh = ineigh;
				}
			}

			if (mg_level <= tree.minAmrLevel) { 
				// static refined levels
				assert(faceReg.face_type==FACE_FINE_BC ||
					faceReg.face_type==FACE_FINE_FINE);
			}
		} // end loop faces

		
		// 2nd sweep, take care of intra-level faces
		for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
			BlockFaceRegister &faceReg = getFaceReg(iblock, face);

			const int dir = face.dir();
			const int side = face.side();
			const int ii = dir==0 ? 1 : 0;
			const int jj = dir==1 ? 1 : 0;
			const int kk = dir==2 ? 1 : 0;

			const int ineigh = faceReg.face_neigh;

			if (faceReg.face_type == FACE_FINE_CRSE) {
				ChildIndex child_index = block.whichChild;

				bool isSameBlock = true;
				ChildIndex child_index_neigh = child_index.neighbor(face, isSameBlock);
				assert(!isSameBlock);

				faceReg.fine_crse_offset = 
					ChildIndex::ChildIndexOffset(validbox, child_index_neigh);
			} else if (faceReg.face_type == FACE_CRSE_FINE) {
				// loop neighbor's children
				int count = 0;
				for (ChildIndex child=0; child<ChildIndex::NumChild; child++) {
					if (child.hasFaceBoundary(face.opposite())) {
						faceReg.crse_fine_offset[count] = 
							ChildIndex::ChildIndexOffset(validbox, child);
						faceReg.crse_fine_subbox[count] =
							ChildIndex::ChildSubIndexRange(validbox, child);
						faceReg.crse_fine_neigh[count] = 
							tree[ineigh].child[child];

						count += 1;
					}
				}
				assert(count == ChildIndex::NumChildOnFace);
			} else if (faceReg.face_type == FACE_FINE_BC) {
				// TODO

				//LOGPRINTF("%s: BC not implemented\n", __FUNCTION__);
				//exit(1);
			}
		} // end loop faces

		//
		return true;
	} else {
		return false;
	}
}


void MGLevelTower::maskTowerCell(TreeData &mask, 
	int comp, int nlayer) const
{
	// mask must be cell-centered
	assert(mask.isCellData());
	assert(mask.numGrow() >= nlayer);
	assert(mask.numComp() > comp);

	// currently ghost layer must be one
	if (nlayer != 1) {
		LOGPRINTF("%s: nlayer=%d must be 1\n", __FUNCTION__, nlayer);
		exit(1);
	}

	typedef double mask_t;

	const AmrTree &tree = getTree();
	FillPatch fillmask(const_cast<AmrTree&>(tree), mask);

	const IndexBox &validbox = mask.validBox();

	const int mg_level = getMGLevel();
	const MGLevelTower &mg_tower = *this;

	for (int ilevel=mg_level; ilevel>=1; ilevel--) {
		const int begin_level = begin(ilevel);
		const int end_level = end(ilevel);
		
		for (int igrid=begin_level; igrid<end_level; igrid++) {
			const int iblock = mg_tower[igrid];
			DoubleBlockData &maskblock = mask[iblock];

			// first set all mask to invalid
			maskblock.setValue(comp, 1, (mask_t) MASK_INVALID);
			
			// set valid part to internal
			maskblock.setValue(comp, 1, validbox, (mask_t) MASK_INTERNAL);

			// loop faces
			for (FaceIndex face=0; face<FaceIndex::NumFace; face++) {
				const BlockFaceRegister &faceReg = mg_tower.getFaceReg(iblock, face);

				int maskval = MASK_INVALID;
				if (faceReg.isFaceFineFine()) {
					maskval = MASK_INTERNAL;
				} else if (faceReg.isFaceBC()) {
					maskval = MASK_PHYSBC;
				} else if (faceReg.isFaceFineCrse()) {
					maskval = MASK_FINECRSE;
				} else if (faceReg.isFaceCrseFine()) {
					maskval = MASK_CRSEFINE;
				}
				assert(maskval != MASK_INVALID);

				//
				SurroundingIndex surr_face(face);
				assert(surr_face.isFaceSurr(face.dir(), face.side()));

				fillmask.setBlockBoundaryValue(iblock, comp, 1, nlayer, 
					(mask_t) maskval, surr_face);
			}
		}
	}
}

void MGLevelTower::fillTowerBndry(
	TreeData &inout, int scomp, int ncomp, int nlayer,
	int fillPhysBC, int fillCorner, int useInjection) 
{
	if (nlayer <= 0) return;

	assert(inout.isCellData());
	assert(0<=scomp && scomp+ncomp<=inout.numComp());
	assert(nlayer <= inout.numGrow());

	const AmrTree &tree = getTree();
	const MGLevelTower &tower = *this;
	if (tower.isEmptyLevel()) return;

	// the current finest level
	const int mg_level = tower.getMGLevel();

	// temporary filling utility
	FillPatch fill(const_cast<AmrTree&>(tree), inout);
	// set direct injection for fine-coarse boundary fill
	if (useInjection) {
		for (int comp=scomp; comp<scomp+ncomp; comp++) {
			fill.setProlongation(comp, InterpPatch::PROLONG_INJECTION);
		}
	}

	// fill from coarse level
	for (int ilevel=1; ilevel<=mg_level; ilevel++) {
		
		if (ilevel<mg_level && tower.isSubLevelEmpty(ilevel+1)) {
			assert(tower.isSubLevelEmpty(ilevel));
			// skip if next finer level is empty
			continue;
		}

		const int onlySameLevel = 0; // must get data from coarse level
		fill.fillLevelBoundary(ilevel, scomp, ncomp, nlayer,
			onlySameLevel, fillPhysBC, fillCorner);
	}
}



void MGLevelTower::syncSubLevelCrseFineFlux(
	int fine_ilevel, int crse_ilevel,
	int dir, TreeData &fluxdata, int fluxcomp,
	TowerFluxSyncMode sync_mode) const
{
	assert(0<=dir && dir<NDIM);
	assert(0<=fluxcomp && fluxcomp<fluxdata.numComp());
	assert(fluxdata.isFaceData(dir));

	const AmrTree &tree = getTree();
	const MGLevelTower &tower = *this;
	
	assert(2<=fine_ilevel && fine_ilevel<=tower.getMGLevel());
	assert(1<=crse_ilevel && crse_ilevel<=tower.getMGLevel());
	assert(fine_ilevel = crse_ilevel+1);

	// return if no blocks on fine level, i.e. no fine-crse face 
	if (tower.isEmptyLevel()) return;
	if (tower.isSubLevelEmpty(fine_ilevel)) return;
	if (tower.isSubLevelEmpty(crse_ilevel)) return;

	// object used to restrict fine flux -> coarse flux
	// the restriction method is set according to synchronization mode
	InterpPatch flux_interp(const_cast<AmrTree&>(tree), fluxdata, fluxdata.validBox().type());
	if (sync_mode == FLUX_SYNC_AVG) {
		flux_interp.setRestriction(fluxcomp, InterpPatch::RESTRICT_AVERAGE);
	}
	else if (sync_mode == FLUX_SYNC_SUM) {
		flux_interp.setRestriction(fluxcomp, InterpPatch::RESTRICT_SUM);
	}
	else {
		LOGPRINTF("Invalid sync_mode\n"); std::abort();
	}

	//
	const IndexBox &fluxbox = fluxdata.validBox();
	const IndexBox &cellbox = tree.validBlockCellBox();

	//
	FaceIndex face[2];
	IndexBox facebox[2]; 
	for (int s=0; s<=1; s++) {
		face[s] = FaceIndex(dir, s);
		facebox[s] = IndexBox::BoundaryFace(fluxbox, dir, s);
		assert(facebox[s].type() == fluxbox.type() && 
			facebox[s].lo(dir) == fluxbox.endpoint(dir,s) &&
			facebox[s].hi(dir) == fluxbox.endpoint(dir,s));
	}

	//
	Vector3i fineoffset[ChildIndex::NumChild];
	IndexBox crsesubbox[ChildIndex::NumChild][2];
	for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
		//
		fineoffset[ichild] = ChildIndex::ChildIndexOffset(cellbox, ichild);

		for (int side=0; side<=1; side++) {
			const IndexBox childbox = ChildIndex::ChildSubIndexRange(cellbox, ichild);
			crsesubbox[ichild][side] = IndexBox::BoundaryFace(childbox, dir, side);
		}
	}

	// loop fine level
	for (int igrid=tower.begin(fine_ilevel); igrid<tower.end(fine_ilevel); igrid++) {
		const int iblock = tower[igrid];
		assert(iblock >= 0);

		const AmrTreeNode &fine_block = tree[iblock];
		assert(fine_block.getLevel() == fine_ilevel);

		for (int side=0; side<=1; side++) {
			const BlockFaceRegister &face_reg = tower.getFaceReg(iblock, face[side]);

			if (face_reg.isFaceFineCrse()) {
				// on fine-crse face, restrict flux 
				const int iparent = fine_block.parent;
				assert(iparent >= 0);

				const int &child_idx = fine_block.whichChild;

				flux_interp.restrictFineToCrseBox(
					iparent, fluxdata[iparent],
					iblock, fluxdata[iblock],
					crsesubbox[child_idx][side],
					fineoffset[child_idx], 
					fluxcomp, 1);
			}
		}
	} // end loop blocks on fine level

	// loop coarse level
	for (int igrid=tower.begin(crse_ilevel); igrid<tower.end(crse_ilevel); igrid++) {
		const int iblock = tower[igrid];
		assert(iblock >= 0);

		const AmrTreeNode &crse_block = tree[iblock];
		assert(crse_block.getLevel() == crse_ilevel);

		for (int side=0; side<=1; side++) {
			const BlockFaceRegister &face_reg = tower.getFaceReg(iblock, face[side]);

			if (face_reg.isFaceCrseFine()) {
				// on crse-fine face, collect flux from neighbor
				// this is ok because neighbor must be covered by fine level
				// NOTE the face flips to opposite side for the neighbor
				const int ineigh = crse_block.neighbor[face[side]];
				assert(ineigh >= 0);
				assert(tree[ineigh].isNotLeaf());

				DoubleGridDataUtil.CopyData(
					fluxdata[iblock], fluxdata[ineigh],
					fluxcomp, fluxcomp, 1,
					facebox[side], facebox[1-side]);
			}
		}
	} // end loop blocks on coarse level

}


void MGLevelTower::setTowerDataValue(TreeData &data, 
	int scomp, int ncomp, int ngrow, double value) const
{
	assert(0<=scomp && scomp+ncomp<=data.numComp());
	assert(ngrow <= data.numGrow());

	TreeData::SetValue(data, value, 
		scomp, ncomp, ngrow,
		this->levelBlock, 0, this->numLevelBlock);

	//const AmrTree &tree = getTree();
	//const MGLevelTower &tower = *this;
	//if (tower.isEmptyLevel()) return;

	//const IndexBox &validbox = data.validBox();
	//const IndexBox grownbox = IndexBox::Extend(validbox, ngrow);

	//for (int igrid=0; igrid<tower.numLevelBlock; igrid++) {
	//	const int iblock = tower[igrid];

	//	data[iblock].setValue(scomp, ncomp, grownbox, value);
	//}
}
void MGLevelTower::setSubLevelDataValue(
	int ilevel, TreeData &data, 
	int scomp, int ncomp, int ngrow, 
	double value) const
{
	assert(0<=scomp && scomp+ncomp<=data.numComp());
	assert(ngrow <= data.numGrow());

	const AmrTree &tree = getTree();
	const MGLevelTower &tower = *this;
	if (tower.isEmptyLevel()) return;

	assert(1<=ilevel && ilevel<=tower.getMGLevel());
	if (tower.size(ilevel) <= 0) return;

	const int level_begin = tower.begin(ilevel);
	const int level_end = tower.end(ilevel);

	TreeData::SetValue(data, value,
		scomp, ncomp, ngrow, 
		this->levelBlock, level_begin, level_end);

	//const IndexBox &validbox = data.validBox();
	//const IndexBox grownbox = IndexBox::Extend(validbox, ngrow);

	//for (int igrid=level_begin; igrid<level_end; igrid++) {
	//	const int iblock = tower[igrid];
	//	assert(tree[iblock].getLevel() == ilevel);

	//	data[iblock].setValue(scomp, ncomp, grownbox, value);
	//}
}

SAYAKA_NS_END;



