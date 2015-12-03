

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <memory>

#include "log.h"

#include "SayakaCommons.h"
#include "SayakaFillPatch.h"


namespace sayaka
{


void FillPatch::fillBoundary(int scomp, int ncomp, int nlayer)
{
	if (nlayer <= 0) return;

	//const int finest_level = tree.currentFinestLevel();
	const int finest_level = tree.maxRefineLevel();
	assert(finest_level >= tree.minRefineLevel());
	assert(finest_level <= tree.maxRefineLevel());

	for (int level=1; level<=finest_level; level++) {
		int in_level = 0;
		fillLevelBoundary(level, scomp, ncomp, nlayer, in_level);
	}
}

void FillPatch::fillLevelBoundary(int level, 
	int scomp, int ncomp, int nlayer,
	int onlySameLevel, int applyPhysBC, int fillCorner) 
{
	if (nlayer <= 0) return;

	const AmrTreeLevel &tree_level = tree.getTreeLevel(level);
	if (tree_level.isEmptyLevel()) return;

	for (int igrid=0; igrid<tree_level.numLevelBlock; igrid++) {
		// current block
		int iblock = tree_level[igrid];
		assert(tree[iblock].getLevel() == level);

		//
		fillBlockBoundary(iblock, scomp, ncomp, nlayer, 
			onlySameLevel, applyPhysBC, fillCorner);
	}
}

void FillPatch::fillBlockBoundary(int iblock, 
	int scomp, int ncomp, int nlayer,
	int onlySameLevel, int applyPhysBC, int fillCorner)
{
	if (nlayer <= 0) return;

	DiagFillFlag corner_flag = (fillCorner ? DIAG_FILL : DIAG_IGNORE);

	// current block
	const AmrTreeNode &block = tree[iblock];
	const SurroundingBlocks &surr = tree.surrBlocks[iblock];

	// in-level block/block boundary
	fillBlockBndryAtSameLevel(
		iblock, tree_data[iblock], 
		scomp, ncomp, nlayer,
		COORD_ALL, corner_flag);

	// boundary with coarser level
	if (!onlySameLevel && surr.hasCoarseLevelSurrounding()) {
		int iparent = tree[iblock].parent;
		fillBlockBndryFromCoarseLevel(
			iblock, tree_data[iblock],
			iparent, tree_data[iparent],
			scomp, ncomp, nlayer, 
			COORD_ALL, corner_flag);
	}

	// block/BC boundary
	if (applyPhysBC && surr.hasExternalBoundarySurrounding()) {
		fillBlockBndryCond(
			iblock, tree_data[iblock],
			scomp, ncomp, nlayer,
			COORD_ALL, corner_flag);
	}
}


void FillPatch::fillBlockBndryAtSameLevel(
	int iblock, 
	DoubleBlockData &dstdata, // to be boundary-filled
	int scomp, int ncomp, int nlayer, 
	int coordflag, DiagFillFlag cornerflag)
{
	// range of true data only
	const IndexBox &validBox = tree.validBlockCellBox();
	//const int nbx = validBox.size(0);
	//const int nby = validBox.size(1);
	//const int nbz = validBox.size(2);

	// range including ghost data
	//const IndexBox &grownBox = tree.grownBlockCellBox();
	//IndexBox layerBox = validBox;
	//layerBox.extend(nlayer);

	// 
	const bool ignore_corner = (cornerflag == DIAG_IGNORE);

	const SurroundingBlocks &surr = tree.surrBlocks[iblock];

	//
	const FillPatchRange &range_cache = getRangeCache(nlayer);
	assert(range_cache.nlayer == nlayer);

	for (int ksurr=-ZDIM; ksurr<=ZDIM; ksurr++) {
		for (int jsurr=-YDIM; jsurr<=YDIM; jsurr++) {
			for (int isurr=-XDIM; isurr<=XDIM; isurr++) {
				// do not touch the block itself
				if (isurr==0 && jsurr==0 && ksurr==0) continue;
				// skip corners if not needed
				// NOTE that corner is indicated by |ii|+|jj|+|kk|>1
				if (ignore_corner && abs(isurr)+abs(jsurr)+abs(ksurr)>1) continue;

				if (surr.isSameLevelSurrounding(isurr,jsurr,ksurr)) { 
					// in-level surrounding block exists

					/*
					// construct buffer copy index
					// this is the (logical) translation of index
					Vector3i vTrans = Vector3i::VecMake(
						isurr*nbx, jsurr*nby, ksurr*nbz);

					// first destination box for this block
					IndexBox dstbox = validBox;
					dstbox.shift(vTrans);
					dstbox = IndexBox::Intersection(dstbox, layerBox);
					assert(dstbox.isValid());

					// then source box from surrounding block
					IndexBox srcbox = layerBox;
					srcbox.shift(-vTrans);
					srcbox = IndexBox::Intersection(srcbox, validBox);
					assert(srcbox.isValid());

					//
					assert(IndexBox(srcbox).shift(vTrans) == dstbox);
					*/

					// get ranges from cache
					const SurroundingIndex surr_index(isurr,jsurr,ksurr);
					const IndexBox &dstbox = range_cache.dstRanges[surr_index];
					const IndexBox &srcbox = range_cache.srcRanges[surr_index];

					const int jblock = surr(isurr,jsurr,ksurr);
					assert(jblock >= 0);
					const DoubleBlockData &srcdata = tree_data[jblock];

					if (0) {
						if (var_loc == VARLOC_CELL) { 
							// the simplest case
							// cell-centered data can be directly copied
							// by using the target index boxes
							//DoubleBlockData::CopyData(
							DoubleGridDataUtil.CopyData(
								dstdata, srcdata, 
								scomp, scomp, ncomp,
								dstbox, srcbox);
						} else {
							LOGPRINTF("only cell-centered data!\n");
							exit(1);
						}
					} else {
						// in-level boundary is filled by simply copying
						// from the neighbor block
						// NOTE this will cause overlap for face-centered data
						// Is that ok?
						//DoubleBlockData::CopyData(
						DoubleGridDataUtil.CopyData(
							dstdata, srcdata, 
							scomp, scomp, ncomp,
							dstbox, srcbox);
					}
				}
			}
		}
	}
} // fillpatch_fillblockbndryatsamelevel

void FillPatch::fillBlockBndryFromCoarseLevel(
	int iblock, DoubleBlockData &blockdata,
	int iparent, DoubleBlockData &parentdata,
	int scomp, int ncomp, int nlayer, 
	int coordflag, DiagFillFlag cornerflag)
{
	assert(iparent >= 0);
	// needed only if has coarse surrounding
	if (!tree.surrBlocks[iblock].hasCoarseLevelSurrounding()) return;
	
	// cached range of boundary
	const FillPatchRange &range_cache = getRangeCache(nlayer);
	assert(range_cache.nlayer == nlayer);

	const IndexBox &validBox = tree.validBlockCellBox();
	//const int nbx = validBox.size(0);
	//const int nby = validBox.size(1);
	//const int nbz = validBox.size(2);

	//const IndexBox layerBox = IndexBox::Extend(validBox, nlayer);

	//
	const ChildIndex ichild = tree[iblock].whichChild;
	//Vector3i voffset = ChildIndex::ChildIndexOffset(validBox, ichild);
	const Vector3i &voffset = getChildOffsetCache(ichild);
	//const int ioff = voffset.x;
	//const int joff = voffset.y;
	//const int koff = voffset.z;

	// 
	const bool ignore_corner = (cornerflag == DIAG_IGNORE);

	//
	SurroundingBlocks &surr = tree.surrBlocks[iblock];

	for (int ksurr=-ZDIM; ksurr<=ZDIM; ksurr++) {
		for (int jsurr=-YDIM; jsurr<=YDIM; jsurr++) {
			for (int isurr=-XDIM; isurr<=XDIM; isurr++) {
				// skip central block
				if (isurr==0 && jsurr==0 && ksurr==0) continue;
				// skip corner if not needed
				if (ignore_corner && abs(isurr)+abs(jsurr)+abs(ksurr)>1) continue;

				// if block has coarse surrounding in this direction
				// then interpolate boundary values from parent data
				if (surr.isCoarseLevelSurrounding(isurr,jsurr,ksurr)) {
					/*
					Vector3i vshift;
					vshift.x = isurr * nbx;
					vshift.y = jsurr * nby;
					vshift.z = ksurr * nbz;

					// range of boundary data to be filled
					IndexBox bndryBox = validBox;
					bndryBox.shift(vshift);
					bndryBox = IndexBox::Intersection(bndryBox, layerBox);
					assert(bndryBox.isValid());
					*/

					// range of boundary data to be filled
					const SurroundingIndex surr_index(isurr,jsurr,ksurr);
					const IndexBox &bndryBox = range_cache.dstRanges[surr_index];

					if (var_loc.isCellVar()) {
						// cell-centered data
						prolongCrseToFineBox(
							iblock, blockdata,
							iparent, parentdata,
							bndryBox, voffset,
							scomp, ncomp);
					} else if (var_loc.isFaceVar()) {
						// for face-centered data, it is necessary to 
						// adjust the boundary range,
						// so that the fine face data on coarse/fine block boundary
						// is not affected
						/*
						const int stag_dir = var_loc.getFaceVarDir();
						IndexBox targetBox = bndryBox;
						const int ijksurr = SelectDirIndex(stag_dir, isurr,jsurr,ksurr);
						if (ijksurr == -1) { 
							// at low side, shrink high end
							targetBox.extendHigh(stag_dir, -1);
						} else if (ijksurr == 1) {
							// at high side, shrink low end
							targetBox.extendLow(stag_dir, -1);
						}
						*/
						// get non-overlapping boundary range from cache
						const IndexBox &targetBox = range_cache.dstNonOverlapRanges[surr_index];
						
						assert(targetBox.isValid());
						assert(bndryBox.contains(targetBox));

						prolongCrseToFineBox(
							iblock, blockdata,
							iparent, parentdata, 
							targetBox, voffset,
							scomp, ncomp);
					}
					else {
						LOGPRINTF("%s: unsupported data index type=%d!\n", __FUNCTION__, static_cast<int>(var_loc));
						exit(1);
					}
				}
			}
		}
	}
} // fillpatch_fillblockbndryfromcoarselevel



void FillPatch::fillBlockBndryCond(
	int iblock, 
	DoubleBlockData &data, // to be boundary-filled
	int scomp, int ncomp, int nlayer, 
	int coordflag, DiagFillFlag cornerflag)
{
	// range of true data only
	//const IndexBox &validBox = tree.validBlockCellBox();
	//const int nbx = validBox.size(0);
	//const int nby = validBox.size(1);
	//const int nbz = validBox.size(2);

	// range including ghost data
	//const IndexBox &grownBox = tree.grownBlockCellBox();
	//IndexBox layerBox = validBox;
	//layerBox.extend(nlayer);
	
	// 
	const bool ignore_corner = (cornerflag == DIAG_IGNORE);

	const SurroundingBlocks &surr = tree.surrBlocks[iblock];
	if (!surr.hasExternalBoundarySurrounding()) return;

	int basedir = 0;
	if (tree_data.isCellData()) {
		basedir = 0;
	} else {
		LOGPRINTF("%s: cell-centered data only\n", __FUNCTION__);
		exit(1);
	}

	int count = 0;
	for (int sweep=0; sweep<NDIM; sweep++) {
		int sweepdir = (basedir+sweep) % NDIM;

		int bndlo[MAX_DIM] = { 0, 0, 0 };
		int bndhi[MAX_DIM] = { 0, 0, 0 };
		int incr[MAX_DIM] = { 1, 1, 1};

		for (int dir=0; dir<=sweep; dir++) {
			bndlo[(dir+sweep)%NDIM] = -1;
			bndhi[(dir+sweep)%NDIM] = 1;
		}
		incr[sweepdir] = 2;

		for (int ksurr=bndlo[2]; ksurr<=bndhi[2]; ksurr+=incr[2]) {
			for (int jsurr=bndlo[1]; jsurr<=bndhi[1]; jsurr+=incr[1]) {
				for (int isurr=bndlo[0]; isurr<=bndhi[0]; isurr+=incr[0]) {
					count += 1;

					if (ignore_corner && abs(isurr)+abs(jsurr)+abs(ksurr)>1) continue;

					if (surr.isExternalBoundarySurrounding(isurr,jsurr,ksurr)) {
						// hits a physical boundary here
						int phys_bc = surr(isurr,jsurr,ksurr);
						if (bc_patch != NULL) {
							bc_patch->fillBlockBC(phys_bc,
								sweep, sweepdir,
								isurr, jsurr, ksurr, 
								iblock, data, scomp, ncomp);
						} else {
							LOGPRINTF("%s: bc_patch=NULL\n", __FUNCTION__);
							exit(1);
						}

					}
				}
			}
		}
	}
	assert(count == SurroundingIndex::NumSurr-1);
} // fillpatch_fillblockbndrycond



void FillPatch::restrictAll(int scomp, int ncomp)  
{
	// TODO
	// In case that restriction operation need ghost cells
	// we should fill boundary first.
	// As only averaging is implemented, this is not necessary.
	if (0) {
		fillBoundary(scomp, ncomp, tree_data.numGrow());
	}

	// 
	const int finestLevel = tree.currentFinestLevel();

	// NOTE the restriction order must be fine -> coarse
	for (int level=finestLevel-1; level>=1; level--) {
		restrictLevel(level, scomp, ncomp);
	}
}

void FillPatch::restrictLevel(int clevel, int scomp, int ncomp) {
	assert(clevel <= tree.maxRefineLevel());
	if (clevel == tree.maxRefineLevel()) return;

	const AmrTreeLevel &tree_level_crse = tree.getTreeLevel(clevel);

	for (int igrid=0; igrid<tree_level_crse.numLevelBlock; igrid++) {
		//
		const int i = tree_level_crse[igrid];
		const AmrTreeNode &block = tree[i];
		assert(block.getLevel() == clevel);

		if (!block.isLeaf()) { // block covered by children
			// perform restriction taking data from its children
			restrictBlock(i, scomp, ncomp);
		}
	}
}

void FillPatch::restrictBlock(int iblock, int scomp, int ncomp) {
	if (tree[iblock].isLeaf()) return;

	const AmrTreeNode &block = tree[iblock];

	// perform restriction by looping its children
	for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
		//
		const int jblock = block.child[ichild];
		assert(jblock >= 0);

		const IndexBox &validBox = tree.validBlockCellBox();
		//
		const Vector3i fineOffset = ChildIndex::ChildIndexOffset(validBox, ichild);
		IndexBox subCrseBox = ChildIndex::ChildSubIndexRange(validBox, ichild);

		if (var_loc.isFaceVar()) {
			int var_dir = var_loc.getFaceVarDir();
			subCrseBox.staggerInDir(var_dir);
		}

		restrictFineToCrseBox(
			iblock, tree_data[iblock],
			jblock, tree_data[jblock],
			subCrseBox, fineOffset,
			scomp, ncomp);
	} 
}

void FillPatch::prolongLevel(int flevel, int scomp, int ncomp) {
	assert(flevel <= tree.maxAmrLevel);
	if (flevel == 1) return;

	const AmrTreeLevel &tree_level_fine = tree.getTreeLevel(flevel);
	if (tree_level_fine.isEmptyLevel()) return;

	for (int igrid=0; igrid<tree_level_fine.numLevelBlock; igrid++) {
		const int i = tree_level_fine.levelBlock[igrid];
		assert(tree[i].parent >= 0);

		prolongBlock(i, scomp, ncomp);
	}
}

void FillPatch::prolongBlock(int iblock, int scomp, int ncomp) {
	if (!tree[iblock].hasParent()) return;

	//const IndexBox &validBox = tree.validBlockCellBox();

	const ChildIndex ichild = tree[iblock].whichChild;
	//const Vector3i fineOffset = ChildIndex::ChildIndexOffset(validBox, ichild);
	const Vector3i &fineOffset = getChildOffsetCache(ichild);

	const int jparent = tree[iblock].parent;

	const IndexBox &fineBox = tree_data.validBox();

	// interpolate data inside valid range
	prolongCrseToFineBox(
		iblock, tree_data[iblock],
		jparent, tree_data[jparent],
		fineBox, fineOffset,
		scomp, ncomp);
}

void FillPatch::correctFaceBlockData(
	int iface, 
	int iblock, DoubleBlockData &dstdata,
	int jblock, const DoubleBlockData &srcdata,
	int dstcomp, int srccomp, int ncomp) 
{
	const FaceIndex face_block(iface);;
	const FaceIndex face_neigh = face_block.opposite();
	const int face_dir = face_block.dir();
	assert(tree_data.validBox().isStaggeredBox(face_dir));
	assert(dstdata.box().isStaggeredBox(face_dir));
	assert(srcdata.box().isStaggeredBox(face_dir));

	// face-centered valid box
	const IndexBox &validBox = tree_data.validBox();

	// destination box is the face
	const IndexBox dstbox = IndexBox::AdjacentFace(validBox, 
		face_block.dir(), face_block.side());
	// source box is the opposite face
	const IndexBox srcbox = IndexBox::AdjacentFace(validBox, 
		face_neigh.dir(), face_neigh.side());

	//DoubleBlockData::CopyData(
	DoubleGridDataUtil.CopyData(
		dstdata, srcdata,
		dstcomp, srccomp, ncomp,
		dstbox, srcbox);
}

void FillPatch::syncFaceBlockData(int _iface, 
	int iblock, DoubleBlockData &idata,
	int jblock, DoubleBlockData &jdata,
	int scomp, int ncomp)
{
	const FaceIndex iface(_iface);
	const FaceIndex jface = iface.opposite();
	const int face_dir = iface.dir();
	assert(tree_data.validBox().isStaggeredBox(face_dir));
	assert(idata.box().isStaggeredBox(face_dir));
	assert(jdata.box().isStaggeredBox(face_dir));

	// face-centered valid box
	const IndexBox &validBox = tree_data.validBox();

	// destination box is the face
	const IndexBox ibox = IndexBox::AdjacentFace(validBox, 
		iface.dir(), iface.side());
	// source box is the opposite face
	const IndexBox jbox = IndexBox::AdjacentFace(validBox, 
		jface.dir(), jface.side());

	const int &ilo = ibox.ilo();
	const int &ihi = ibox.ihi();
	const int &jlo = ibox.jlo();
	const int &jhi = ibox.jhi();
	const int &klo = ibox.klo();
	const int &khi = ibox.khi();

	const int ioff = jbox.ilo() - ilo;
	const int joff = jbox.jlo() - jlo;
	const int koff = jbox.klo() - klo;

	for (int comp=scomp; comp<scomp+ncomp; comp++) {
		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					// set common face to 0.5*iface + 0.5*jface
					double avg = 0.5 * (idata(i,j,k,comp) + 
						jdata(i+ioff,j+joff,k+koff,comp));
					idata(i,j,k,comp) = avg;
					jdata(i+ioff,j+joff,k+koff,comp) = avg;
				}
			}
		}
	}
}


void FillPatch::setBlockBoundaryValue(
	int iblock, /*DoubleBlockData &out,*/
	int scomp, int ncomp, int nlayer, 
	double bndryValue)
{
	if (nlayer <= 0) return;
	
	assert(0<=scomp && scomp+ncomp<=tree_data.numComp());
	assert(nlayer <= tree_data.numGrow());
	// currently only support cell data
	assert(tree_data.isCellData());

	const FillPatchRange &range_cache = getRangeCache(nlayer);
	assert(range_cache.nlayer == nlayer);

	//
	DoubleBlockData &out = tree_data[iblock];

	for (SurroundingIndex surr_index=0; surr_index<SurroundingIndex::NumSurr; surr_index++) {
		// skip center
		if (surr_index.isCenter()) continue;

		const IndexBox &bndryBox = range_cache.dstRanges[surr_index];

		out.setValue(scomp, ncomp, bndryBox, bndryValue);
	}
}
void FillPatch::setBlockBoundaryValue(
	int iblock, /*DoubleBlockData &out,*/
	int scomp, int ncomp, int nlayer, 
	double bndryValue,
	const SurroundingIndex &surr_index)
{
	if (nlayer <= 0) return;
	
	assert(0<=scomp && scomp+ncomp<=tree_data.numComp());
	assert(nlayer <= tree_data.numGrow());
	// currently only support cell data
	assert(tree_data.isCellData());
	// should not apply on internal part
	assert(!surr_index.isCenter());

	const FillPatchRange &range_cache = getRangeCache(nlayer);
	assert(range_cache.nlayer == nlayer);

	//
	DoubleBlockData &out = tree_data[iblock];

	const IndexBox &bndryBox = range_cache.dstRanges[surr_index];

	out.setValue(scomp, ncomp, bndryBox, bndryValue);
}




void FillPatch::cacheRanges(int nlayer_max) {
	assert(0<=nlayer_max && nlayer_max<=tree_data.numGrow());

	// clear cache
	range_cache.clear();

	// cell box covering the valid range
	// this might be different from the data
	// e.g. the data can be face-centered
	// but in FillPatch we stick to cell-based perspective
	const IndexBox &validBox = tree.validBlockCellBox();
	const Vector3i validSize = validBox.size();

	//
	const IndexBox &dataBox = tree_data.validBox();

	// build cache for all possible numbers of layers
	for (int nlayer=1; nlayer<=nlayer_max; nlayer++) {
		const IndexBox layerBox = IndexBox(validBox).extend(nlayer);

		FillPatchRange range;
		range.nlayer = nlayer;

		// loop surrounding 3x3x3 positions
		for (SurroundingIndex isurr=0; isurr<SurroundingIndex::NumSurr; isurr++) {
			const int ii = isurr.ipos(); assert(-1<=ii && ii<=1);
			const int jj = isurr.jpos(); assert(-1<=jj && jj<=1);
			const int kk = isurr.kpos(); assert(NDIM<3 ? (kk==0) : (-1<=kk && kk<=1));

			Vector3i vshift;
			vshift.x = ii * validSize.x;
			vshift.y = jj * validSize.y;
			vshift.z = kk * validSize.z;

			// boundary ranges of central block
			IndexBox bndryBox = validBox;
			bndryBox.shift(vshift);
			bndryBox = IndexBox::Intersection(bndryBox, layerBox);
			assert(bndryBox.isValid());

			// neighbor ranges covering boundary of central block
			IndexBox neighBox = layerBox;
			neighBox.shift(-vshift); // NOTE the minus sign
			neighBox = IndexBox::Intersection(neighBox, validBox);
			assert(neighBox.isValid());

			// the relationship between DST and SRC boxes
			assert(IndexBox(neighBox).shift(vshift) == bndryBox);

			// save to range
			// cell-based ranges are held for reference
			range.dstCellRanges[isurr] = bndryBox;
			range.srcCellRanges[isurr] = neighBox;
			
			// the exact range, e.g. face-centered
			range.dstRanges[isurr] = bndryBox;
			range.srcRanges[isurr] = neighBox;
			for (int dir=0; dir<NDIM; dir++) {
				// if the exact data is staggered
				// set the target/source to the same
				if (dataBox.type().isStaggered(dir)) {
					range.dstRanges[isurr].staggerInDir(dir);
					range.srcRanges[isurr].staggerInDir(dir);
				}
			}

			// the exact range with no overlapping with valid part
			range.dstNonOverlapRanges[isurr] = range.dstRanges[isurr];
			range.srcNonOverlapRanges[isurr] = range.srcRanges[isurr];
			// adjust the ranges, not adjusted for the center range
			if (!isurr.isCenter()) {
				assert(ii!=0 || jj!=0 || kk!=0);

				for (int dir=0; dir<NDIM; dir++) {
					if (dataBox.isStaggeredBox(dir)) {
						const int ijk = sayaka::SelectDirIndex(dir, ii,jj,kk);

						if (ijk == -1) {
							// at low side, shrink high end
							range.dstNonOverlapRanges[isurr].extendHigh(dir, -1);
							range.srcNonOverlapRanges[isurr].extendHigh(dir, -1);
						} else if (ijk == 1) {
							// at high side, shrink low end
							range.dstNonOverlapRanges[isurr].extendLow(dir, -1);
							range.srcNonOverlapRanges[isurr].extendLow(dir, -1);
						}

						assert(range.dstNonOverlapRanges[isurr].isValid());
						assert(range.srcNonOverlapRanges[isurr].isValid());
						assert(range.dstRanges[isurr].contains(range.dstNonOverlapRanges[isurr]));
						assert(range.srcRanges[isurr].contains(range.srcNonOverlapRanges[isurr]));
					}
				}
			}

		} // end looping surrounding positions

		// put in cache
		range_cache.push_back(range);
	}
}

void FillPatch::cacheOffsets() {
	const IndexBox &validBox = tree.validBlockCellBox();
	for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
		Vector3i fineOffset = ChildIndex::ChildIndexOffset(validBox, ichild);
		child_offset_cache[ichild] = fineOffset;
	}
}

} // namespace_sayaka

