
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <vector>
#include <algorithm>
#include <memory>

#include "log.h"

#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaBoundaryPatch.h"
#include "SayakaInterpPatch.h"
#include "SayakaFillPatch.h"

namespace sayaka
{



/**
 *
 */
void AmrTree::regrid() {

	// check refine level valid
	if (minAmrLevel<1) {
		LOGPRINTF("%s: minRefineLevel<1\n", __FUNCTION__);
		exit(1);
	}
	if (maxAmrLevel<minAmrLevel) {
		LOGPRINTF("%s: maxRefineLevel<minRefineLevel\n", __FUNCTION__);
		exit(1);
	}

	// 
	std::vector<int> refine_save(maxTreeBlockNum, 0);

	//
	for (int i=0; i<numBlocks; i++) {
		const int level = blocks[i].getLevel();

		if (level >= maxAmrLevel) {
			refineFlag[i] = 0;
		}
		if (level <= minAmrLevel) {
			coarsenFlag[i] = 0;
		}
		if (level<minAmrLevel && blocks[i].isLeaf()) {
			refineFlag[i] = 1;
		}

		refine_save[i] = refineFlag[i];
	}

	fillFlags(newChildFlag, 0);

	// check coarsening on leaf blocks
	checkCoarsen();

	// check refine/coarsen consistency
	// if a leaf is to be coarsened, but its parent is to be refined
	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].isLeaf() && coarsenFlag[i]) {
			int iparent = blocks[i].parent;
			if (refine_save[iparent]) {
				coarsenFlag[i] = 0;
			}
		}
	}

	// check if refinement is to be performed
	int needRegrid = 0;
	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].isLeaf()) {
			if (refineFlag[i] || coarsenFlag[i]) {
				needRegrid = 1;
				break;
			}
		}
	}

	if (needRegrid) {
		int nrefine = 0;
		for (int i=0; i<numBlocks; i++) {
			if (blocks[i].isLeaf() && refineFlag[i]) {
				nrefine += 1;
			}
		}

		// save current 
		int nblock_save = numBlocks;

		//
		if (nrefine > 0) {
			refineBlocks();
		}
		//
		nblock_save = coarsenBlocks(nblock_save);
	}
	LOGPRINTF("%s: regrid=%d, nblock=%d\n", __FUNCTION__, needRegrid, numBlocks);

	// build surrounding blocks list
	{
		int hasCorner = 1;
		buildSurroundingBlocks(hasCorner);
	}

	// cache block indices on each level
	cacheTreeLevels();

	//
	this->gridChanged = 1;


	// Now regrid operation has been finished
	// Generate AMR state data on new blocks
	// NOTE this cannot be put elsewhere above
	// because the data interpolation (prolongation) 
	// depends on a full tree structure


} // amrtree_regrid

void AmrTree::checkCoarsen() {
	//
	std::vector<int> nodetype2(maxTreeBlockNum, 0);
	std::vector<int> refine_tmp(maxTreeBlockNum, 0);
	std::vector<int> refine_tmp2(maxTreeBlockNum, 0);
	std::vector<int> refine_tmp3(maxTreeBlockNum, 0);

	// check that if a block is marked for coarsening 
	// but it is not leaf
	// then do not coarsen 
	for (int i=0; i<numBlocks; i++) {
		if (coarsenFlag[i] && !blocks[i].isLeaf()) {
			coarsenFlag[i] = 0;
		}
	}

	/*
	 * check neighbors
	 */
	const int node_leaf = 1;
	const int node_nonleaf = 2;
	// set node type
	for (int i=0; i<numBlocks; i++) {
		nodetype2[i] = node_leaf;
		// if has children, then type=2
		if (blocks[i].hasChildren()) {
			nodetype2[i] = node_nonleaf;
		}
	}
	// check neighbor is non-leaf
	for (int i=0; i<numBlocks; i++) {
		for (int iface=0; iface<FaceIndex::NumFace; iface++) {
			int ineigh = blocks[i].getNeighbor(iface);
			if (ineigh >= 0) { // has neighbor
				int neigh_type = nodetype2[ineigh];
				if (neigh_type == node_nonleaf) {
					refine_tmp[i] = 1;
					break;
				}
			}
		}
	}
	// check edges
	if (NDIM >= 2) {
		for (int i=0; i<numBlocks; i++) {
			for (int iface=0; iface<FaceIndex::NumFace; iface++) {
				int ineigh = blocks[i].neighbor[iface];
				if (ineigh >= 0) {
					if (refine_tmp[ineigh]) {
						refine_tmp2[i] = 1;
						break;
					}
				}
			}
		}
		for (int i=0; i<numBlocks; i++) {
			if (refine_tmp[i] || refine_tmp2[i]) {
				refine_tmp2[i] = 1;
			}
		}
	}
	// check corners
	if (NDIM == 3) {
		for (int i=0; i<numBlocks; i++) {
			for (int iface=0; iface<FaceIndex::NumFace; iface++) {
				int ineigh = blocks[i].neighbor[iface];
				if (ineigh >= 0) {
					if (refine_tmp2[ineigh]) {
						refine_tmp3[i] = 1;
						break;
					}
				}
			}
		}
	}

	// set coarsen flags based on refine_tmp
	// if leaf node has 
	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].isLeaf()) {
			if (refine_tmp[i] || refine_tmp2[i] || refine_tmp3[i]) {
				coarsenFlag[i] = 0;
			}
		}
	}

	/*
	 * do not coarsen base nodes (no parent)
	 */
	for (int i=0; i<numBlocks; i++) {
		if (blocks[i].parent<0 && coarsenFlag[i]) {
			coarsenFlag[i] = 0;
		}
	}

	/*
	 * check all siblings also marked for coarsening
	 * if not, do not coarsen
	 */
	for (int i=0; i<numBlocks; i++) {
		if (coarsenFlag[i]) {
			int iparent = blocks[i].parent;
			for (int ichild=0; ichild<ChildIndex::NumChild; ichild++) {
				int j = blocks[iparent].child[ichild];
				if (!coarsenFlag[j]) {
					coarsenFlag[i] = 0;
					break;
				}
			}
		}
	}
} // amrtree_checkcoarsen

void AmrTree::refineBlocks() {
	
	std::vector<int> refine_tmp(maxTreeBlockNum, 0);
	std::vector<int> refine_tmp2(maxTreeBlockNum, 0);
	std::vector<int> refine_tmp3(maxTreeBlockNum, 0);


	//
	fillFlags(newChildFlag, 0);

	for (int cycle=0; cycle<99999; cycle++) {
		// final check, don't refine non-leaf block
		for (int i=0; i<numBlocks; i++) {
			if (refineFlag[i] && !blocks[i].isLeaf()) {
				refineFlag[i] = 0;
			}
		}

		int nblock = numBlocks;
		for (int i=0; i<numBlocks; i++) {
			if (refineFlag[i]) { 
				// refine this block
				AmrTreeNode &parentBlock = blocks[i];

				for (TreeChildIndex ichild=0; ichild<TreeChildIndex::NumChild; ichild++) { 
					// create children
					// first check capacity
					if (nblock>= maxTreeBlockNum) {
						LOGPRINTF("Block number overflow\n"); 
						exit(1);
					}

					// sieze child ID
					int j = nblock;
					nblock += 1;

					parentBlock.child[ichild] = j;
					// parent-child relationship
					AmrTreeNode &childBlock = blocks[j];
					childBlock.refineLevel = parentBlock.refineLevel+1;
					childBlock.parent = i;
					childBlock.whichChild = ichild;

					// physical size
					childBlock.boundBox = ChildIndex::ChildBoundBox(parentBlock.boundBox, ichild);
					childBlock.blockCenter = childBlock.boundBox.center();
					childBlock.blockLength = childBlock.boundBox.length();

					// extend flags
					std::copy(parentBlock.blockFlags, parentBlock.blockFlags+MAX_FLAG, childBlock.blockFlags);
					newChildFlag[j] = 1;
				} // end loop children
			} // refine block
		} // loop blocks

		// now connect neighbors for new blocks
		for (int i=0; i<numBlocks; i++) {
			if (refineFlag[i]) {
				// this block has new children
				for (TreeChildIndex ii=0; ii<TreeChildIndex::NumChild; ii++) {

					const int ichild = blocks[i].getChild(ii);

					for (FaceIndex jface=0; jface<FaceIndex::NumFace; jface++) {

						bool isSibling;
						TreeChildIndex jj = ii.neighbor(jface, isSibling);

						if (isSibling) {
							// ii and jj are in the same parent block, easy
							const int jchild = blocks[i].getChild(jj);
							blocks[ichild].neighbor[jface] = jchild;
						} else {
							// ii and jj are in different blocks
							// we have to go through parent's neighbor to locate jchild
							int ineigh = blocks[i].getNeighbor(jface);
							if (ineigh >= 0) {
								// parent has neighbor across this face, then try to find neighbor's child
								int jchild = blocks[ineigh].getChild(jj);
								// save neighbor in this block
								blocks[ichild].neighbor[jface] = jchild;
								// set self in neighbor (if valid)
								if (jchild >= 0) {
									// NOTE that the face changes to opposite side from neighbor's viewpoint
									blocks[jchild].neighbor[jface.opposite()] = ichild;
								}
							} else {
								// parent has no neighbor across this face
								// I think this is an error case 
								// because this means a jump>1 in AMR level
								// possibly on boundary
								//LOGPRINTF("%s: parent=%d child=%d: illegal level jump\n", __FUNCTION__, i, ii);
								//exit(1);
							}
						}

					} // loop face
				} // loop children
			} // newly refined
		} // loop blocks

		// update block number
		this->numBlocks = nblock;

		// set special node type 2
		//for (int i=0; i<maxTreeBlockNum; i++) {
		//	for (int j=0; j<FaceIndex::NumFace; j++) {
		//		blocks[i].nodeType2[j] = 1;
		//	}
		//}
		for (int i=0; i<numBlocks; i++) {
			int *nodeType2 = blocks[i].nodeType2;

			for (int j=0; j<FaceIndex::NumFace; j++) {
				nodeType2[j] = 1;
			}

			// loop its children
			for (ChildIndex ii=0; ii<ChildIndex::NumChild; ii++) {
				int j = blocks[i].child[ii];
				if (j >= 0) {
					// this is a parent(i), tag its faces
					for (int k=0; k<FaceIndex::NumFace; k++) {
						if (nodeType2[k] != 3) {
							nodeType2[k] = 2;
						}
					}

					// children of child(j)
					for (int jj=0; jj<ChildIndex::NumChild; jj++) {
						int k = blocks[j].child[jj];
						if (k >= 0) {
							for (int dir=0; dir<NDIM; dir++) {
								// face in the same direction, but on another side
								FaceIndex ff(dir, 1-ii.side(dir));
								//FaceIndex ff(dir, ii.side(dir));
								nodeType2[ff] = 3;
							}
						}
					}
				}
			} // end loop children of i
		} // end loop nodes

		// reset node types
		setNodeType();

		// check neighbors more than one level difference
		for (int i=0; i<maxTreeBlockNum; i++) {
			refine_tmp[i] = 0;
			refine_tmp2[i] = 0;
			refine_tmp3[i] = 0;
			refineFlag[i] = 0;
		}

		int nref = 0;
		int nref2 = 0;
		int nref3 = 0;

		// 1. check neighbors and set
		for (int i=0; i<numBlocks; i++) {
			for (int iface=0; iface<FaceIndex::NumFace; iface++) {
				int ineigh = blocks[i].neighbor[iface];
				if (ineigh >= 0) {
					if (blocks[ineigh].nodeType2[iface] == 3) {
						refine_tmp[i] = 1;
						nref += 1;
						break;
					}
				}
			}
		}
		if (NDIM >= 2) {
			for (int i=0; i<numBlocks; i++) {
				int ix = 0;
				int iy = 0;
				int iz = 0;

				for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
					int ineigh = blocks[i].neighbor[iface];
					if (ineigh >= 0) {
						if (refine_tmp[ineigh]) {
							if (iface.dir() == 0) {
								ix += 1;
							} else if (iface.dir() == 1) {
								iy += 1;
							} else if (iface.dir() == 2) {
								iz += 1;
							}
						}
					}
				}

				if ((ix>=1 && iy>=1) || (ix>=1 && iz>=1) || (iy>=1 && iz>=1)) {
					refine_tmp2[i] = 1;
					nref2 += 1;
				}
			}
			for (int i=0; i<numBlocks; i++) {
				if (refine_tmp2[i] || refine_tmp[i]) {
					refine_tmp2[i] = 1;
				}
			}
		}
		if (NDIM == 3) {
			for (int i=0; i<numBlocks; i++) {
				int ix = 0;
				int iy = 0;
				int iz = 0;

				for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
					int ineigh = blocks[i].neighbor[iface];
					if (ineigh >= 0) {
						if (refine_tmp2[ineigh]) {
							if (iface.dir() == 0) {
								ix += 1;
							} else if (iface.dir() == 1) {
								iy += 1;
							} else if (iface.dir() == 2) {
								iz += 1;
							}
						}
					}
				}

				if (ix>=1 && iy>=1 && iz>=1) {
					refine_tmp3[i] = 1;
					nref3 += 1;
				}
			}
		}

		// set refine flags
		int repeat = 0;
		for (int i=0; i<numBlocks; i++) {
			refineFlag[i] = 0;
			if (blocks[i].isLeaf()) {
				if (refine_tmp[i] || refine_tmp2[i] || refine_tmp3[i]) {
					repeat = 1;
					refineFlag[i] = 1;
				}
			}
		}

		//LOGPRINTF("%s: cycle=%d\n", __FUNCTION__, cycle);

		if (!repeat) break;
	} 

	// refinement finished, set neighbor for new nodes that are on domain boundary
	for (int i=0; i<numBlocks; i++) {
		if (newChildFlag[i]) { // this is a new child
			// its parent
			int iparent = blocks[i].parent;

			for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
				if (blocks[i].neighbor[iface] < 0) { 
					// may possibly on boundary
					// go up to its parent's neighbor
					int ineigh_parent = blocks[iparent].neighbor[iface];
					if (ineigh_parent <= NeighborType::PhysBndry) {
						blocks[i].neighbor[iface] = ineigh_parent;
					}
				}
			}
		}
	}

	// clear refine flags
	fillFlags(refineFlag, 0);
} // amrtree_refineblocks

/**
 * 
 */
int AmrTree::coarsenBlocks(int oldBlockNum) {
	int nblockOld = oldBlockNum;

	// check
	checkCoarsen();

	//
	std::vector<int> newloc(maxTreeBlockNum, -1);

	for (int i=0; i<numBlocks; i++) {
		if (coarsenFlag[i]) {
			blocks[i].fillBlockFlags(-1);
		}
	}

	// pack blocks 
	int nblock = numBlocks;
	{
		int ipack = 0;

		for (int i=0; i<numBlocks; i++) {
			if (!coarsenFlag[i]) { 
				// save this block
				assert(ipack <= i); 

				// pack to position
				if (ipack != i) {
					blocks[ipack] = blocks[i];
					newChildFlag[ipack] = newChildFlag[i];

					// copy state data associated with the tree
					for (int idata=0; idata<treeStateData.size(); idata++) {
						TreeStateData &state = treeStateData[idata];
						if (0) {
						// old data
						TreeData &sold = state.prevData();
						sold[ipack] = sold[i];
						}
						// new data
						TreeData &snew = state.currData();
						snew[ipack] = snew[i];
					}
				}

				// save packed position
				newloc[i] = ipack;
				ipack += 1;
			} else { 
				// discard this block
				// simply leave it here so it will be overwritten by others
				nblock -= 1;
				nblockOld -= 1;
			}
		}
	}

	// overwrite old locations
	// Since used blocks have been packed,
	// the rest of them between [nb,nb_old) are left unused
	for (int i=nblock; i<numBlocks; i++) {
		coarsenFlag[i] = 0;

		blocks[i].setUnused();

		// set unused data to zero
		for (int idata=0; idata<treeStateData.size(); idata++) {
			TreeStateData &state = treeStateData[idata];
			if (0) {
			// old data
			TreeData &sold = state.prevData();
			sold[i].setValue(0.0);
			}
			// new data
			TreeData &snew = state.currData();
			snew[i].setValue(0.0);
		}
	}

	// update block number
	this->numBlocks = nblock;

	// update pointers to parent/child/neighbor
	for (int i=0; i<numBlocks; i++) {
		//
		if (blocks[i].parent >= 0) {
			int oldpar = blocks[i].parent;
			int newpar = newloc[oldpar];
			assert(newpar >= 0);
			blocks[i].parent = newpar;
		}

		//
		for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
			if (blocks[i].child[ichild] >= 0) {
				int oldchild = blocks[i].child[ichild];
				int newchild = newloc[oldchild];
				blocks[i].child[ichild] = newchild;
			}
		}

		//
		for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
			if (blocks[i].neighbor[iface] >= 0) {
				int oldneigh = blocks[i].neighbor[iface];
				int newneigh = newloc[oldneigh];
				blocks[i].neighbor[iface] = newneigh;
			}
		}
	}
	
	//
	setNodeType();

	//
	fillFlags(coarsenFlag, 0);

	return nblockOld;
} // amrtree_coarsenblocks

void AmrTree::setNodeType() {
	for (int i=0; i<numBlocks; i++) {
		blocks[i].nodeType = TreeNodeType_NoLeafChild;

		for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
			int j = blocks[i].child[ichild];
			if (j < 0) {
				// has no child
				blocks[i].nodeType = TreeNodeType_IsLeaf;
			} else {
				// seek child's child
				if (!blocks[j].hasChildren()) {
					// child(j) has no child
					// so block(i) 
					blocks[i].nodeType = TreeNodeType_HasLeafChild;
				}
			}
		}
	}
} // amrtree_setnodetype

void AmrTree::buildSurroundingBlocks(int buildCorner) {

	for (int i=0; i<numBlocks; i++) {
		SurroundingBlocks &surrblocks = surrBlocks[i];
		surrblocks.setUnused();

		// 1 center: current block
		surrblocks(0,0,0) = i;

		// 6 faces: direct neighbors
		for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
			int j = blocks[i].neighbor[iface];
			int isurr=0, jsurr=0, ksurr=0;
			if (iface == 0) {
				isurr = -1;
				//surrblocks(-1,0,0) = j;
			} else if (iface == 1) {
				isurr = 1;
				//surrblocks(1,0,0) = j;
			} else if (iface == 2) {
				jsurr = -1;
				//surrblocks(0,-1,0) = j;
			} else if (iface == 3) {
				jsurr = 1;
				//surrblocks(0,1,0) = j;
			} else if (iface == 4) {
				ksurr = -1;
				//surrblocks(0,0,-1) = j;
			} else if (iface == 5) {
				ksurr = 1;
				//surrblocks(0,0,1) = j;
			}
			
			surrblocks(isurr,jsurr,ksurr) = j;
			//
			surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(0,0,0);
		}

		// 12 edges: neighbor's neighbor
		if (0) {
		for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
			int j = blocks[i].neighbor[iface];
			int cneigh[MAX_FACE];

			if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
				if (iface == 0) {
					surrblocks(-1,-1,0) = cneigh[2];
					surrblocks(-1,1,0) = cneigh[3];
					surrblocks(-1,0,-1) = cneigh[4];
					surrblocks(-1,0,1) = cneigh[5];
					//
					surrblocks.bcref(-1,-1,0) = SurroundingIndex(-1,0,0);
					surrblocks.bcref(-1,1,0) = SurroundingIndex(-1,0,0);
					surrblocks.bcref(-1,0,-1) = SurroundingIndex(-1,0,0);
					surrblocks.bcref(-1,0,1) = SurroundingIndex(-1,0,0);
				} else if (iface == 1) {
					surrblocks(1,-1,0) = cneigh[2];
					surrblocks(1,1,0) = cneigh[3];
					surrblocks(1,0,-1) = cneigh[4];
					surrblocks(1,0,1) = cneigh[5];
					//
					surrblocks.bcref(1,-1,0) = SurroundingIndex(1,0,0);
					surrblocks.bcref(1,1,0) = SurroundingIndex(1,0,0);
					surrblocks.bcref(1,0,-1) = SurroundingIndex(1,0,0);
					surrblocks.bcref(1,0,1) = SurroundingIndex(1,0,0);
				} else if (iface == 2) {
					surrblocks(-1,-1,0) = cneigh[0];
					surrblocks(1,-1,0) = cneigh[1];
					surrblocks(0,-1,-1) = cneigh[4];
					surrblocks(0,-1,1) = cneigh[5];
					//
					surrblocks.bcref(-1,-1,0) = SurroundingIndex(0,-1,0);
					surrblocks.bcref(1,-1,0) = SurroundingIndex(0,-1,0);
					surrblocks.bcref(0,-1,-1) = SurroundingIndex(0,-1,0);
					surrblocks.bcref(0,-1,1) = SurroundingIndex(0,-1,0);
				} else if (iface == 3) {
					surrblocks(-1,1,0) = cneigh[0];
					surrblocks(1,1,0) = cneigh[1];
					surrblocks(0,1,-1) = cneigh[4];
					surrblocks(0,1,1) = cneigh[5];
					//
					surrblocks.bcref(-1,1,0) = SurroundingIndex(0,1,0);
					surrblocks.bcref(1,1,0) = SurroundingIndex(0,1,0);
					surrblocks.bcref(0,1,-1) = SurroundingIndex(0,1,0);
					surrblocks.bcref(0,1,1) = SurroundingIndex(0,1,0);
				} else if (iface == 4) {
					surrblocks(-1,0,-1) = cneigh[0];
					surrblocks(1,0,-1) = cneigh[1];
					surrblocks(0,-1,-1) = cneigh[2];
					surrblocks(0,1,-1) = cneigh[3];
					//
					surrblocks.bcref(-1,0,-1) = SurroundingIndex(0,0,-1);
					surrblocks.bcref(1,0,-1) = SurroundingIndex(0,0,-1);
					surrblocks.bcref(0,-1,-1) = SurroundingIndex(0,0,-1);
					surrblocks.bcref(0,1,-1) = SurroundingIndex(0,0,-1);
				} else if (iface == 5) {
					surrblocks(-1,0,1) = cneigh[0];
					surrblocks(1,0,1) = cneigh[1];
					surrblocks(0,-1,1) = cneigh[2];
					surrblocks(0,1,1) = cneigh[3];
					//
					surrblocks.bcref(-1,0,1) = SurroundingIndex(0,0,1);
					surrblocks.bcref(1,0,1) = SurroundingIndex(0,0,1);
					surrblocks.bcref(0,-1,1) = SurroundingIndex(0,0,1);
					surrblocks.bcref(0,1,1) = SurroundingIndex(0,0,1);
				}
			}
		}
		} else {
			for (int ksurr=-ZDIM; ksurr<=ZDIM; ksurr++) {
				for (int jsurr=-YDIM; jsurr<=YDIM; jsurr++) {
					for (int isurr=-XDIM; isurr<=XDIM; isurr++) {
						if (abs(isurr)+abs(jsurr)+abs(ksurr) == 2) {
							int f0, f1;
							int i0, j0, k0;
							int i1, j1, k1;
							if (isurr == 0) {
								f0 = jsurr==-1 ? 2 : 3;
								i0 = isurr; j0 = 0; k0 = ksurr;
								f1 = ksurr==-1 ? 4 : 5;
								i1 = isurr; j1 = jsurr; k1 = 0;
							} else if (jsurr == 0) {
								f0 = isurr==-1 ? 0 : 1;
								i0 = 0; j0 = jsurr; k0 = ksurr;
								f1 = ksurr==-1 ? 4 : 5;
								i1 = isurr; j1 = jsurr; k1 = 0;
							} else { // ksurr==0
								f0 = isurr==-1 ? 0 : 1;
								i0 = 0; j0 = jsurr; k0 = ksurr;
								f1 = jsurr==-1 ? 2 : 3;
								i1 = isurr; j1 = 0; k1 = ksurr;
							}

							int lbref = 0;
							if (surrblocks(i0,j0,k0)>=0) { // true cell
								int lb = surrblocks(i0,j0,k0);
								lbref = 1;
								surrblocks(isurr,jsurr,ksurr) = blocks[lb].neighbor[f0];
								surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i0,j0,k0);
							}
							else if (surrblocks(i0, j0, k0) <= NeighborType::PhysBndry) { // BC
								surrblocks(isurr,jsurr,ksurr) = surrblocks(i0,j0,k0);
								surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i0,j0,k0);
							} else { // coarse cell
								lbref = 2;
								//surrblocks(isurr,jsurr,ksurr) = -1;
								surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i0,j0,k0);
							}

							if (surrblocks(i1,j1,k1)>=0) { // true cell
								int lb = surrblocks(i1,j1,k1);
								lbref = 1;
								surrblocks(isurr,jsurr,ksurr) = blocks[lb].neighbor[f1];
								surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i1,j1,k1);
							} else if (surrblocks(i1,j1,k1)<= NeighborType::PhysBndry) { // BC
								// check we do not have true cell record
								if (lbref==0) {
									surrblocks(isurr,jsurr,ksurr) = surrblocks(i1,j1,k1);
									surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i1,j1,k1);
								}
							} else {
								if (lbref == 0) {
									lbref = 2;
									surrblocks.bcref(isurr,jsurr,ksurr) = SurroundingIndex(i1,j1,k1);
								}
							}
						}
					}
				}
			}
		}
		
		// by now we have 9 blocks (2D) or 19 blocks (3D)
		// 8 corners: next treat 3D corners
		if (NDIM==3 && buildCorner) {
			for (int kk=-1; kk<=1; kk+=2) {
				for (int jj=-1; jj<=1; jj+=2) {
					for (int ii=-1; ii<=1; ii+=2) {
						std::pair<int,int> corner = surrblocks.collectCorner(ii,jj,kk, *this);
						surrblocks(ii,jj,kk) = corner.first;
						surrblocks.bcref(ii,jj,kk) = corner.second;
					}
				}
			}


			/*
			int cneigh[MAX_FACE];

			{
				int j = surrblocks(-1,0,-1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(-1,-1,-1) = cneigh[2];
					surrblocks(-1,1,-1) = cneigh[3];
				}
			}{
				int j = surrblocks(-1,0,1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(-1,-1,1) = cneigh[2];
					surrblocks(-1,1,1) = cneigh[3];
				}
			}{
				int j = surrblocks(1,0,-1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(1,-1,-1) = cneigh[2];
					surrblocks(1,1,-1) = cneigh[3];
				}
			}{
				int j = surrblocks(1,0,1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(1,-1,1) = cneigh[2];
					surrblocks(1,1,1) = cneigh[3];
				}
			}

			// fix corner(000)
			if (surrblocks(-1,-1,-1) == -1) { // not found, try another direction
				int j = surrblocks(0,-1,-1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(-1,-1,-1) = cneigh[0];
				}

				if (surrblocks(-1,-1,-1) == -1) { // still not found, try another direction
					int j = surrblocks(-1,-1,0);
					if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
						surrblocks(-1,-1,-1) = cneigh[4];
					}
				}
			}

			// fix corner(010)
			if (surrblocks(-1,1,-1) == -1) {
				int j = surrblocks(0,1,-1);
				if (surrblocks.collectNeighborOrFillBoundary(cneigh, *this, j)) {
					surrblocks(-1,1,-1) = cneigh[0];
				}

				if (surrblocks(-1,1,-1) == -1) {

				}
			}

			// fix corner(001)
			if (surrblocks(-1,-1,1) == -1) {
				
			}

			// fix corner(011)
			if (surrblocks(-1,1,1) == -1) {
				
			}

			// fix corner(100)
			if (surrblocks(1,-1,-1) == -1) {
				
			}

			// fix corner(110)
			if (surrblocks(1,1,-1) == -1) {
				
			}

			// fix corner(101)
			if (surrblocks(1,-1,1) == -1) {
				
			}

			// fix corner(111)
			if (surrblocks(1,1,1) == -1) {
				
			}
			*/
		} // 3D with corner

		// finally collect surrounding node types
		for (int kk=-ZDIM; kk<=ZDIM; kk++) {
			for (int jj=-YDIM; jj<=YDIM; jj++) {
				for (int ii=-XDIM; ii<=XDIM; ii++) {
					if (surrblocks(ii,jj,kk) >= 0) {
						int iblock = surrblocks(ii,jj,kk);
						surrblocks.type(ii,jj,kk) = blocks[iblock].nodeType;
					} else {
						surrblocks.type(ii,jj,kk) = -1;
					}
				}
			}
		}

	} // end loop blocks

} // amrtree_buildsurroundingblocks

void AmrTree::cacheTreeLevel(int ilevel) {
	assert(1<=ilevel && ilevel<=maxAmrLevel);

	AmrTreeLevel &tree_level = getTreeLevel(ilevel);

	// reset block index cached
	tree_level.clear();

	for (int iblock=0; iblock<numBlocks; iblock++) {
		if (blocks[iblock].getLevel() == ilevel) {
			tree_level.numLevelBlock += 1;
			tree_level.levelBlock.push_back(iblock);
		}
	}
}
void AmrTree::cacheTreeLevels() {

	for (int ilevel=1; ilevel<=maxAmrLevel; ilevel++) {
		// reset block index cached
		levels[ilevel].clear();
	}

	for (int iblock=0; iblock<numBlocks; iblock++) {
		int ilevel = blocks[iblock].getLevel();
		assert(1<=ilevel && ilevel<=maxAmrLevel);

		levels[ilevel].numLevelBlock += 1;
		levels[ilevel].levelBlock.push_back(iblock);
	}
}



/**
 * NOTE 
 * Only the new data is processed!
 * Old data is not touched.
 * It is expected in next time step,
 * the old data gets overwritten by the new data!
 */
void AmrTree::initNewBlockDataAfterRegrid() {
	
	//
	int countNew = 0;
	int level_begin = 99;
	int level_end = 0;

	for (int i=0; i<numBlocks; i++) {
		if (newChildFlag[i]) { // find new block
			countNew += 1;

			int levelNew = blocks[i].getLevel();
			level_begin = std::min(level_begin, levelNew);
			level_end = std::max(level_end, levelNew);
		}
	}

	// no new block
	if (countNew <= 0) return;

	assert(1<=level_begin && level_end<=maxAmrLevel);
	assert(level_begin <= level_end);
	assert(level_begin > 1);

	const int ngrow = this->blockNumGrow;

	//
	for (int level=level_begin; level<=level_end; level++) {
		const int level_crse = level - 1;
		assert(level_crse >= 1);

		const AmrTreeLevel &tree_level_crse = getTreeLevel(level_crse);

		// fill coarse level parent block
		for (int igrid=0; igrid<tree_level_crse.numLevelBlock; igrid++) {
			int iblock = tree_level_crse[igrid];
			assert(blocks[iblock].getLevel() == level_crse);
			
			if (!blocks[iblock].isLeaf()) {
				int need_fill = 0;
				
				for (ChildIndex ichild=0; ichild<ChildIndex::NumChild; ichild++) {
					int jblock = blocks[iblock].child[ichild];
					assert(jblock >= 0);
					if (newChildFlag[jblock]) {
						need_fill = 1;
						break;
					}
				}

				if (need_fill) {
					for (int istate=0; istate<treeStateData.size(); istate++) {
						TreeStateData &state = treeStateData[istate];
						TreeData &sdata = state.currData();
						const int scomp = 0;
						const int ncomp = sdata.numComp();
						//const int ngrow = sdata.numGrow();

						FillPatch *fillpatch = state.fill_new_patch;
						assert(fillpatch);

						int inlevel_fill_only = 1;
						fillpatch->fillBlockBoundary(iblock, 
							scomp, ncomp, ngrow, 
							inlevel_fill_only);
					} 
				}
			}
		} // end loop blocks at parent level

		//
		const AmrTreeLevel &tree_level = getTreeLevel(level);

		// begin fill this level
		for (int igrid=0; igrid<tree_level.numLevelBlock; igrid++) {
			int iblock = tree_level[igrid];
			assert(blocks[iblock].getLevel() == level);

			if (newChildFlag[iblock]) {
				for (int istate=0; istate<treeStateData.size(); istate++) {
					TreeStateData &state = treeStateData[istate];
					TreeData &sdata = state.currData();

					const int scomp = 0;
					const int ncomp = sdata.numComp();

					// interpolated from parent level
					FillPatch *fillpatch = state.fill_new_patch;
					assert(fillpatch);

					fillpatch->prolongBlock(iblock, scomp, ncomp);

					// if this is a face-centered data, we make our best to
					// ensure that data on the block face is consistent with 
					// any existing neighbor block
					for (int dir=0; dir<NDIM; dir++) {
						if (sdata.validBox().isStaggeredBox(dir)) {
							// data is staggered in current direction
							// try to fill the block by existing neighbor blocks
							for (int side=0; side<=1; side++) {
								FaceIndex iface(dir, side);
								int jblock = blocks[iblock].neighbor[iface];
								if (jblock>=0 && !newChildFlag[jblock]) {
									assert(blocks[jblock].neighbor[iface.opposite()] == iblock);
									// we find neighbor across this face
									// and the neighbor is not a new block
									fillpatch->correctFaceBlockData(iface,
										iblock, sdata[iblock],
										jblock, sdata[jblock],
										scomp, scomp, ncomp);
								}
							}
						}
					} // end loop directions
				} // end loop states
			} // new block
		} // end loop blocks at this level
	} // end loop levels containing new blocks
}

} // namespace sayaka

