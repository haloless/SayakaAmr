#pragma once


#include <cmath>

#include <iostream>

#include "log.h"
#include "vector3d.h"

#include "SayakaCommons.h"
#include "Sayaka.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"
#include "SayakaBoundaryPatch.h"
#include "SayakaFillPatch.h"
#include "SayakaVarLocInterpolator.h"
#include "SayakaMacSolver.h"


using namespace sayaka;

static void test1() {
	const int nx = 4;
	const int ny = 2;
	const int nz = 1;

	const int ngx = 1;
	const int ngy = 1;
	const int ngz = 0;

	Vector3i domlo = Vector3i::VecMake(0, 0, 0);
	Vector3i domhi = Vector3i::VecMake(nx-1, ny-1, nz-1);
	IndexBox dom(domlo, domhi);

	IndexBox cellbox = dom;
	cellbox.extend(ngx, ngy, ngz);
	LOGPRINTF("cellbox.capacity = %d\n", cellbox.capacity());

	IndexBox xfacebox = cellbox;
	xfacebox.staggerInDir(0);
	LOGPRINTF("xfacebox.capacity = %d\n", xfacebox.capacity());

	IndexBox nodebox = cellbox;
	nodebox.stagger(1, 1, 0);
	LOGPRINTF("nodebox.capacity = %d\n", nodebox.capacity());

	Vector3d poslo = Vector3d::VecMake(0.0, 0.0, 0.0);
	Vector3d poshi = Vector3d::VecMake(1.0, 1.0, 1.0);
	RealBox bound(poslo, poshi);

	const int ncomp = 2;
	GridData<double> g(cellbox, ncomp);
	g.setValue(0);

	{
		GridData<double> g2 = g;
		g2.setValue(1);
	}

	
	for (ChildIndex i=0; i<ChildIndex::NumChild; i++) {
		//LOGPRINTF("Child%d (%d,%d,%d)\n", static_cast<int>(i), i.iside(), i.jside(), i.kside());
		LOGPRINTF("Child%d (%d,%d,%d)\n", static_cast<int>(i), i.side(0), i.side(1), i.side(2));
		RealBox childBox = ChildIndex::ChildBoundBox(bound, i);
		LOGPRINTF("Bound (%f,%f,%f) (%f,%f,%f)\n", 
			childBox.lo().x, childBox.lo().y, childBox.lo().z, 
			childBox.hi().x, childBox.hi().y, childBox.hi().z);

		for (int dir=0; dir<NDIM; dir++) {
			for (int side=0; side<2; side++) {
				bool isSameBlock = true;
				ChildIndex ineigh = i.neighbor(dir, side, isSameBlock);
				LOGPRINTF("Neigh(dir=%d,side=%d) %d inblock=%d\n", 
					dir, side, ineigh, (int)isSameBlock);
			}
		}
	}
	for (FaceIndex i=0; i<FaceIndex::NumFace; i++) {
		LOGPRINTF("Face%d (dir=%d,side=%d)\n", i, i.dir(), i.side());
	}
}


struct ProbParam
{
	Vector3d circ_cen0;
	Vector3d circ_cen;
	double circ_rad;
};

static void test_refinement(AmrTree &tree, const ProbParam &param) {
	for (int i=0; i<tree.numBlocks; i++) {
		AmrTreeNode &block = tree.blocks[i];
		int level = block.getLevel();

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			Vector3d cen = tree.blocks[i].blockCenter;
			double rcen = Vec3Length(cen, param.circ_cen);
			double dist = rcen - param.circ_rad;

			double dh = block.blockLength.x;
			int needRefine = (abs(dist) <= dh*2.0);
			int needCoarsen = (abs(dist) >= dh*4.0);

			if (level < tree.maxAmrLevel) {
				if (needRefine) {
					tree.refineFlag[i] = 1;
				}
			}

			if (level>tree.minAmrLevel && !tree.refineFlag[i]) {
				if (needCoarsen) {
					tree.coarsenFlag[i] = 1;
				}
			}
		}
	}
}

void test2() {
	const int minLevel = 3;
	const int maxLevel = 8;

	const double xmin = -4.0;
	const double xmax = 4.0;
	const double ymin = -4.0;
	const double ymax = 4.0;
	const double zmin = -1.0;
	const double zmax = 1.0;
	Vector3d probmin = Vector3d::VecMake(xmin, ymin, zmin);
	Vector3d probmax = Vector3d::VecMake(xmax, ymax, zmax);

	//
	ProbParam param;
	param.circ_cen0 = Vector3d::VecMake(xmin, ymin, 0.0);
	param.circ_cen = param.circ_cen0;
	param.circ_rad = 2.0;

	
	AmrTree tree;
	//
	tree.minAmrLevel = minLevel;
	tree.maxAmrLevel = maxLevel;
	// allocate tree memory
	tree.init();

	
	{
		// setup base level
		// a single node on the base level
		tree.numBlocks = 1;
		AmrTreeNode &root = tree.blocks[0];

		RealBox bbox(probmin, probmax);
		root.boundBox = bbox;
		root.blockCenter = bbox.center();
		root.blockLength = bbox.length();

		root.nodeType = TreeNodeType_IsLeaf;
		root.setLevel(1);

		// periodic by referring to itself
		for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
			root.neighbor[iface] = 0;
		}
	}
	if (1) {
		// static refinement up to min level
		for (int level=1; level<minLevel; level++) {
			for (int i=0; i<tree.numBlocks; i++) {
				tree.refineFlag[i] = 1;
			}
			tree.regrid();
		}
	}
	if (1) {
		// dynamic refinement up to max level
		for (int level=minLevel; level<maxLevel; level++) {
			test_refinement(tree, param);
			tree.regrid();
		}
		{
			char filename[64];
			sprintf(filename, "output/amr_leaf_%04d.tec.dat", 0);
			tree.writeTreeLeafState(filename, 1, 0);
		}
	}

	const int nstep = 100;
	for (int step=1; step<=nstep; step++) {
		Vector3d circ_sft = Vector3d::VecMake((xmax-xmin)/nstep*step, (ymax-ymin)/nstep*step, 0.0);
		param.circ_cen = param.circ_cen0 + circ_sft;

		//
		//const int ncycle = maxLevel-minLevel+1;
		const int ncycle = 1;
		//for (int cycle=0; cycle<ncycle; cycle++) {
		//	for (int i=0; i<tree.numBlocks; i++) {
		//		//tree.refineFlag[i] = 1;
		//		Vector3d cen = tree.blocks[i].blockCenter;
		//		double rcen = Vec3Length(&cen, &circ_cen);
		//		double dist = rcen - circ_rad;
		//		if (abs(dist) <= tree.blocks[i].blockLength.x*2.0) {
		//			tree.refineFlag[i] = 1;
		//		} else if (abs(dist) >= tree.blocks[i].blockLength.x*4.0) {
		//			tree.coarsenFlag[i] = 1;
		//		}
		//	}
		//	tree.regrid();
		//}

		for (int cycle=0; cycle<ncycle; cycle++) {
			test_refinement(tree, param);
			tree.regrid();
		}
		

		{
			char filename[64];
			sprintf(filename, "output/amr_leaf_%04d.tec.dat", step);
			tree.writeTreeLeafState(filename, step+1, step);
		}

	}
}


#define TESTDIR(path) ("test2" path)
class Test2
{
	enum {
		PhysBC_Fixed = -100,
		PhysBC_Symmetry,
	};

	enum {
		Comp_Phi = 0,
		NumComp,
	};

	AmrTree *m_tree;
	TreeData *m_data;
	BCPatch *m_bcpatch;
	FillPatch *m_fillpatch;
	TreeStateData *m_state;

	TreeData *m_node_data;

	int m_step;
	double m_time;

	Vector3d probmin, probmax;
	RealBox probbox;

	int n0x, n0y, n0z; // tile at root level

	int nb;
	int nbx, nby, nbz;
	int ng;
	int ngx, ngy, ngz;
	IndexBox validBox, grownBox;

	int minlevel, maxlevel;

	// 
	Vector3d circ_cen;
	double circ_rad;

public:

	void init(int argc, char *argv[]) {


		m_step = 0;
		m_time = 0;

		probmin = Vector3d::VecMake(-4.0, -4.0, 0.0);
		probmax = Vector3d::VecMake(4.0, 4.0, 0.0);
		probbox = RealBox(probmin, probmax);

		// problem
		circ_cen = Vector3d::VecMake(0.0, 0.0, 0.0);
		circ_rad = 2.0;

		//
		n0x = 2;
		n0y = 2;
		n0z = 1;

		minlevel = 3;
		maxlevel = 6;

		nb = 4;
		nbx = nb;
		nby = nb;
		nbz = 1;
		ng = 2;
		ngx = ng;
		ngy = ng;
		ngz = 0;

		validBox = IndexBox(
			Vector3i::VecMake(0,0,0),
			Vector3i::VecMake(nbx-1,nby-1,nbz-1));
		grownBox = IndexBox(
			IndexBox(validBox).extend(ngx,ngy,ngz));


		//
		this->init_tree();
		if (0) {
			// special refinement for coarsening test
			for (int level=minlevel; level<std::min(maxlevel,minlevel+1); level++) {
				for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
					m_tree->refineFlag[iblock] = 1;
				}
				m_tree->regrid();
			}
		}

		//
		VarLocInterpolator varloc_interp(*m_tree);

		//
		this->init_data();
		//
		this->init_patch();
		// 
		this->init_state();
		////
		//m_fillpatch->fillBoundary(0, NumComp, ng);

		m_fillpatch->fillBoundary(0, NumComp, 1);
		// cell->node
		varloc_interp.interp_cell_to_node(*m_node_data, *m_data, Comp_Phi, Comp_Phi, NumComp, 0);

		// save data
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid_before.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data_before.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata_before.tec.dat"), m_step+1, m_time);

		// dynamic refinement
		for (int level=minlevel; level<maxlevel; level++) {
			for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
				this->test_block_refine(iblock);
			}

			m_tree->regrid();

			//
			m_tree->initNewBlockDataAfterRegrid();

			m_fillpatch->restrictAll(Comp_Phi, NumComp);
		}

		m_fillpatch->fillBoundary(0, NumComp, 1);
		// cell->node
		varloc_interp.interp_cell_to_node(*m_node_data, *m_data, Comp_Phi, Comp_Phi, NumComp, 0);

		// save tree structure
		m_tree->writeTreeLeafState(TESTDIR("/amr_leaf.tec.dat"), m_step+1, m_time);
		// save data
		m_tree->writeTreeLeafBlockGrid(TESTDIR("/amr_leaf_grid.tec.dat"), m_step+1, m_time);
		m_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_data.tec.dat"), m_step+1, m_time);
		m_node_data->writeTreeLeafBlockData(TESTDIR("/amr_leaf_nodedata.tec.dat"), m_step+1, m_time);

		//for (int i=0; i<m_tree->numBlocks; i++) {
		//	//if ((*m_tree)[i].isLeaf()) 
		//	{
		//		char filename[64];
		//		sprintf(filename, "test3/amr_block%04d.tec.dat", i);
		//		m_data->writeBlockData(i, filename, m_step+1, m_time);
		//	}
		//}

	}


protected:

	void init_tree() {
		m_tree = new AmrTree;
		//
		m_tree->minAmrLevel = minlevel;
		m_tree->maxAmrLevel = maxlevel;
		//
		m_tree->blockValidIndex = validBox;
		m_tree->blockGrownIndex = grownBox;
		m_tree->blockNumGrow = ng;
		
		// allocate memory
		m_tree->init();

		{ // setup root level
			const int root_level = 1;
			const int nroot = n0x * n0y * n0z;

			m_tree->numBlocks = nroot;

			for (int k=0; k<n0z; k++) {
				for (int j=0; j<n0y; j++) {
					for (int i=0; i<n0x; i++) {
						int iroot = i + j*n0x + k*n0x*n0y;

						AmrTreeNode &root = m_tree->blocks[iroot];

						Vector3d block_lo = probbox.lo();
						Vector3d block_hi = probbox.lo();
						block_lo.x += probbox.length(0) / n0x * i;
						block_lo.y += probbox.length(1) / n0y * j;
						block_hi.x += probbox.length(0) / n0x * (i+1);
						block_hi.y += probbox.length(1) / n0y * (j+1);

						RealBox rootbox(block_lo, block_hi);
						root.boundBox = rootbox;
						root.blockCenter = rootbox.center();
						root.blockLength = rootbox.length();

						root.nodeType = TreeNodeType_IsLeaf;
						root.setLevel(root_level);

						// set BC
						for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
							if (iface == 0) {
								if (i == 0) {
									root.neighbor[iface] = PhysBC_Symmetry;
								} else {
									root.neighbor[iface] = (i-1) + j*n0x + k*n0x*n0y;
								}
							} else if (iface == 1) {
								if (i == n0x-1) {
									root.neighbor[iface] = PhysBC_Symmetry;
								} else {
									root.neighbor[iface] = (i+1) + j*n0x + k*n0x*n0y;
								}
							} else if (iface == 2) {
								if (j == 0) {
									root.neighbor[iface] = PhysBC_Symmetry;
								} else {
									root.neighbor[iface] = i + (j-1)*n0x + k*n0x*n0y;
								}
							} else if (iface == 3) {
								if (j == n0y-1) {
									root.neighbor[iface] = PhysBC_Symmetry;
								} else {
									root.neighbor[iface] = i + (j+1)*n0x + k*n0x*n0y;
								}
							}
						} // end loop face
					}
				}
			}

		}

		// static refine up to min level
		if (1) {
			m_tree->initRefineToMinLevel();
		} else {
			// build surrounding list
			m_tree->buildSurroundingBlocks(1);
			// build level cache for root level
			m_tree->cacheTreeLevels();

			for (int level=1; level<minlevel; level++) {
				for (int i=0; i<m_tree->numBlocks; i++) {
					m_tree->refineFlag[i] = 1;
				}
				m_tree->regrid();
			}
		}
	}

	// must be called after tree initialization
	void init_data() {
		m_data = new TreeData(*m_tree, validBox, grownBox, ng, NumComp);
		m_data->setValue(0);
		
		// set value for exising blocks
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			calc_block_value(iblock);
		}

		// allocate node-based data
		IndexBox nodeBox(validBox);
		nodeBox.staggerAll();

		m_node_data = new TreeData(*m_tree, nodeBox, nodeBox, 0, NumComp);
		m_node_data->setValue(0);

	}

	void init_patch() {
		m_bcpatch = new BCPatch(*m_tree, *m_data, VARLOC_CELL);
		for (int comp=0; comp<NumComp; comp++) {
			BCRegister& bcreg = m_bcpatch->boundaryRegister(comp);
			for (SurroundingIndex isurr=0; isurr<SurroundingIndex::NumSurr; isurr++) {
				int ipos = isurr.ipos();
				int jpos = isurr.jpos();
				int kpos = isurr.kpos();

				if (1) {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, Dirichlet);
				} else {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, SimpleFill);
				}

				bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Symmetry, Neumann);
			}
		}

		m_fillpatch = new FillPatch(*m_tree, *m_data/*, VARLOC_CELL*/);
		for (int comp=0; comp<NumComp; comp++) {
			//m_fillpatch->setProlongation(comp, InterpPatch::PROLONG_INJECTION);
			m_fillpatch->setProlongation(comp, InterpPatch::PROLONG_CENTERED_LINEAR);
			//m_fillpatch->setProlongation(comp, InterpPatch::PROLONG_CENTERED_LINEAR_LIMITED);
			m_fillpatch->setRestriction(comp, InterpPatch::RESTRICT_AVERAGE);
		}
		// bind BC
		m_fillpatch->bc_patch = m_bcpatch;
	}

	void init_state() {
		m_state = new TreeStateData;
		// set variables
		m_state->new_data = m_data;
		m_state->bc_patch = m_bcpatch;
		m_state->fill_new_patch = m_fillpatch;

		// add to tree
		m_tree->treeStateData.push_back(*m_state);
	}

	void test_block_refine(int iblock) {
		AmrTreeNode &block = (*m_tree)[iblock];
		int level = block.getLevel();

		if (block.isLeaf() || block.isParentOfLeafChild()) {
			const IndexBox &vbox = m_tree->validBlockCellBox();
			const int nbx = vbox.size(0);
			const int nby = vbox.size(1);
			const int nbz = vbox.size(2);

			const int &ilo = vbox.ilo();
			const int &jlo = vbox.jlo();
			const int &klo = vbox.klo();
			const int &ihi = vbox.ihi();
			const int &jhi = vbox.jhi();
			const int &khi = vbox.khi();

			const RealBox &bbox = (*m_tree)[iblock].boundBox;
			const double dx = bbox.length(0) / nbx;
			const double dy = bbox.length(1) / nby;
			const double dz = bbox.length(2) / nbz;
			const double dh = dx;
			//const double &xlo = bbox.lo().x;
			//const double &ylo = bbox.lo().y;
			//const double &zlo = bbox.lo().z;

			DoubleBlockData &blockdata = (*m_data)[iblock];

			int needRefine = 0;
			int needCoarsen = 1;

			for (int k=klo; k<=khi; k++) {
				//double zz = zlo + ((double) k + 0.5) * dz;
				for (int j=jlo; j<=jhi; j++) {
					//double yy = ylo + ((double) j + 0.5) * dy;
					for (int i=ilo; i<=ihi; i++) {
						//double xx = xlo + ((double) i + 0.5) * dx;

						//Vector3d pos = Vector3d::VecMake(xx,yy,zz);
						double val = blockdata(i,j,k,Comp_Phi);
						
						needRefine = needRefine || (abs(val) <= dh*2.0);
						needCoarsen = needCoarsen && (abs(val) >= dh*4.0);
					}
				}
			}

			if (level < m_tree->maxRefineLevel()) {
				if (needRefine) {
					m_tree->refineFlag[iblock] = 1;
				}
			}

			if (level>m_tree->minAmrLevel && !m_tree->refineFlag[iblock]) {
				if (needCoarsen) {
					m_tree->coarsenFlag[iblock] = 1;
				}
			}
		}
	}

	double calc_value(int comp, const Vector3d &pos) const {
		Vector3d diff = pos - circ_cen;
		double dist = Vec3Length(&diff);
		double value = dist - circ_rad;
		return value;
	}

	void calc_block_value(int iblock) {
		const IndexBox &vbox = m_tree->validBlockCellBox();
		const int nbx = vbox.size(0);
		const int nby = vbox.size(1);
		const int nbz = vbox.size(2);

		const int &ilo = vbox.ilo();
		const int &jlo = vbox.jlo();
		const int &klo = vbox.klo();
		const int &ihi = vbox.ihi();
		const int &jhi = vbox.jhi();
		const int &khi = vbox.khi();

		const RealBox &bbox = (*m_tree)[iblock].boundBox;
		const double dx = bbox.length(0) / nbx;
		const double dy = bbox.length(1) / nby;
		const double dz = bbox.length(2) / nbz;
		const double &xlo = bbox.lo().x;
		const double &ylo = bbox.lo().y;
		const double &zlo = bbox.lo().z;

		DoubleBlockData &blockdata = (*m_data)[iblock];

		for (int comp=0; comp<NumComp; comp++) {
			for (int k=klo; k<=khi; k++) {
				double zz = zlo + ((double) k + 0.5) * dz;
				for (int j=jlo; j<=jhi; j++) {
					double yy = ylo + ((double) j + 0.5) * dy;
					for (int i=ilo; i<=ihi; i++) {
						double xx = xlo + ((double) i + 0.5) * dx;

						Vector3d pos = Vector3d::VecMake(xx,yy,zz);
						double val = calc_value(comp, pos);
						blockdata(i,j,k,comp) = val;
					}
				}
			}
		}
	}
};
#ifdef TESTDIR
#undef TESTDIR
#endif



class Test3
{
	
	enum {
		PhysBC_Fixed = -100,
		PhysBC_Symmetry,
	};

	enum {
		Comp_Phi = 0,
		NumComp,
	};


	AmrTree *m_tree;
	TreeData *m_data;
	BCPatch *m_bcpatch;
	FillPatch *m_fillpatch;

	int m_step;
	double m_time;

	Vector3d probmin, probmax;
	RealBox probbox;

	int nb;
	int nbx, nby, nbz;
	int ng;
	int ngx, ngy, ngz;
	IndexBox validBox, grownBox;

	int minlevel, maxlevel;

public:

	void init(int argc, char *argv[]) {
		m_step = 0;
		m_time = 0;

		probmin = Vector3d::VecMake(0.0, 0.0, 0.0);
		probmax = Vector3d::VecMake(4.0, 4.0, 0.0);
		probbox = RealBox(probmin, probmax);

		minlevel = 2;
		maxlevel = 3;

		nb = 4;
		nbx = nb;
		nby = nb;
		nbz = 1;
		ng = 2;
		ngx = ng;
		ngy = ng;
		ngz = 0;

		validBox = IndexBox(
			Vector3i::VecMake(0,0,0),
			Vector3i::VecMake(nbx-1,nby-1,nbz-1));
		grownBox = IndexBox(
			IndexBox(validBox).extend(ngx,ngy,ngz));

		//
		this->init_tree();
		// save tree structure
		m_tree->writeTreeLeafState("test3/amr_leaf.tec.dat", m_step+1, m_time);
		m_tree->writeTreeLeafBlockGrid("test3/amr_leaf_grid.tec.dat", m_step+1, m_time);

		//
		this->init_data();
		
		//
		this->init_patch();
		//
		m_fillpatch->fillBoundary(0, NumComp, ng);

		// save 
		m_data->writeTreeLeafBlockData("test3/amr_leaf_data.tec.dat", m_step+1, m_time);

		for (int i=0; i<m_tree->numBlocks; i++) {
			//if ((*m_tree)[i].isLeaf()) 
			{
				char filename[64];
				sprintf(filename, "test3/amr_block%04d.tec.dat", i);
				m_data->writeBlockData(i, filename, m_step+1, m_time);
			}
		}

	}


protected:

	int isHot(const Vector3d &vpos) const {
		//return vpos.x < 2.0;
		return (vpos.x<2.0 && vpos.y<2.0) || (vpos.x>2.0 && vpos.y>2.0);
	}

	void init_tree() {
		m_tree = new AmrTree;
		//
		m_tree->minAmrLevel = minlevel;
		m_tree->maxAmrLevel = maxlevel;
		//
		m_tree->blockValidIndex = validBox;
		m_tree->blockGrownIndex = grownBox;
		m_tree->blockNumGrow = ng;
		//
		m_tree->init();

		{ // setup base level
			m_tree->numBlocks = 1;
			AmrTreeNode &root = m_tree->blocks[0];

			root.boundBox = probbox;
			root.blockCenter = probbox.center();
			root.blockLength = probbox.length();

			root.nodeType = TreeNodeType_IsLeaf;
			root.setLevel(1);

			// set BC
			for (FaceIndex iface=0; iface<FaceIndex::NumFace; iface++) {
				if (iface.dir() == 0) {
					root.neighbor[iface] = PhysBC_Fixed;
				} else {
					root.neighbor[iface] = PhysBC_Symmetry;
				}
			}
		}
		// static refine up to min level
		for (int level=1; level<minlevel; level++) {
			for (int i=0; i<m_tree->numBlocks; i++) {
				m_tree->refineFlag[i] = 1;
			}
			m_tree->regrid();
		}
		// specified refine
		for (int level=minlevel; level<maxlevel; level++) {
			for (int i=0; i<m_tree->numBlocks; i++) {
				if (isHot(m_tree->blocks[i].blockCenter)) {
					m_tree->refineFlag[i] = 1;
				}
			}
			m_tree->regrid();
		}
	}

	// must be called after tree initialization
	void init_data() {
		m_data = new TreeData(*m_tree, validBox, grownBox, ng, NumComp);
		m_data->setValue(0);
		
		for (int iblock=0; iblock<m_tree->numBlocks; iblock++) {
			if (isHot(m_tree->blocks[iblock].blockCenter)) {
				m_data->m_data[iblock].setValue(2.0);
			} else {
				m_data->m_data[iblock].setValue(1.0);
			}
		}
	}

	void init_patch() {
		m_bcpatch = new BCPatch(*m_tree, *m_data, VARLOC_CELL);
		for (int comp=0; comp<NumComp; comp++) {
			BCRegister& bcreg = m_bcpatch->boundaryRegister(comp);
			for (SurroundingIndex isurr=0; isurr<SurroundingIndex::NumSurr; isurr++) {
				int ipos = isurr.ipos();
				int jpos = isurr.jpos();
				int kpos = isurr.kpos();

				if (1) {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, Dirichlet);
				} else {
					bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Fixed, SimpleFill);
				}

				bcreg.setBCMap(ipos, jpos, kpos, PhysBC_Symmetry, Neumann);
			}
		}

		m_fillpatch = new FillPatch(*m_tree, *m_data, VARLOC_CELL);
		for (int comp=0; comp<NumComp; comp++) {
			m_fillpatch->setProlongation(comp, InterpPatch::PROLONG_INJECTION);
			//m_fillpatch->setProlongation(comp, InterpPatch::PROLONG_CENTERED_LINEAR);
		}
		// bind BC
		m_fillpatch->bc_patch = m_bcpatch;
	}
};

