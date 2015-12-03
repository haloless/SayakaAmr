
#include "SayakaWriter.h"

//#define SAYAKA_USE_VTK

#ifdef SAYAKA_USE_VTK
// our helper header
#	define VTKHELPER_LINKER_PRAGMA
#	include <vtkHelper.hpp>
// VTK
#	include <vtkCellData.h>
#	include <vtkPointData.h>
#	include <vtkUniformGrid.h>
#	include <vtkAMRBox.h>
#	include <vtkHierarchicalBoxDataSet.h>
#	include <vtkXMLHierarchicalBoxDataWriter.h>
#	include <vtkMultiBlockDataSet.h>
#	include <vtkXMLMultiBlockDataWriter.h>
//
#	include <vtkAppendFilter.h>
#	include <vtkUnstructuredGrid.h>
#	include <vtkXMLUnstructuredGridWriter.h>
#endif


namespace sayaka
{


#ifdef SAYAKA_USE_VTK
static vtkUniformGrid* GetBlockGrid(int iblock,
	const TreeData &thedata, const WriterMeta &meta)
{
	const AmrTree &tree = thedata.tree();
	const IndexBox &cellbox = tree.validBlockCellBox();
	// cell 
	const int ncx = cellbox.size(0);
	const int ncy = cellbox.size(1);
	const int ncz = cellbox.size(2);
	const int ncell = ncx * ncy * ncz;
	// node
	const int npx = ncx + XDIM;
	const int npy = ncy + YDIM;
	const int npz = ncz + ZDIM;
	const int nnode = npx * npy * npz;

	const IndexBox &databox = thedata.validBox();
	decl_box_range(databox, i,j,k);
	const int ndata = databox.stride();

	const int ncomp = thedata.numComp();

	const AmrTreeNode &block = tree[iblock];
	const DoubleBlockData &bdata = thedata[iblock];

	//
	const int ilevel = block.getLevel();

	// low corner of the block
	const RealBox &bbox = block.boundBox;
	const Vector3d &bxlo = bbox.lo();
	const Vector3d dh = tree.getBlockCellSize(iblock);

	//
	//NEW_VTKOBJ(vtkUniformGrid, vtkgrid);
	vtkUniformGrid *vtkgrid = vtkUniformGrid::New();
	vtkgrid->Initialize();
	vtkgrid->SetOrigin(bxlo.x, bxlo.y, bxlo.z);
	vtkgrid->SetSpacing(dh.x, dh.y, dh.z);
	vtkgrid->SetDimensions(npx, npy, npz);

	//
	NEW_VTKOBJ(vtkIntArray, levelArray);
	levelArray->SetName("AmrLevel");
	levelArray->SetNumberOfComponents(1);
	levelArray->SetNumberOfTuples(ncell);
	vtkgrid->GetCellData()->AddArray(levelArray);
	//
	for (int k=cellbox.klo(); k<=cellbox.khi(); k++) {
		for (int j=cellbox.jlo(); j<=cellbox.jhi(); j++) {
			for (int i=cellbox.ilo(); i<=cellbox.ihi(); i++) {
				int idx = cellbox.offset(i,j,k);
				assert(0<=idx && idx<ncell);
				levelArray->SetTuple1(idx, ilevel);
			}
		}
	}

	for (int comp=0; comp<ncomp; comp++) {
		//vtkDoubleArray *phiArray = vtkDoubleArray::New();
		NEW_VTKOBJ(vtkDoubleArray, phiArray);
		phiArray->SetName(meta.name[comp].c_str());
		phiArray->SetNumberOfComponents(1);
		phiArray->SetNumberOfTuples(ndata);
		if (thedata.isCellData()) {
			vtkgrid->GetCellData()->AddArray(phiArray);
		} else if (thedata.isNodeData()) {
			vtkgrid->GetPointData()->AddArray(phiArray);
		} else {
			LOGPRINTF("%s: unsupported data\n", __FUNCTION__);
		}
		// collect this array
		//parrays.push_back(phiArray);

		//
		for (int k=klo; k<=khi; k++) {
			for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					int idx = databox.offset(i,j,k);
					assert(0<=idx && idx<ndata);
					phiArray->SetTuple1(idx, bdata(i,j,k,comp));
				}
			}
		}
	}

	//
	vtkgrid->Update();

	// 
	//for (int comp=0; comp<ncomp; comp++) {
	//	parrays[comp]->Delete();
	//}
	return vtkgrid;
}
#endif

int WriteTreeDataVtk(const TreeData &thedata, 
	const char *filename, WriterMeta &meta) 
{
	int ret = 0;

#ifdef SAYAKA_USE_VTK
	const AmrTree &tree = thedata.tree();
	const int finest_ilevel = tree.currentFinestLevel();
	//const int finest_ilevel = 1;

	// prepare
	meta.guessName(thedata);
	//
	//const int ncomp = thedata.numComp();

	const IndexBox &cellbox = tree.validBlockCellBox();
	//// cell 
	//const int ncx = cellbox.size(0);
	//const int ncy = cellbox.size(1);
	//const int ncz = cellbox.size(2);
	//const int ncell = ncx * ncy * ncz;
	//// node
	//const int npx = ncx + XDIM;
	//const int npy = ncy + YDIM;
	//const int npz = ncz + ZDIM;
	//const int nnode = npx * npy * npz;

	//const IndexBox &databox = thedata.validBox();
	//decl_box_range(databox, i,j,k);
	//const int ndata = databox.stride();

	// find the low corner of the problem
	double problo[MAX_DIM] = { 0 };
	for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		const RealBox &bbox = tree[iblock].boundBox;
		if (iblock == 0) {
			for (int dir=0; dir<NDIM; dir++) {
				problo[dir] = bbox.lo()(dir);
			}
		} else {
			for (int dir=0; dir<NDIM; dir++) {
				problo[dir] = std::min(problo[dir], bbox.lo()(dir));
			}
		}
	}

	//
	NEW_VTKOBJ(vtkHierarchicalBoxDataSet, hierBoxDataSet);
	hierBoxDataSet->Initialize();
	
	// set total level number
	hierBoxDataSet->SetNumberOfLevels(finest_ilevel);
	// set block number on each level
	for (int ilevel=1; ilevel<=finest_ilevel; ilevel++) {
		const AmrTreeLevel &level = tree.getTreeLevel(ilevel);
		assert(!level.isEmptyLevel());

		const int ngrid = level.numLevelBlock;
		hierBoxDataSet->SetNumberOfDataSets(ilevel-1, ngrid);
		hierBoxDataSet->SetRefinementRatio(ilevel-1, 2);
	}

	// 
	for (int ilevel=1; ilevel<=finest_ilevel; ilevel++) {
		const AmrTreeLevel &level = tree.getTreeLevel(ilevel);
		assert(!level.isEmptyLevel());

		for (int igrid=0; igrid<level.numLevelBlock; igrid++) {
			const int iblock = level[igrid];
			const Vector3d dh = tree.getBlockCellSize(iblock);
			// calculate the logical index range of block
			int ilo[MAX_DIM] = { 0 };
			int ihi[MAX_DIM] = { 0 };

			const RealBox &bbox = tree[iblock].boundBox;
			const Vector3d &bxlo = bbox.lo();

			for (int dir=0; dir<NDIM; dir++) {
				double xc = bxlo(dir) + 0.5*dh(dir);
				int idx = (int) ((xc-problo[dir]) / dh(dir));
				ilo[dir] = idx;
				ihi[dir] = ilo[dir] + cellbox.size(dir) - 1;
			}

			vtkUniformGrid *vtkgrid = GetBlockGrid(iblock, thedata, meta);

			hierBoxDataSet->SetDataSet(ilevel-1, igrid, 
				ilo, ihi, vtkgrid);

			vtkgrid->Delete();
		}
	}
	
	hierBoxDataSet->Update();

	//
	
	NEW_VTKOBJ(vtkXMLHierarchicalBoxDataWriter, writer);
	writer->SetFileName(filename);
	writer->SetInput(hierBoxDataSet);
	if (!writer->Write()) {
		LOGPRINTF("%s: failed to write %s\n", __FUNCTION__, filename);
		ret = 1;
	} else {
		LOGPRINTF("%s: saved %s\n", __FUNCTION__, filename);
		ret = 0;
	}
#else
	LOGPRINTF("%s: not available\n", __FUNCTION__);
	//exit(1);
#endif

	return ret;
}

int WriteLeafDataVtk(const TreeData &thedata, 
	const char *filename, WriterMeta &meta) 
{
	int ret = 0;

#ifdef SAYAKA_USE_VTK
	const AmrTree &tree = thedata.tree();

	// prepare
	meta.guessName(thedata);
	//
	const int ncomp = thedata.numComp();
	//const int ncomp = 1;

	const IndexBox &cellbox = tree.validBlockCellBox();
	// cell 
	const int ncx = cellbox.size(0);
	const int ncy = cellbox.size(1);
	const int ncz = cellbox.size(2);
	const int ncell = ncx * ncy * ncz;
	// node
	const int npx = ncx + XDIM;
	const int npy = ncy + YDIM;
	const int npz = ncz + ZDIM;
	const int nnode = npx * npy * npz;

	const IndexBox &databox = thedata.validBox();
	decl_box_range(databox, i,j,k);
	const int ndata = databox.stride();


	//
	NEW_VTKOBJ(vtkAppendFilter, vtkappender);
	vtkappender->MergePointsOn();
	//vtkappender->set

	for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		const AmrTreeNode &block = tree[iblock];
		const DoubleBlockData &bdata = thedata[iblock];

		if (block.isLeaf()) {
			//
			const int ilevel = block.getLevel();

			// low corner of the block
			const RealBox &bbox = block.boundBox;
			const Vector3d &bxlo = bbox.lo();
			const Vector3d dh = tree.getBlockCellSize(iblock);
			
			//
			NEW_VTKOBJ(vtkUniformGrid, vtkgrid);
			vtkgrid->Initialize();
			vtkgrid->SetOrigin(bxlo.x, bxlo.y, bxlo.z);
			vtkgrid->SetSpacing(dh.x, dh.y, dh.z);
			vtkgrid->SetDimensions(npx, npy, npz);

			//
			NEW_VTKOBJ(vtkIntArray, levelArray);
			levelArray->SetName("AmrLevel");
			levelArray->SetNumberOfComponents(1);
			levelArray->SetNumberOfTuples(ncell);
			vtkgrid->GetCellData()->AddArray(levelArray);
			//
			for (int k=cellbox.klo(); k<=cellbox.khi(); k++) {
			for (int j=cellbox.jlo(); j<=cellbox.jhi(); j++) {
			for (int i=cellbox.ilo(); i<=cellbox.ihi(); i++) {
				int idx = cellbox.offset(i,j,k);
				assert(0<=idx && idx<ncell);
				//levelArray->SetValue(idx, ilevel);
				levelArray->SetTuple1(idx, ilevel);
			}
			}
			}
			

			std::vector<vtkDoubleArray*> parrays;

			for (int comp=0; comp<ncomp; comp++) {
				vtkDoubleArray *phiArray = vtkDoubleArray::New();
				phiArray->SetName(meta.name[comp].c_str());
				phiArray->SetNumberOfComponents(1);
				phiArray->SetNumberOfTuples(ndata);
				if (thedata.isCellData()) {
					vtkgrid->GetCellData()->AddArray(phiArray);
				} else if (thedata.isNodeData()) {
					vtkgrid->GetPointData()->AddArray(phiArray);
				} else {
					LOGPRINTF("%s: unsupported data\n", __FUNCTION__);
				}
				// collect this array
				parrays.push_back(phiArray);

				//
				for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					int idx = databox.offset(i,j,k);
					assert(0<=idx && idx<ndata);
					//phiArray->SetValue(idx, bdata(i,j,k,comp));
					phiArray->SetTuple1(idx, bdata(i,j,k,comp));
				}
				}
				}
			}

			//
			vtkgrid->Update();
			vtkappender->AddInput(vtkgrid);

			// 
			for (int comp=0; comp<ncomp; comp++) {
				parrays[comp]->Delete();
			}
		} // is leaf block
	} // end loop blocks

	if (vtkappender->GetNumberOfInputConnections(0) > 0) {
		vtkappender->Update();
	}

	//
	NEW_VTKOBJ(vtkUnstructuredGrid, output);
	output->ShallowCopy(vtkappender->GetOutput());

	//
	vtkHelper_setGridTime(output, meta.time);

	//
	
	NEW_VTKOBJ(vtkXMLUnstructuredGridWriter, writer);
	writer->SetFileName(filename);
	writer->SetInput(output);
	if (!writer->Write()) {
		LOGPRINTF("%s: failed to write %s\n", __FUNCTION__, filename);
		ret = 1;
	} else {
		LOGPRINTF("%s: saved %s\n", __FUNCTION__, filename);
		ret = 0;
	}
#else
	LOGPRINTF("%s: not available\n", __FUNCTION__);
	//exit(1);
#endif
	return ret;
}

} // namespace_sayaka



