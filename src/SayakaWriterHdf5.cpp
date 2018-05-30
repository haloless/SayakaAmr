

#include "SayakaWriter.h"


#define SAYAKA_USE_HDF

#ifdef SAYAKA_USE_HDF
#	define H5_USE_16_API
#	include <hdf5.h>
#	pragma comment(lib, "libhdf5.lib")
//#	pragma comment(lib, "hdf5_cpp.lib")
//#	pragma comment(lib, "hdf5_tools.lib")
#endif

SAYAKA_NS_BEGIN;


//static const int verbose = 1;
static const int verbose = 0;

struct realvect_t {
	double x;
	double y;
#if (SAYAKA_SPACEDIM == 3)
	double z;
#endif
};

struct prob_domain_t {
	int lo_i;
	int lo_j;
#if (SAYAKA_SPACEDIM == 3)
	int lo_k;
#endif
	int hi_i;
	int hi_j;
#if (SAYAKA_SPACEDIM == 3)
	int hi_k;
#endif
};

static hid_t create_file_hdf5(const char *filename) {
	hid_t fileid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	return fileid;
}

static hid_t create_prob_domain_type_id() {
	hid_t prob_domain_type_id = H5Tcreate(H5T_COMPOUND, sizeof(prob_domain_t));

	H5Tinsert(prob_domain_type_id, "lo_i", HOFFSET(prob_domain_t, lo_i), H5T_NATIVE_INT);
	H5Tinsert(prob_domain_type_id, "lo_j", HOFFSET(prob_domain_t, lo_j), H5T_NATIVE_INT);
#if (SAYAKA_SPACEDIM == 3)
	H5Tinsert(prob_domain_type_id, "lo_k", HOFFSET(prob_domain_t, lo_k), H5T_NATIVE_INT);
#endif

	H5Tinsert(prob_domain_type_id, "hi_i", HOFFSET(prob_domain_t, hi_i), H5T_NATIVE_INT);
	H5Tinsert(prob_domain_type_id, "hi_j", HOFFSET(prob_domain_t, hi_j), H5T_NATIVE_INT);
#if (SAYAKA_SPACEDIM == 3)
	H5Tinsert(prob_domain_type_id, "hi_k", HOFFSET(prob_domain_t, hi_k), H5T_NATIVE_INT);
#endif

	return prob_domain_type_id;
}

static hid_t create_reat_vect_type_id() {
	hid_t realvect_id = H5Tcreate(H5T_COMPOUND, sizeof(realvect_t));
	H5Tinsert(realvect_id, "x", HOFFSET(realvect_t, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(realvect_id, "y", HOFFSET(realvect_t, y), H5T_NATIVE_DOUBLE);
#if (SAYAKA_SPACEDIM == 3)
	H5Tinsert(realvect_id, "z", HOFFSET(realvect_t, z), H5T_NATIVE_DOUBLE);
#endif
	return realvect_id;
}

static void write_header_hdf5(hid_t file_id,
	int num_components, 
	int num_levels, 
	const int nbox_level[],
	const Vector3d dx_level[],
	const double problo[], const double probhi[],
	const AmrTree &tree, const TreeData &data,
	const WriterMeta &meta) 
{
	//
	const IndexBox &cellbox = tree.validBlockCellBox();
	//const int ncell = cellbox.stride();
	const IndexBox &databox = data.validBox();
	const int ndata = databox.stride();

	//
	hid_t group_id = H5Gopen(file_id, "/");

	{ // number of components
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "num_components", 
			H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
		H5Awrite(attribute_id, H5T_NATIVE_INT, &num_components);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}
	{ // number of levels
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "num_levels", 
			H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
		H5Awrite(attribute_id, H5T_NATIVE_INT, &num_levels);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}
	{ // iteration
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "iteration", 
			H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
		int temp[1] = { meta.step };
		H5Awrite(attribute_id, H5T_NATIVE_INT, temp);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}
	{ // time (dummy)
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "time", 
			H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
		double temp[1] = { meta.time };
		H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, temp);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}

#if (1)
	{ // anisotropic stretch
		hid_t realvect_id = create_reat_vect_type_id();

		// 
		realvect_t ansi_stretch;
		ansi_stretch.x = 1.0;
		ansi_stretch.y = 1.0;
#if (SAYAKA_SPACEDIM == 3)
		ansi_stretch.z = 1.0;
#endif

		// we use block 0 as reference
		const Vector3d dh = tree.getBlockCellSize(0);
		ansi_stretch.x = 1.0;
		ansi_stretch.y = dh.y / dh.x;
#if (SAYAKA_SPACEDIM == 3)
		ansi_stretch.z = dh.z / dh.x;
#endif

		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "aspect_ratio", 
			realvect_id, dataspace_id, H5P_DEFAULT);
		H5Awrite(attribute_id, realvect_id, &ansi_stretch);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);

		H5Tclose(realvect_id);
	}
#endif

	{ // cell- or node-centered data?
		int data_centering = 0;
		if (data.isCellData()) {
			data_centering = 0;
		} else if (data.isNodeData()) {
			data_centering = 7;
		} else {
			LOGPRINTF("%s: only cell or node data!!\n", __FUNCTION__);
			exit(1);
		}

		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "data_centering", 
			H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
		int temp[1] = { data_centering };
		H5Awrite(attribute_id, H5T_NATIVE_INT, temp);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}
	{ // problem low corner
		hid_t realvect_id = create_reat_vect_type_id();

		// 
		realvect_t prob_lo;
		prob_lo.x = problo[0];
		prob_lo.y = problo[1];
#if (SAYAKA_SPACEDIM == 3)
		prob_lo.z = problo[2];
#endif

		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "prob_lo", 
			realvect_id, dataspace_id, H5P_DEFAULT);
		H5Awrite(attribute_id, realvect_id, &prob_lo);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);

		H5Tclose(realvect_id);
	}

	// components names
	for (int i=0; i<num_components; i++) {
		const int namelen = 64;
		char name[namelen];
		sprintf(name, "component_%i", i);

		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t string_type = H5Tcopy(H5T_C_S1);
		H5Tset_size(string_type, namelen);
		hid_t attribute_id = H5Acreate(group_id, name, 
			string_type, dataspace_id, H5P_DEFAULT);

		const std::string &compname = meta.name[i];
		H5Awrite(attribute_id, string_type, compname.c_str());

		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
		H5Tclose(string_type);
	}

	H5Gclose(group_id); // group "/"

	// 
	group_id = H5Gcreate(file_id, "/Chombo_global", 0);

	{ // space dimension
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "SpaceDim",
			H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
		int ndim[1] = { sayaka::NDIM };
		H5Awrite(attribute_id, H5T_NATIVE_INT, ndim);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}
	{ // test real (dummy)
		hid_t dataspace_id = H5Screate(H5S_SCALAR);
		hid_t attribute_id = H5Acreate(group_id, "testReal",
			H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
		double temp[1] = { 0 };
		H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, temp);
		H5Sclose(dataspace_id);
		H5Aclose(attribute_id);
	}

	H5Gclose(group_id); // group "/Chombo_global"

	// groups for each level
	for (int ilevel=0; ilevel<num_levels; ilevel++) {
		//const int level = ilevel + 1; // exact level for KyokoAmr

		char name[32];
		sprintf(name, "/level_%i", ilevel);
		group_id = H5Gcreate(file_id, name, 0);

		{ // refine ratio
			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			hid_t attribute_id = H5Acreate(group_id, "ref_ratio",
				H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
			int temp[1] = { 2 };
			H5Awrite(attribute_id, H5T_NATIVE_INT, temp);
			H5Sclose(dataspace_id);
			H5Aclose(attribute_id);
		}
		{ // timestep 
			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			hid_t attribute_id = H5Acreate(group_id, "dt",
				H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
			double temp[1] = { meta.dt };
			H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, temp);
			H5Sclose(dataspace_id);
			H5Aclose(attribute_id);
		}
		{ // cell size
			const double *dx = dx_level[ilevel].data();
			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			hid_t attribute_id = H5Acreate(group_id, "dx",
				H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
			H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, dx);
			H5Sclose(dataspace_id);
			H5Aclose(attribute_id);
		}

		{ // problem domain
			hid_t prob_domain_type_id = create_prob_domain_type_id();

			{ // logical indexing space
				prob_domain_t domain;
				const Vector3d dh = dx_level[ilevel];
				domain.lo_i = 0; 
				domain.hi_i = (int) ((probhi[0]-problo[0]-0.5*dh.x) / dh.x);
				domain.lo_j = 0; 
				domain.hi_j = (int) ((probhi[1]-problo[1]-0.5*dh.y) / dh.y);
#if (SAYAKA_SPACEDIM == 3)
				domain.lo_k = 0; 
				domain.hi_k = (int) ((probhi[2]-problo[2]-0.5*dh.z) / dh.z);
#endif

				hid_t dataspace_id = H5Screate(H5S_SCALAR);
				hid_t attribute_id = H5Acreate(group_id, "prob_domain", 
					prob_domain_type_id, dataspace_id, H5P_DEFAULT);

				H5Awrite(attribute_id, prob_domain_type_id, &domain);

				H5Sclose(dataspace_id);
				H5Aclose(attribute_id);

				if (verbose) {
					LOGPRINTF("level=%d, ilo=%d, jlo=%d, ihi=%d, jhi=%d\n",
						ilevel, domain.lo_i, domain.lo_j, domain.hi_i, domain.hi_j);
				}
			}

			{ // boxes
				hsize_t dims[1] = { nbox_level[ilevel] };
				hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
				hid_t dset_id = H5Dcreate(group_id, "boxes", 
					prob_domain_type_id, dataspace_id, H5P_DEFAULT);
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
			}

			H5Tclose(prob_domain_type_id);
		}
		{ // data:datatype=0
			hsize_t dims[1] = { ndata * num_components * nbox_level[ilevel] };
			hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
			hid_t dset_id = H5Dcreate(group_id, "data:datatype=0", 
				H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
			H5Dclose(dset_id);
			H5Sclose(dataspace_id);
		}

		{ // group2
			hid_t group_id2 = H5Gcreate(group_id, "data_attributes", 0);

			hid_t dataspace_id = H5Screate(H5S_SCALAR);
			hid_t attribute_id = H5Acreate(group_id2, "comps", 
				H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
			//
			H5Awrite(attribute_id, H5T_NATIVE_INT, &num_components);
			//
			H5Sclose(dataspace_id);
			H5Aclose(attribute_id);

			H5Gclose(group_id2);
		}

		H5Gclose(group_id); // group level i
	} // end loop levels
}

int WriteDataHdf5(const TreeData &data, 
	const char *filename, WriterMeta &meta) 
{
	const AmrTree &tree = data.tree();

	const IndexBox &databox = data.validBox();
	decl_box_range(databox, i,j,k);
	//
	const int ndata = databox.stride();
	const int ncomp = data.numComp();
	const int finest_level = tree.currentFinestLevel();
	//const int num_levels = std::min(finest_level, 4);
	const int num_levels = finest_level;

	// prepare data
	meta.guessName(data);

	// find level cell size
	std::vector<Vector3d> dx_level(num_levels);
	dx_level[0] = tree.getBlockCellSize(0);
	for (int ilevel=1; ilevel<num_levels; ilevel++) {
		double ratio = 2.0;
		dx_level[ilevel] = dx_level[ilevel-1];
		for (int dir=0; dir<NDIM; dir++) {
			dx_level[ilevel](dir) /= ratio;
		}
	}

	// find level number
	std::vector<int> nbox_level;
	for (int ilevel=0; ilevel<num_levels; ilevel++) {
		int ngrid = tree.getTreeLevel(ilevel+1).numLevelBlock;
		nbox_level.push_back(ngrid);
		if (verbose) {
			LOGPRINTF("level=%d, ngrid=%d\n", ilevel, ngrid);
		}
	}

	// find the problem size
	double problo[NDIM] = { 0 };
	double probhi[NDIM] = { 0 };
	for (int iblock=0; iblock<tree.numBlocks; iblock++) {
		const RealBox &bbox = tree[iblock].boundBox;
		if (iblock == 0) {
			for (int dir=0; dir<NDIM; dir++) {
				problo[dir] = bbox.lo()(dir);
				probhi[dir] = bbox.hi()(dir);
			}
		} else {
			for (int dir=0; dir<NDIM; dir++) {
				problo[dir] = std::min(problo[dir], bbox.lo()(dir));
				probhi[dir] = std::max(probhi[dir], bbox.hi()(dir));
			}
		}
	}


	// open HDF5
	hid_t file_id = create_file_hdf5(filename);

	// header
	write_header_hdf5(file_id, 
		ncomp, num_levels, 
		&nbox_level[0], &dx_level[0],
		problo, probhi,
		tree, data, meta);
	if (verbose) {
		LOGPRINTF("%s: write header\n", __FUNCTION__);
	}

	// used to create box structure
	hid_t prob_domain_type_id = create_prob_domain_type_id();


	// loop over levels
	for (int level=1; level<=num_levels; level++) {
		const AmrTreeLevel &treelevel = tree.getTreeLevel(level);

		const int ilevel = level - 1;
		const double *dh = dx_level[level-1].data();

		// open group for this level
		char name[32]; sprintf(name, "/level_%i", ilevel);
		hid_t group_id = H5Gopen(file_id, name);

		for (int igrid=0; igrid<treelevel.numLevelBlock; igrid++) {
			const int iblock = treelevel[igrid];
			const AmrTreeNode &block = tree[iblock];
			const DoubleBlockData &db = data[iblock];
			
			{// open boxes dataset
				hsize_t dims[1] = { 1 };
				hid_t space_id = H5Screate_simple(1, dims, NULL);
				hid_t dset_id = H5Dopen(group_id, "boxes");

				// single grid
				hsize_t count[1] = { 1 };
				hsize_t offset[1] = { igrid };
				hid_t filespace = H5Dget_space(dset_id);
				H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
				
				// grid extent
				const RealBox &bbox = block.boundBox;
				const Vector3d &bxlo = bbox.lo();
				const Vector3d &bxhi = bbox.hi();

				int box_lo[NDIM], box_hi[NDIM];
				for (int dir=0; dir<NDIM; dir++) {
					box_lo[dir] = (int) ((bxlo(dir)+0.5*dh[dir]-problo[dir]) / dh[dir]);
					box_hi[dir] = (int) ((bxhi(dir)-0.5*dh[dir]-problo[dir]) / dh[dir]);
				}

				prob_domain_t box;
				box.lo_i = box_lo[0]; box.hi_i = box_hi[0];
				box.lo_j = box_lo[1]; box.hi_j = box_hi[1];
#if (SAYAKA_SPACEDIM == 3)
				box.lo_k = box_lo[2]; box.hi_k = box_hi[2];
#endif

				H5Dwrite(dset_id, prob_domain_type_id, space_id, filespace, H5P_DEFAULT, &box);

				H5Sclose(space_id);
				H5Dclose(dset_id);
				H5Sclose(filespace);
			}

			// begin treat data
			hsize_t count[1] = { ndata };
			hsize_t offset[1] = { igrid * ncomp * ndata };

			// 
			for (int ivar=0; ivar<ncomp; ivar++) {
				// open 'data:datatype=0' dataset
				hsize_t dims[1] = { ndata };
				hid_t space_id = H5Screate_simple(1, dims, NULL);
				hid_t dset_id = H5Dopen(group_id, "data:datatype=0");

				hid_t filespace = H5Dget_space(dset_id);
				H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

				// increment for next component
				offset[0] += count[0];

				// collect variable
				std::vector<double> unk(ndata, 0.0);
				for (int k=klo; k<=khi; k++) {
				for (int j=jlo; j<=jhi; j++) {
				for (int i=ilo; i<=ihi; i++) {
					int idx = databox.offset(i,j,k);
					unk[idx] = db(i,j,k,ivar);
				}
				}
				}

				H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, space_id, filespace, H5P_DEFAULT, &unk[0]);

				H5Sclose(space_id);
				H5Dclose(dset_id);
				H5Sclose(filespace);
			} // end data
		} // end loop grids at level

		// close level group
		H5Gclose(group_id);

		if (verbose) {
			LOGPRINTF("%s: level=%d/%d\n", __FUNCTION__, level, num_levels);
		}
	}

	//
	H5Tclose(prob_domain_type_id);

	// remember to close
	H5Fclose(file_id);

	LOGPRINTF("%s: write %s\n", __FUNCTION__, filename);

	return 0;
}


SAYAKA_NS_END;


