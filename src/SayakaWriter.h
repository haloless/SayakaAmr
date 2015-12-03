#pragma once

#include <sstream>

#include "SayakaCommons.h"
#include "SayakaTree.h"
#include "SayakaTreeData.h"

namespace sayaka
{

struct WriterMeta
{
	double time;
	double dt;
	int step;

	std::vector<std::string> name;

public:
	WriterMeta() 
	{
		reset();
	}

	void guessName(const TreeData &data) {
		const int ncomp = data.numComp();
		guessName(ncomp);
	}
	void guessName(const int ncomp) {
		name.resize(ncomp);
		for (int comp=0; comp<ncomp; comp++) {
			if (name[comp].empty()) {
				std::stringstream ss;
				ss << "comp" << comp;
				name[comp] = ss.str();
			}
		}
	}

	void reset() {
		time = 0;
		dt = 1.0;
		step = 0;
		name.clear();
	}
};

int WriteDataHdf5(const TreeData &data, const char *filename, WriterMeta &meta);

int WriteTreeDataVtk(const TreeData &data, const char *filename, WriterMeta &meta);
// Merge leaf blocks to a whole unstructured mesh
int WriteLeafDataVtk(const TreeData &data, const char *filename, WriterMeta &meta);


} // namespace_sayaka

