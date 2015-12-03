#pragma once

#include <cstdio>

#include <string>
#include <map>

#define CONFIG_ERROR_FILE_OPEN 1
#define CONFIG_ERROR_PARSE 2

class CONFIG{
protected:
	std::map<std::string,std::string> m_Data;
public:
	CONFIG(){}
	virtual ~CONFIG(){}

	//---Reads config file
	int load_from_file(const char *pFile);
	int load_from_fp(FILE* fp,const char *read_until=NULL);

	//---Accessors
	const char* getData(const char *key)const;
};

//----
int config_from_file(CONFIG *conf, const char *fname);
int config_from_fp(CONFIG *conf, FILE* fp, const char *read_until=NULL);

//getXXX returns parameter; if the parameter is not written in file, default value is used
//getxxxA requires parameter in the file. No default value is allowed.
int getInt(const CONFIG *config, const char* prefix, const char* pKey, int defaultValue);
int getIntA(const CONFIG *config, const char* prefix, const char* pKey);
void getInt3(const CONFIG *config, const char* prefix, const char* pKey, int *a,int *b,int *c, int defA, int defB, int defC);
void getInt3A(const CONFIG *config, const char* prefix, const char* pKey, int *a,int *b,int *c);

bool getOnOff(const CONFIG *config, const char* prefix, const char* pKey, bool defaultValue);
bool getOnOffA(const CONFIG *config, const char* prefix, const char* pKey);

double getDouble(const CONFIG *config, const char* prefix, const char* pKey, double defaultValue);
double getDoubleA(const CONFIG *config, const char* prefix, const char* pKey);
void getDouble3(const CONFIG *config, const char* prefix, const char* pKey, double *a, double *b, double *c, double defA, double defB, double defC);
void getDouble3A(const CONFIG *config, const char* prefix, const char* pKey, double *a, double *b, double *c);

const char* getText(const CONFIG *config, const char* prefix, const char* pKey, const char* defaultValue="");
const char* getTextA(const CONFIG *config, const char* prefix, const char* pKey);


void stripString(std::string &str);
bool startswith(const std::string &str, const char* head);

//removes "\n" at the end of "buf"
char* trimCRLF(char *buf);

bool beginsWith(const char* line, const char* keyword);
bool endsWith(const char* line, const char* keyword);
