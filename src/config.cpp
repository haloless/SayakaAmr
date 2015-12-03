
#include <utility>
#include <algorithm>

//#include "stdafx.h"	//every .cpp/.c file must include this FIRST
#include "config.h"
#include "log.h"


bool beginsWith(const char* line, const char* keyword)
{
	return 0==strncmp(line,keyword,strlen(keyword));
}

bool endsWith(const char* line, const char* keyword)
{
	int length = (int)strlen(line);
	int offset = std::max(0,length - (int)strlen(keyword));
	return 0==strncmp(line + offset,keyword,strlen(keyword));
}

int CONFIG::load_from_file(const char *pFile)
{
	FILE *fp=fopen(pFile,"r");
	if(!fp)return CONFIG_ERROR_FILE_OPEN;

	int ret = load_from_fp(fp);
	fclose(fp);
	return ret;
}

int CONFIG::load_from_fp(FILE* fp,const char *read_until)
{
	//clear old data
	m_Data.clear();

	//read each line
	char buf[256];
	std::string linebuf;
	for(int lineIdx=1;;lineIdx++){
		char *pRes=fgets(buf,sizeof(buf)/sizeof(buf[0])-1,fp);
		//trim CR & LF
		bool bEOL=false;	//true if buffer reaches end of line
		for(;;){
			size_t len=strlen(buf);
			if( buf[len-1]=='\n' || buf[len-1]=='\r'){
				buf[len-1]='\0';
				bEOL=true;
			}else break;
		}
		linebuf+=buf;
		if(!pRes)bEOL=true;

		if(bEOL && !linebuf.empty()){
			if(read_until && beginsWith(linebuf.c_str(),read_until)){
				return 0;
			}
			//split line by '=' or ':'

			if(linebuf[0]!=';'){	//line begins with ';' is a comment line
				std::string strKey;
				size_t len=linebuf.length();

				size_t idx=0;
				for(;idx<len;idx++){
					if(linebuf[idx]=='=' || linebuf[idx]==':'){
						idx++;
						break;
					}else{
						strKey+=linebuf[idx];
					}
				}
				if(idx>=len){
					LOGPRINTF("Parse error: line #%d: \"%s\" missing '=' or ':'\n",lineIdx,linebuf.c_str());
					fclose(fp);
					return CONFIG_ERROR_PARSE;
				}

				//trim spaces in key
				stripString(strKey);
				//make strKey lower
				std::transform(strKey.begin(), strKey.end(), strKey.begin(), tolower);

				//trim spaces in data
				std::string data=linebuf.c_str()+idx;
				stripString(data);
				//store data
				m_Data[strKey]=data;
			}

			linebuf.clear();
		}
		if(!pRes)break;
	}

	return 0;
}

//---Accessors
const char* CONFIG::getData(const char *key)const
{
	std::map<std::string,std::string>::const_iterator ite=m_Data.find(key);
	if(ite==m_Data.end())return NULL;
	else return (*ite).second.c_str();
}

//-----

void stripString(std::string &str)
{
	str.erase(0, str.find_first_not_of(' '));	//prefixing spaces
	str.erase(0, str.find_first_not_of('\t'));	//prefixing tabs
	str.erase(str.find_last_not_of(' ')+1);		//surfixing spaces
	str.erase(str.find_last_not_of('\t')+1);	//surfixing tabs
}

char* trimCRLF(char *buf)
{
	for(;;){
		size_t len=strlen(buf);
		if( buf[len-1]=='\n' || buf[len-1]=='\r'){
			buf[len-1]='\0';
		}else{
			break;
		}
	}
	return buf;
}

int config_from_file(CONFIG *conf, const char *fname)
{
	return conf->load_from_file(fname);
}

int config_from_fp(CONFIG *conf, FILE* fp, const char *read_until)
{
	return conf->load_from_fp(fp,read_until);
}

//getXXX returns parameter; if the parameter is not written in file, default value is used
//getxxxA requires parameter in the file. No default value is allowed.

int getInt(const CONFIG *config, const char* prefix, const char* pKey, int defaultValue)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return atoi(data);
	}else{
		return defaultValue;
	}
}

int getIntA(const CONFIG *config, const char* prefix, const char* pKey)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return atoi(data);
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

void getInt3(const CONFIG *config, const char* prefix, const char* pKey, int *a,int *b,int *c, int defA, int defB, int defC)
{
	*a = defA;
	*b = defB;
	*c = defC;

	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		int ret = sscanf(data, "%d %d %d", a, b, c);
		if(ret!=3){
			ret = sscanf(data, "%d,%d,%d", a, b, c);
			if(ret!=3){
				LOGPRINTF("ERROR!: in '%s', 3 parameters are required\n",(strPrefix + pKey).c_str());
				exit(-1);
			}
		}
	}
}

void getInt3A(const CONFIG *config, const char* prefix, const char* pKey, int *a,int *b,int *c)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		int ret = sscanf(data, "%d %d %d", a, b, c);
		if(ret!=3){
			ret = sscanf(data, "%d,%d,%d", a, b, c);
			if(ret!=3){
				LOGPRINTF("ERROR!: in '%s', 3 parameters are required\n",(strPrefix + pKey).c_str());
				exit(-1);
			}
		}
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

bool getOnOff(const CONFIG *config, const char* prefix, const char* pKey, bool defaultValue)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		if(stricmp(data,"yes")==0 || stricmp(data,"true")==0 || stricmp(data,"on")==0){
			return true;
		}else if(stricmp(data,"no")==0 || stricmp(data,"false")==0 || stricmp(data,"off")==0){
			return false;
		}else{
			// 0 = false, otherwise true
			return atoi(data) != 0;
		}
	}else{
		return defaultValue;
	}
}

bool getOnOffA(const CONFIG *config, const char* prefix, const char* pKey)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		if(stricmp(data,"yes")==0 || stricmp(data,"true")==0 || stricmp(data,"on")==0){
			return true;
		}else if(stricmp(data,"no")==0 || stricmp(data,"false")==0 || stricmp(data,"off")==0){
			return false;
		}else{
			// 0 = false, otherwise true
			return atoi(data) != 0;
		}
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

double getDouble(const CONFIG *config, const char* prefix, const char* pKey, double defaultValue)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return atof(data);
	}else{
		return defaultValue;
	}
}

double getDoubleA(const CONFIG *config, const char* prefix, const char* pKey)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return atof(data);
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

void getDouble3(const CONFIG *config, const char* prefix, const char* pKey, double *a, double *b, double *c, double defA, double defB, double defC)
{
	*a = defA;
	*b = defB;
	*c = defC;

	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		int ret = sscanf(data, "%lf %lf %lf", a, b, c);
		if(ret!=3){
			ret = sscanf(data, "%lf,%lf,%lf", a, b, c);
			if(ret!=3){
				LOGPRINTF("ERROR!: in '%s', 3 parameters are required\n",(strPrefix + pKey).c_str());
				exit(-1);
			}
		}
	}
}

void getDouble3A(const CONFIG *config, const char* prefix, const char* pKey, double *a, double *b, double *c)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		int ret = sscanf(data, "%lf %lf %lf", a, b, c);
		if(ret!=3){
			ret = sscanf(data, "%lf,%lf,%lf", a, b, c);
			if(ret!=3){
				LOGPRINTF("ERROR!: in '%s', 3 parameters are required\n",(strPrefix + pKey).c_str());
				exit(-1);
			}
		}
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

const char* getText(const CONFIG *config, const char* prefix, const char* pKey, const char* defaultValue)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return data;
	}else{
		return defaultValue;
	}
}


const char* getTextA(const CONFIG *config, const char* prefix, const char* pKey)
{
	std::string strPrefix = prefix ? prefix : "";
	const char* data=config->getData((strPrefix + pKey).c_str());
	if(data){
		return data;
	}else{
		LOGPRINTF("ERROR!: '%s' is missing in the file\n",(strPrefix + pKey).c_str());
		exit(-1);
	}
}

