#pragma once

//#include "stdafx.h"


#include <stdio.h>
#include <stdarg.h>


//---printing and logging
struct CLogPrintf
{
	static FILE* ms_fp;
	static bool open(const char* fname){
		ms_fp=fopen(fname,"w");
		setvbuf(ms_fp,NULL,_IONBF,0);
		return ms_fp!=NULL;
	}
	static void close(){
		if(ms_fp){
			fclose(ms_fp);
		}
		ms_fp=NULL;
	}
	static int LogPrintf(const char* fmt,...){
		va_list ap;
		va_start(ap,fmt);
		// print on screen
		char buffer[512];
		int num=vsprintf(buffer,fmt,ap);
		fputs(buffer,stdout);
		fflush(stdout);

		// save to file
		if(ms_fp){
			fputs(buffer,ms_fp);
			fflush(ms_fp);
		}

		va_end(ap);
		return num;
	}
	static int writeToLog(const char* fmt,...){
		if(!ms_fp)return 0;
		va_list ap;
		va_start(ap,fmt);
		int num=vfprintf(ms_fp,fmt,ap);
		fflush(ms_fp);
		return num;
	}
	static void includeFileInLog(const char* fname){
		if(ms_fp){
			fprintf(ms_fp,"[INSERTED FILE %s]\n---BEGIN---\n",fname);
			FILE *fp=fopen(fname,"r");
			if(fp){
				const int buflen = 256;
				char buffer[buflen];
				for(;;){
					size_t read=fread(buffer,sizeof(buffer[0]),buflen,fp);
					if(read==0)break;
					else fwrite(buffer,sizeof(buffer[0]),read,ms_fp);
				}
				fclose(fp);
			}
			fprintf(ms_fp,"---END---\n");
			fflush(ms_fp);
		}
	}
};

#define LOGPRINTF CLogPrintf::LogPrintf
#define WriteLog CLogPrintf::writeToLog


