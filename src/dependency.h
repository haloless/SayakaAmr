
#pragma once

#ifndef INT64
 typedef __int64 INT64;
#endif
#ifndef INT32
 typedef __int32 INT32;
#endif

#if defined(WIN32) || defined(WIN64)

#include <direct.h>	//mkdir

//increments value, returns OLD value
template<typename T>
inline T AtomicIncrement(volatile T* ptr)
{
	return _InterlockedIncrement(ptr)-1;
}

//increments value, returns OLD value
template<>
inline INT64 AtomicIncrement<INT64>(volatile INT64* ptr)
{
	return _InterlockedIncrement64(ptr)-1;
}

//exchanges value, returns OLD value
template<typename T,typename T2>
inline T AtomicExchange(volatile T* ptr,T2 newValue)
{
	return _InterlockedExchange(ptr,newValue);
}

inline long AtomicAdd(volatile long *ptr,long value)
{
	return _InterlockedExchangeAdd(ptr,value);
}

inline __int64 AtomicAdd(volatile __int64 *ptr,__int64 value)
{
	return _InterlockedExchangeAdd64(ptr,value);
}

inline double AtomicAdd(volatile double *target,double d)
{
	double prev;
	double tmp;
	prev=*target;
	tmp=prev+d;
	__int64 result;
	while((result=_InterlockedCompareExchange64((__int64*)target,*((__int64*)&tmp),*((__int64*)&prev))) != (*(__int64*)&prev)){
		prev=*((double*)&result);
		tmp=prev+d;
	}
	return *((double*)&result);
}

inline float AtomicAdd(volatile float *target,float d)
{
	float prev;
	float tmp;
	prev=*target;
	tmp=prev+d;
	long result;
	while((result=_InterlockedCompareExchange((long*)target,*((long*)&tmp),*((long*)&prev))) != (*(long*)&prev)){
		prev=*((float*)&result);
		tmp=prev+d;
	}
	return *((float*)&result);
}




//exchanges value, returns OLD value
template<>
inline INT64 AtomicExchange<INT64,INT64>(volatile INT64* ptr,INT64 newValue)
{
	return _InterlockedExchange64(ptr,newValue);
}

inline void create_dir(const char* dirname)
{
	_mkdir(dirname);
}


inline double AtomicAbsMin(volatile double *target,double d)
{
	double prev;
	double tmp;
	prev=*target;
	if(fabs(prev)>fabs(d)){
		tmp=d;
		__int64 result;
		while((result=_InterlockedCompareExchange64((__int64*)target,*((__int64*)&tmp),*((__int64*)&prev))) != (*(__int64*)&prev)){
			prev=*((double*)&result);
			if(fabs(prev)>fabs(d)){
				tmp=d;
			}else{
				return prev;
			}
		}
		return *((double*)&result);
	}else{
		return prev;
	}
}

inline double AtomicAbsMax(volatile double *target,double d)
{
	double prev;
	double tmp;
	prev=*target;
	if(fabs(prev)<fabs(d)){
		tmp=d;
		__int64 result;
		while((result=_InterlockedCompareExchange64((__int64*)target,*((__int64*)&tmp),*((__int64*)&prev))) != (*(__int64*)&prev)){
			prev=*((double*)&result);
			if(fabs(prev)<fabs(d)){
				tmp=d;
			}else{
				return prev;
			}
		}
		return *((double*)&result);
	}else{
		return prev;
	}
}

#else

#define AtomicIncrement(x)		(__sync_fetch_and_add(x,1))

#define AtomicExchange(ptr,newValue) (__sync_lock_test_and_set(ptr,newValue))

template <typename T,typename T2>
inline T AtomicAdd(volatile T *ptr,T2 value)
{
	return __sync_fetch_and_add(ptr,value);
}
template<>
inline double AtomicAdd<double,double>(volatile double *target,double d)
{
	double prev;
	double tmp;
	prev=*target;
	tmp=prev+d;
	__int64 result;
	while((result=__sync_val_compare_and_swap((__int64*)target,*((__int64*)&tmp),*((__int64*)&prev))) != (*(__int64*)&prev)){
		prev=*((double*)&result);
		tmp=prev+d;
	}
	return *((double*)&result);
}

template<>
inline float AtomicAdd<float,float>(volatile float *target,float d)
{
	float prev;
	float tmp;
	prev=*target;
	tmp=prev+d;
	__int32 result;
	while((result=__sync_val_compare_and_swap((__int32*)target,*((__int32*)&tmp),*((__int32*)&prev))) != (*(__int32*)&prev)){
		prev=*((float*)&result);
		tmp=prev+d;
	}
	return *((float*)&result);
}

inline void create_dir(const char* dirname)
{
	mkdir(dirname,S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
}

#endif
