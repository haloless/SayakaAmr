#pragma once

#include <vector>

#include "SayakaCommons.h"


SAYAKA_NS_BEGIN;

template<typename T>
using Array = std::vector<T>;

//
//template<typename T>
//class Array : public std::vector<T>
//{
//	using super_type = std::vector<T>;
//
//public:
//
//	using super_type::super_type;
//
//	////
//	//Array() 
//	//	: super_type() 
//	//{}
//	////
//	//explicit Array(size_t len) 
//	//	: super_type(len) 
//	//{}
//	////
//	//Array(size_t len, const T &val)
//	//	: super_type(len, val)
//	//{}
//
//	
//	//
//	T* data_ptr() { return &this->operator[](0); }
//	const T* data_ptr() const { return &this->operator[](0); }
//};

template<typename T>
class PtrArray : public Array<T*>
{
	using super_type = Array<T*>;
	using value_type = T;
	using ptr_type = T * ;

protected:
	bool m_managed;

public:
	PtrArray(bool managed = true) 
		: super_type(), m_managed(managed) {}

	PtrArray(const PtrArray &rhs) = delete;

	//PtrArray(PtrArray &&rhs) {
	//	this->swap(rhs);

	//}

	~PtrArray() {
		if (m_managed) {
			for (auto it = begin(); it != end(); ++it) {
				delete *it;
				*it = nullptr;
			}
		}
	}

	bool isManaged() const { return m_managed; }



};



SAYAKA_NS_END;


