#pragma once

#include <vector>

namespace sayaka
{

template<typename T>
class Array : public std::vector<T>
{
	typedef std::vector<T> super_type;

public:
	//
	Array() 
		: super_type() 
	{}
	//
	explicit Array(size_t len) 
		: super_type(len) 
	{}
	//
	Array(size_t len, const T &val)
		: super_type(len, val)
	{}

	
	//
	T* data_ptr() { return &this->operator[](0); }
	const T* data_ptr() const { return &this->operator[](0); }
};

} // namespace_sayaka



