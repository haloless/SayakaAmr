#pragma once

#include <math.h>
#include <float.h>

#include <utility>
#include <iostream>

#ifdef USE_ATOMIC_VEC_OP
#include "dependency.h"
#endif

inline double pow2(double x) { return x*x; }
inline double pow3(double x) { return x*x*x; }

inline double fsign(double x){return x>0 ? 1.0 : -1.0;}

inline double fabsmax(double a,double b) {
	if(fabs(a)>fabs(b)){
		return a;
	}else{
		return b;
	}
}
inline double fabsmin(double a,double b) {
	if(fabs(a)<fabs(b)){
		return a;
	}else{
		return b;
	}
}

template<typename T>
struct Vector3t
{
	T x;
	T y;
	T z;

public:
	const T* data() const { return &x; }
	T* data() { return &x; }

	const T& operator()(int dir) const { return data()[dir]; }
	T& operator()(int dir) { return data()[dir]; }

	//
	void setZero() {
		this->x = (T) 0;
		this->y = (T) 0;
		this->z = (T) 0;
	}

	//
	Vector3t operator-() const {
		Vector3t neg;
		neg.x = -this->x;
		neg.y = -this->y;
		neg.z = -this->z;
		return neg;
	}

	//
	Vector3t& operator+=(const Vector3t &vec) {
		this->x += vec.x;
		this->y += vec.y;
		this->z += vec.z;
		return *this;
	}
	Vector3t& operator-=(const Vector3t &vec) {
		this->x -= vec.x;
		this->y -= vec.y;
		this->z -= vec.z;
		return *this;
	}
	Vector3t& operator*=(const T& s) {
		this->x *= s;
		this->y *= s;
		this->z *= s;
		return *this;
	}
	Vector3t& operator/=(const T& s) {
		this->x /= s;
		this->y /= s;
		this->z /= s;
		return *this;
	}

	// compare by all components
	bool allLT(const Vector3t &rhs) const { return x < rhs.x && y < rhs.y && z < rhs.z; }
	bool allLE(const Vector3t &rhs) const { return x <= rhs.x && y <= rhs.y && z <= rhs.z; }
	bool allGT(const Vector3t &rhs) const { return x > rhs.x && y > rhs.y && z > rhs.z; }
	bool allGE(const Vector3t &rhs) const { return x >= rhs.x && y >= rhs.y && z >= rhs.z; }

	Vector3t minVec(const Vector3t &rhs) const {
		return Vector3t{ std::min(x, rhs.x), std::min(y, rhs.y), std::min(z,rhs.z) };
	}

	Vector3t maxVec(const Vector3t &rhs) const {
		return Vector3t{ std::max(x, rhs.x), std::max(y, rhs.y), std::max(z,rhs.z) };
	}

	auto asTuple() const {
		return std::make_tuple(x, y, z);
	}
public:

	//
	static inline Vector3t VecMake(const T &vx, const T &vy, const T &vz) {
		Vector3t v;
		v.x = vx; v.y = vy; v.z = vz;
		return v;
	}
	static inline Vector3t VecZero() {
		Vector3t vzero;
		vzero.setZero();
		return vzero;
	}
	static inline Vector3t VecUnit() {
		Vector3t v;
		v.x = (T) 1;
		v.y = (T) 1;
		v.z = (T) 1;
		return v;
	}

	//static inline Vector3t VecMin(const Vector3t &a, const Vector3t &b) {
	//	return Vector3t{ std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) };
	//}
	//static inline Vector3t VecMax(const Vector3t &a, const Vector3t &b) {
	//	return Vector3t{ std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) };
	//}
};
//
template<typename T>
inline Vector3t<T> operator+(const Vector3t<T> &a, const Vector3t<T> &b) {
	Vector3t<T> c = a; 
	c += b;
	return c;
}
template<typename T>
inline Vector3t<T> operator-(const Vector3t<T> &a, const Vector3t<T> &b) {
	Vector3t<T> c = a;
	c -= b;
	return c;
}
template<typename T>
inline Vector3t<T> operator*(const Vector3t<T> &v, const T &s) {
	Vector3t<T> vs = v;
	vs *= s;
	return vs;
}
template<typename T>
inline Vector3t<T> operator*(const T &s, const Vector3t<T> &v) {
	Vector3t<T> vs = v;
	vs *= s;
	return vs;
}
template<typename T>
inline bool operator==(const Vector3t<T> &a, const Vector3t<T> &b) {
	return a.x==b.x && a.y==b.y && a.z==b.z;
}
template<typename T>
inline bool operator!=(const Vector3t<T> &a, const Vector3t<T> &b) {
	return !(a==b);
}

// IO
template<typename T>
std::ostream& operator<<(std::ostream &os, const Vector3t<T> &v) {
	os << '[' << v.x << ',' << v.y << ',' << v.z << ']';
	return os;
}

template<typename T>
std::istream& operator>>(std::istream &is, Vector3t<T> &v) {
	char c;
	is >> c >> c >> v.x >> c >> v.y >> c >> v.z >> c;
	return is;
}

//int
/*
struct Vector3i
{
	int ix;
	int iy;
	int iz;

public:
	inline const int* data() const { return &ix; }
	inline int* data() { return &ix; }

	inline int operator()(int dir) const { return data()[dir]; }
	int& operator()(int dir) { return data()[dir]; }
};

inline void Vec3ZeroI(Vector3i *pV)
{
	pV->ix=0;
	pV->iy=0;
	pV->iz=0;
}
*/
typedef Vector3t<int> Vector3i;

//double
/*
struct Vector3d{
	double x;
	double y;
	double z;

public:
	const double* data() const { return &x; }
	double* data() { return &x; }

	double operator()(int dir) const { return data()[dir]; }
	double& operator()(int dir) { return data()[dir]; }
};
*/
typedef Vector3t<double> Vector3d;

//vOut = a + b
inline void Vec3Add(Vector3d* vOut,const Vector3d* a,const Vector3d* b)
{
	double x=a->x+b->x;
	double y=a->y+b->y;
	double z=a->z+b->z;
	vOut->x=x;
	vOut->y=y;
	vOut->z=z;
}

#ifdef USE_ATOMIC_VEC_OP
//vOut = a + b; for parallel computing
inline void AtomicVec3Add(volatile Vector3d* vOut,const Vector3d* a)
{
	double x=a->x;
	double y=a->y;
	double z=a->z;

	AtomicAdd(&vOut->x,x);
	AtomicAdd(&vOut->y,y);
	AtomicAdd(&vOut->z,z);
}
#endif

//vOut = a - b
inline void Vec3Sub(Vector3d* vOut,const Vector3d* a,const Vector3d* b)
{
	double x=a->x-b->x;
	double y=a->y-b->y;
	double z=a->z-b->z;
	vOut->x=x;
	vOut->y=y;
	vOut->z=z;
}

#ifdef USE_ATOMIC_VEC_OP
//vOut = a - b; for parallel computing
inline void AtomicVec3Sub(volatile Vector3d* vOut,const Vector3d* a)
{
	double x=a->x;
	double y=a->y;
	double z=a->z;
	AtomicAdd(&vOut->x,-x);
	AtomicAdd(&vOut->y,-y);
	AtomicAdd(&vOut->z,-z);
}
#endif

//returns dot(a, b)
inline double Vec3InnerProd(const Vector3d* a,const Vector3d* b)
{
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

//cross product; vOut = cross(a, b)
inline void Vec3OuterProd(Vector3d* vOut,const Vector3d* a,const Vector3d* b)
{
	double x=(a->y * b->z)-(a->z * b->y);
	double y=(a->z * b->x)-(a->x * b->z);
	double z=(a->x * b->y)-(a->y * b->x);

	vOut->x=x;
	vOut->y=y;
	vOut->z=z;
}

//length^2
inline double Vec3LengthSq(const Vector3d* a)
{
	return a->x*a->x + a->y*a->y + a->z*a->z;
}

//length^2
inline double Vec3LengthSq(const Vector3d* a,const Vector3d* b)
{
	return pow2(a->x-b->x) + pow2(a->y-b->y) + pow2(a->z-b->z);
}

//length
inline double Vec3Length(const Vector3d* a)
{
	return sqrt(a->x*a->x + a->y*a->y + a->z*a->z);
}
inline double Vec3Length(const Vector3d &a) {
	return Vec3Length(&a);
}

//length
inline double Vec3Length(const Vector3d* a,const Vector3d* b)
{
	return sqrt(pow2(a->x-b->x) + pow2(a->y-b->y) + pow2(a->z-b->z));
}
inline double Vec3Length(const Vector3d& a,const Vector3d& b) {
	return sqrt(pow2(a.x-b.x) + pow2(a.y-b.y) + pow2(a.z-b.z));
}


//normalize
inline void Vec3Normalize(Vector3d* a)
{
	double fLen=Vec3Length(a);
	double x=a->x/fLen;
	double y=a->y/fLen;
	double z=a->z/fLen;
	a->x=x;
	a->y=y;
	a->z=z;
}

//normalize if the vector is not zero; zero vector is left untouched
inline void Vec3NormalizeNonZero(Vector3d* a)
{
	double fLen=Vec3Length(a);
	if(fLen>FLT_MIN){
		double x=a->x/fLen;
		double y=a->y/fLen;
		double z=a->z/fLen;
		a->x=x;
		a->y=y;
		a->z=z;
	}
}

inline void Vec3Invert(Vector3d* a)
{
	double x=-a->x;
	double y=-a->y;
	double z=-a->z;
	a->x=x;
	a->y=y;
	a->z=z;
}

//vOut = a * f
inline void Vec3Mul(Vector3d* vOut,const Vector3d* a,double f)
{
	double x=a->x*f;
	double y=a->y*f;
	double z=a->z*f;
	vOut->x=x;
	vOut->y=y;
	vOut->z=z;
}

//vOut = vOut + a * f
inline void Vec3MulAdd(Vector3d* vOut,const Vector3d* a,double f)
{
	double x=a->x*f;
	double y=a->y*f;
	double z=a->z*f;
	vOut->x+=x;
	vOut->y+=y;
	vOut->z+=z;
}

//vOut = a / f
inline void Vec3Div(Vector3d* vOut,const Vector3d* a,double f)
{
	double x=a->x/f;
	double y=a->y/f;
	double z=a->z/f;
	vOut->x=x;
	vOut->y=y;
	vOut->z=z;
}

//vOut = {0}
inline void Vec3Zero(Vector3d* vOut)
{
	vOut->x=0.0;
	vOut->y=0.0;
	vOut->z=0.0;
}

