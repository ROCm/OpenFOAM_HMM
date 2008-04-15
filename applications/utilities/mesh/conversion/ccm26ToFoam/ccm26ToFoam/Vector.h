#ifndef _CCMVECTOR_H
#define _CCMVECTOR_H

#include <math.h>
#include <iostream>

class Vector
{
 public:
    Vector() {
	_xyz[0] = 0.0;
	_xyz[1] = 0.0;
	_xyz[2] = 0.0;
    }
    Vector(float x, float y, float z) {
	_xyz[0] = x;
	_xyz[1] = y;
	_xyz[2] = z;
    }
    Vector(const Vector& p1, const Vector& p2) {	    
	_xyz[0] = p2._xyz[0] - p1._xyz[0];
	_xyz[1] = p2._xyz[1] - p1._xyz[1];
	_xyz[2] = p2._xyz[2] - p1._xyz[2];
    }
    Vector(const float *xyz) {
	_xyz[0] = xyz[0];
	_xyz[1] = xyz[1];
	_xyz[2] = xyz[2];
    }

    void set(int i, float val) {
	_xyz[i] = val;
    }

    void setXYZ(const float *xyz) {
	_xyz[0] = xyz[0];
	_xyz[1] = xyz[1];
	_xyz[2] = xyz[2];
    }

    void setXYZ(float x, float y, float z) {
	_xyz[0] = x;
	_xyz[1] = y;
	_xyz[2] = z;
    }

    void axpy(float a, const Vector& x, const Vector& y) {
	_xyz[0] = a*x._xyz[0] + y._xyz[0];
	_xyz[1] = a*x._xyz[1] + y._xyz[1];
	_xyz[2] = a*x._xyz[2] + y._xyz[2];
    }


    float x() const {return _xyz[0];}
    float y() const {return _xyz[1];}
    float z() const {return _xyz[2];}
    const float *xyz() const {return _xyz;}
    float& operator[](unsigned int i) { return _xyz[i]; }

    float mag() const {
	return sqrt(_xyz[0]*_xyz[0] + 
		    _xyz[1]*_xyz[1] + 
		    _xyz[2]*_xyz[2]);
    }

    float mag_squared() const {
	return _xyz[0]*_xyz[0] + 
	    _xyz[1]*_xyz[1] + 
	    _xyz[2]*_xyz[2];
    }

    int normalize() {
	float mag = this->mag();
	if (mag < 1.0e-9)
	    return 1;

	float scale = 1.0f / mag;
	for(int i=0; i<3; i++)
	    _xyz[i] *= scale;
	return 0;
    }

    Vector& operator+=(const Vector& v1) {
	_xyz[0] += v1._xyz[0];
	_xyz[1] += v1._xyz[1];
	_xyz[2] += v1._xyz[2];
	return *this;
    }

    Vector& operator-=(const Vector& v1) {
	_xyz[0] -= v1._xyz[0];
	_xyz[1] -= v1._xyz[1];
	_xyz[2] -= v1._xyz[2];
	return *this;
    }

    Vector& operator*=(const float scale) {
	_xyz[0] *= scale;
	_xyz[1] *= scale;
	_xyz[2] *= scale;
	return *this;
    }

    Vector operator-(const Vector &v) const {
	return Vector(
	    this->_xyz[0] - v._xyz[0],
	    this->_xyz[1] - v._xyz[1],
	    this->_xyz[2] - v._xyz[2]);
    }

    Vector operator+(const Vector &v) const {
	return Vector(
	    this->_xyz[0] + v._xyz[0],
	    this->_xyz[1] + v._xyz[1],
	    this->_xyz[2] + v._xyz[2]);
    }

    static float distance(const Vector& v1, const Vector& v2) {
	return sqrt((v1._xyz[0]-v2._xyz[0])*(v1._xyz[0]-v2._xyz[0])+
		    (v1._xyz[1]-v2._xyz[1])*(v1._xyz[1]-v2._xyz[1])+
		    (v1._xyz[2]-v2._xyz[2])*(v1._xyz[2]-v2._xyz[2]));
    }

    static float distance_squared(const Vector& v1, const Vector& v2) {
	return (v1._xyz[0]-v2._xyz[0])*(v1._xyz[0]-v2._xyz[0])+
	    (v1._xyz[1]-v2._xyz[1])*(v1._xyz[1]-v2._xyz[1])+
	    (v1._xyz[2]-v2._xyz[2])*(v1._xyz[2]-v2._xyz[2]);
    }

    static Vector cross(const Vector& v1, const Vector& v2) {
	return Vector(
	    v1._xyz[1] * v2._xyz[2] - v1._xyz[2] * v2._xyz[1],
	    v1._xyz[2] * v2._xyz[0] - v1._xyz[0] * v2._xyz[2],
	    v1._xyz[0] * v2._xyz[1] - v1._xyz[1] * v2._xyz[0]);
    }

    static Vector subtract_component(const Vector& v, const Vector& dir) {
	float dot = Vector::dot(v, dir);
	return Vector(v._xyz[0] - dot*dir._xyz[0],
		      v._xyz[1] - dot*dir._xyz[1],
		      v._xyz[2] - dot*dir._xyz[2]);
    }
		
    static float dot(const Vector& v1, const Vector& v2) {
	return v1._xyz[0]*v2._xyz[0] +
	    v1._xyz[1]*v2._xyz[1] +
	    v1._xyz[2]*v2._xyz[2];
    }

    float distanceToSegment(const Vector& p1, const Vector& p2) const;
 private:
    float _xyz[3];
};

std::ostream& operator<< (std::ostream& s, const Vector& v);

#endif

