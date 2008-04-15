#include "Vector.h"

using namespace std;

float Vector::distanceToSegment(const Vector& p1, const Vector& p2) const
{
    // vector from p1 to p2
    Vector p1p2(p1,p2);

    float length = p1p2.mag();
    if (length <= 1.0e-9)
	return Vector::distance(*this, p1);

    p1p2.normalize();

    // vector from p1 to this
    Vector p1point(p1, *this);
	
    // u is distance along segment from point 1
    float u = Vector::dot(p1point, p1p2);

    if (u <= 0.0)  // closest to point 1
	return p1point.mag();
    else if (u >= length) // closest to point 2
	return Vector::distance(*this, p2);
    else // use pythagorean theorem to compute projected distance
	return sqrt(p1point.mag_squared() - u*u);

    return 0.0;
}

ostream& operator<< (ostream& s, const Vector& v)
{
    s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return s;
}
