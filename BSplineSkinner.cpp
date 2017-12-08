#include <utility>

#include "BSplineSkinner.h"

using namespace glm;
using namespace std;

//Returns angle between vectors v1 and v2, sign corresponding to n
//	normal should be normalized
float getAngle(vec3 v1, vec3 v2, vec3 n) {
	return atan2(dot(n, cross(v1, v2)), dot(v1, v2));
}

struct BasisPair {
	vec3 bx, by1, by2;
};

BasisPair getFrame(vec3 v1, vec3 v2) {
	BasisPair ret;
	ret.bx = cross(v1, v2);
	ret.by1 = cross(-v1, ret.bx);
	ret.by2 = cross(v2, ret.bx);
}

//Performs line intersection by looking for closest point to deal with floating point inaccuracy
template<V> V lineIntersection(V p1, V l1, V p2, V l2) {
	float s =	(dot(l1, p2) - dot(l1, p1) - dot(l1, l2)*(dot(l2, p2) - dot(l2, p1)) /
				(1 - dot(l1, l2)*dot(l1, l2));
	return p1 + l1*s;
}

void BSplineSkinner::generateCurve(Joint *joint, int link) {
	
}