#include <utility>

#include "BSplineSkinner.h"
#include "BezierSpline.h"

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

	return ret;
}

float getUFromS2(float s, float w) {
	float t = 1 - s;
	return (-sqrt(t*t*w*w - t*t + t) + t*w - t + 1) /
		(2 * t*w - 2 * t + 1);
}

//Performs line intersection by looking for closest point to deal with floating point inaccuracy
template<class V> V lineIntersection(V p1, V l1, V p2, V l2) {
	float s = (dot(l1, p2) - dot(l1, p1) - dot(l1, l2)*(dot(l2, p2) - dot(l2, p1))) /
		(1 - dot(l1, l2)*dot(l1, l2));
	return p1 + l1*s;
}

vector<bezier<vec4>> getCurveSet(Joint *joint, int link, float theta) {
	vec3 p = joint->links[link].b->pos;
	vec3 l1 = normalize(joint->links[link].dir());
	int n = joint->links.size();
	BasisPair initialBasis = getFrame(l1, normalize(joint->links[(link + 1) % n].dir()));

	for (int i = (link + 1) % n; i != link; i = (i + 1) % n) {
		vec3 p3 = joint->links[i].b->pos;
		vec3 l3 = normalize(joint->links[link].dir());
	//	vec3 p2 = lineIntersection()
	}

	return vector<bezier<vec4>>();
}

vec3 getCirclePoint(vec3 center, vec3 bx, vec3 by, float radius, float theta) {
	return center + (bx*cos(theta) + by*sin(theta))*radius;
}

vec3 generatePoint(Joint *joint, int link, float s, float theta) {


	return vec3();
}

void BSplineSkinner::generateCurve(Joint *joint, int link) {
	/*vec3 v0 = joint->links[link].dir();
	vec3 c0 = joint->links[link].b->pos;


	vector<vec3> v1;
	vector<vec3> m;

	for (int i = 0; i < joint->links.size(); i++) {
		if (i != link) {
			v1 = joint->links[i].dir();
			vec3 c1 = joint->links[i].b->pos;
			BasisPair bp = getFrame(v0, v1);
		}
	}*/
}