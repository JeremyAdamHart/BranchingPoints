#include <utility>

#include "BSplineSkinner.h"
#include "BezierSpline.h"

#define M_PI 3.1415926535897932384626433832795f


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
	ret.bx = normalize(cross(v1, v2));
	ret.by1 = normalize(cross(-v1, ret.bx));
	ret.by2 = normalize(cross(v2, ret.bx));

	return ret;
}

vec3 getCirclePoint(vec3 center, vec3 bx, vec3 by, float radius, float theta) {
	return center + (bx*cos(theta) + by*sin(theta))*radius;
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

vector<bezier<vec4>> getCurveSet(Joint *joint, int link, float radius, float theta) {
	vector<bezier<vec4>> curveSet;
	
	vec3 center1 = joint->links[link].b->pos;
	vec3 lA = normalize(joint->links[link].dir());
	int n = joint->links.size();
	BasisPair initBasis = getFrame(lA, normalize(joint->links[(link + 1) % n].dir()));
	vec3 p1 = getCirclePoint(center1, initBasis.bx, initBasis.by1, radius, theta);

	for (int i = (link + 1) % n; i != link; i = (i + 1) % n) {
		vec3 center2 = joint->links[i].b->pos;
		vec3 lB = normalize(joint->links[i].dir());
		BasisPair basis = getFrame(lA, lB);
		float thetaDiff = getAngle(initBasis.bx, basis.bx, lA);
		
		vec3 p3 = getCirclePoint(center2, basis.bx, basis.by2, radius, theta + thetaDiff);
		vec3 p2 = lineIntersection(p1, lA, p3, lB);
		float w = 1.f / max(cos(theta + thetaDiff - M_PI*0.5f), 0.0001f);
		curveSet.push_back(bezier<vec4>({ vec4(p1, 1), vec4(p2, 1)*w, vec4(p3, 1) }));
	}

	return curveSet;
}

vec3 generatePoint(Joint *joint, int link, float s, float radius, float theta) {
	vector<bezier<vec4>> curveSet = getCurveSet(joint, link, radius, theta);
	vec3 p1 = curveSet[0].control[0];
	vec3 pEnd = p1 - joint->links[link].dir();
	float baseLength = length(p1 - pEnd);
	vec3 sharedPoint = p1 + s*(pEnd - p1);
	vec3 blendedPoint = sharedPoint;
	float valid = 1.f;
	float totalWeight = 0.f;
	vec3 offsets(0.f, 0.f, 0.f);

	//Test
	vec3 currentDir = normalize(p1 - joint->links[link].b->pos);

	for (int i = 0; i < curveSet.size(); i++) {
		vec4 p2 = curveSet[i].control[1];
		float sRatio = baseLength / length(p1 - vec3(p2)/p2.w);
		float u = getUFromS2(s*sRatio, p2.w);
		if (u > 0.99f)
			valid = 0.f;
		vec4 curvePoint = curveSet[i].getQuadPoint(u);
		blendedPoint += vec3(curvePoint)/curvePoint.w - sharedPoint;
		vec3 offset = vec3(curvePoint) / curvePoint.w - sharedPoint;

		//Test
		vec3 dir = normalize(curveSet[i].control[2] - p2);
		float angle = acos(dot(dir, currentDir));
		float weight = max(dot(dir, currentDir), 0.f);
		offset *= weight;

		totalWeight += weight; //length(offset);
		offsets += offset;	// *length(offset);
	}

	if (totalWeight < 0.001f)
		return sharedPoint;
	return valid*(sharedPoint + offsets / totalWeight);

	return blendedPoint*valid;
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