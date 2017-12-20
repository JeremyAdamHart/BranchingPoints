#include "Skeleton.h"

using namespace glm;

Link::Link(Joint *a, Joint *b) :a(a), b(b) {}

vec3 Link::dir() const {
	return b->pos - a->pos;
}

Joint::Joint(glm::vec3 pos) :pos(pos) {}

void Joint::addLink(Joint *b) {
	links.push_back(Link(this, b));
}

float getAngle(vec3 v1, vec3 v2, vec3 n) {
	return atan2(dot(n, cross(v1, v2)), dot(v1, v2));
}

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

//TODO: Reconsider what space angles are passed
vec3 Skeleton::getEndpoint(int link, float theta) {
//	float newTheta = convertAngle(0, link, theta);
	float newTheta = theta;
	vec3 bx = bases[link].bx;
	vec3 by = (link == 0) ? bases[link].by1 : bases[link].by2;
	return joint->links[link].b->pos + (bx*cos(newTheta) + by*sin(newTheta))*radius;
}

vec3 Skeleton::getDir(int link, float theta) {
//	float newTheta = convertAngle(0, link, theta);
	float newTheta = theta;
	vec3 bx = bases[link].bx;
	vec3 by = (link == 0) ? bases[link].by1 : bases[link].by2;

	return normalize(bx*cos(newTheta) + by*sin(newTheta));
}

vec3 Skeleton::getDir(int link) {
	return normalize(joint->links[link].dir());
}

Skeleton::Skeleton(Joint *joint, float radius) :joint(joint), radius(radius) {
	vec3 center1 = joint->links[0].b->pos;
	vec3 lA = normalize(joint->links[0].dir());
	int n = joint->links.size();
	bases.push_back(getFrame(lA, normalize(joint->links[1].dir())));

	for (int i = 1; i<joint->links.size(); i++) {
		vec3 lB = normalize(joint->links[i].dir());
		bases.push_back(getFrame(lA, lB));
	}
}

std::vector<bezier<vec4>> Skeleton::getCurveSet(int link, float theta) {
	vector<bezier<vec4>> curveSet;
//	float linkTheta = convertAngle(0, link, theta);
	float linkTheta = theta;

	vec3 center1 = joint->links[link].b->pos;
	vec3 lA = normalize(joint->links[link].dir());
	int n = joint->links.size();
	vec3 p1 = getEndpoint(link, linkTheta);

	for (int i = (link + 1) % n; i != link; i = (i + 1) % n) {
		float newTheta = convertAngle(link, i, linkTheta);
		vec3 lB = normalize(joint->links[i].dir());
		vec3 p3 = getEndpoint(i, newTheta);
		vec3 p2 = lineIntersection(p1, lA, p3, lB);

		//Get angle for w
		vec3 perp = normalize(cross(lA, lB));
		vec3 dir = normalize(p3 - joint->links[i].b->pos);
		vec3 perp2 = normalize(cross(perp, lB));

//		float w = 1.f / max(dot(dir, perp2), 0.0001f);
		float w = 1.f / max(dot(getDir(link, theta), getDir(i)), 0.0001f);
		curveSet.push_back(bezier<vec4>({ vec4(p1, 1), vec4(p2, 1)*w, vec4(p3, 1) }));
	}

	return curveSet;
}

float Skeleton::convertAngle(int linkA, int linkB, float theta) {
	if (linkA == linkB)
		return theta;

	vec3 normal = normalize(joint->links[0].dir());
	if (linkA == 0 || linkB == 0) {
		float thetaDiff = getAngle(bases[linkA].bx, bases[linkB].bx, normal);
		return theta + thetaDiff;
	}
	else {
		vec3 v1 = normalize(joint->links[linkA].dir());
		vec3 v2 = normalize(joint->links[linkB].dir());
		BasisPair basis = getFrame(v1, v2);

		float thetaDiff1 = getAngle(basis.bx, bases[linkA].bx, v1);
		float thetaDiff2 = getAngle(basis.bx, bases[linkB].bx, v2);
		return  -theta - thetaDiff1 - thetaDiff2;
	}
}