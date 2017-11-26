#include "BezierSpline.h"

using namespace std;

template<class V> bezier<V>::bezier() {}

template<class V> bezier<V>::bezier(V *points, int numPoints) {
	control.assign(points, points + numPoints);
}

template<class V> bezier<V>::bezier(const vector<V>& points) {
	control = points;
}

template<class V> V bezier<V>::getPoint(float u) {
	float u_n[control.size()];		//u^n
	float u_n2[control.size()];		//(1-u)^n

	u_n[0] = 1;
	u_n2[0] = 1;

	for (int i = 1; i < control.size(); i++) {
		u_n[i] = u_n[i - 1] * u;
		u_n2[i] = u_n2[i - 1] * (u - 1);
	}

	V point = V();		//Works providing contstructor defaults to 0
	int n = control.size();

	for (int i = 0; i < control.size(); i++) {
		point += u_n2[n - 1 - i] * u_n[i]*control[i];
	}

	return point;
}

template<class V> V bezier<V>::getQuadPoint(float u) {
	if (control.size(); >= 3)
		return (1 - u)*(1 - u)*control[0] + (1 - u)*u*control[1] + u*u*control[2];
	else
		return V();
}

template<class V> std::vector<V> getQuadPoints(int divisions, float uStart, float uEnd) {
	float uStep = (uStart - uEnd) / float(divisions - 1);
	float u = 0;

	vector<V> points;

	for (int i = 0; i < divisions; i++) {
		points.push_back(getQuadPoint(u));
	}

	return points;
}




