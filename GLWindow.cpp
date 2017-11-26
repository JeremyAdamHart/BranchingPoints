#include "GLWindow.h"

#include <iostream>

using namespace glm;
using namespace std;

#include "Drawable.h"
#include "SimpleGeometry.h"
#include "SimpleShader.h"
#include "ColorMat.h"
#include "TrackballCamera.h"
#include "BezierSpline.h"
#include <sstream>

#include <glm/gtc/matrix_transform.hpp>

using namespace renderlib;

#define M_PI 3.1415926535897932384626433832795f

TrackballCamera cam(
	vec3(0, 0, -1), vec3(0, 0, 1),
	glm::perspective(80.f*M_PI/180.f, 1.f, 0.1f, 3.f));

bool reloadShaders = false;
bool windowResized = false;
int windowWidth, windowHeight;

void cursorPositionCallback(GLFWwindow *window, double xpos, double ypos) {
	static vec2 lastPos = vec2(0.f, 0.f);
	
	int vp[4];
	glGetIntegerv(GL_VIEWPORT, vp);
	vec2 mousePos = vec2(float(xpos) / float(vp[2]), 
		float(-ypos) / float(vp[3]))*2.f
		- vec2(1.f, 1.f);

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		vec2 diff = mousePos - lastPos;
		cam.trackballRight(-diff.x*3.14159f);
		cam.trackballUp(-diff.y*3.14159f);
	}
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS) {
		vec2 diff = mousePos - lastPos;
		cam.zoom(pow(1.5, -diff.y));
	}

	lastPos = mousePos;
}

void windowResizeCallback(GLFWwindow* window, int width, int height) {
	windowResized = true;
	windowWidth = width;
	windowHeight = height;
}

void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_SPACE && action == GLFW_RELEASE)
		reloadShaders = true;
}

WindowManager::WindowManager() :
window_width(800), window_height(800)
{
	if(!glfwInit()){
		printf("GLFW failed to initialize\n");
	}
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);	
	window = createWindow(window_width, window_height, 
		"Perlin Noise");

	glfwMakeContextCurrent(window);
	initGLExtensions();

	glfwSwapInterval(1);
	glfwSetCursorPosCallback(window, cursorPositionCallback);

	glClearColor(1.f, 1.f, 1.f, 1.f);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glViewport(0, 0, window_width, window_height);
}

WindowManager::WindowManager(int width, int height, std::string name, glm::vec4 color) :
	window_width(width), window_height(height) 
{
	if(!glfwInit()){
		printf("GLFW failed to initialize\n");
	}
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);	
	window = createWindow(window_width, window_height, name);

	glfwMakeContextCurrent(window);
	initGLExtensions();


	glfwSetCursorPosCallback(window, cursorPositionCallback);

	glClearColor(color.r, color.g, color.b, color.a);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glViewport(0, 0, window_width, window_height);
}

void generateCurves(vector<Drawable>& drawables, int numCurves, int numDivisions, float radius) {
	
	float thetaStep = M_PI / float(numCurves - 1);
	float theta = 0;

	vector<bezier<vec4>> curves;

	for (int i = 0; i < numCurves; i++) {
		vec4 circleA = vec4(-cos(theta)*radius, 0.0, sin(theta)*radius, 1) + vec4(0.5, -0.5, 0, 0);
		vec4 circleB = vec4(0, -cos(theta)*radius, sin(theta)*radius, 1) + vec4(-0.5, 0.5, 0, 0);
		float w = 1.f / max(cos(theta), 0.0001f);
		vec4 middle = vec4(circleA.x, circleB.y, circleA.z, 1.0)*w;

		curves.push_back(bezier<vec4>({ circleA, middle, circleB }));

		theta += thetaStep;
	}

	for (int i = 0; i < curves.size(); i++) {
		vector<vec4> points4D = curves[i].getQuadPoints(50);
		vector<vec3> points;
		for (int j = 0; j < points4D.size(); j++) {
			const vec4& p = points4D[j];
			points.push_back(vec3(p.x / p.w, p.y / p.w, p.z / p.w));
		}

		vector<vec3> control;
		for (int j = 0; j < curves[i].control.size(); j++) {
			const vec4& p = curves[i].control[j];
			control.push_back(vec3(p.x / p.w, p.y / p.w, p.z / p.w) + vec3(0.003, 0.003, 0.003));
		}
		drawables.push_back(Drawable(
			new ColorMat(vec3(0.f, 0.f, 1.f)),
			new SimpleGeometry(control.data(), control.size(), GL_LINE_STRIP)));
		drawables.push_back(Drawable(
			new ColorMat(vec3(1.f, 0.f, 0.f)),
			new SimpleGeometry(points.data(), points.size(), GL_LINE_STRIP)));
	}
}

vec4 circleX(float theta, vec3 center, float radius) {
	return vec4(center, 1) + vec4(0, sin(theta), -cos(theta), 0)*radius;
}

vec4 circleY(float theta, vec3 center, float radius) {
	return vec4(center, 1) + vec4(cos(theta), 0, sin(theta), 0)*radius;
}

vec4 circleZ(float theta, vec3 center, float radius) {
	return vec4(center, 1) + vec4(cos(theta), sin(theta), 0, 0)*radius;
}

void generateBlendedCurves(vector<Drawable>& drawables, int numCurves, int numDivisions, float radius) {
	float thetaStep = 2.f*M_PI / float(numCurves - 1);
	float theta = 0;

	vec3 centerY(0, -1, 0);
	vec3 centerX(1, 0, 0);
	vec3 centerZ(0, 0, 1);

	vector<bezier<vec4>> curvesX;
	vector<bezier<vec4>> curvesZ;

	for (int i = 0; i < numCurves; i++) {
		float xTheta = -theta - M_PI*0.5f;
		float zTheta = -theta;
		vec4 yPoint = circleY(theta, centerY, radius);
		vec4 xPoint = circleX(xTheta, centerX, radius);
		vec4 zPoint = circleZ(zTheta, centerZ, radius);
		float wX = 1.f / max(cos(theta), 0.0001f);
		float wY = 1.f / max(cos(theta - M_PI*0.5f), 0.0001f);
		vec4 middleX = vec4(yPoint.x, xPoint.y, xPoint.z, 1.0)*wX;
		vec4 middleZ = vec4(zPoint.x, zPoint.y, xPoint.z, 1.0)*wY;

		curvesX.push_back(bezier<vec4>({ yPoint, middleX, xPoint }));
		curvesZ.push_back(bezier<vec4>({ yPoint, middleZ, zPoint }));

		theta += thetaStep;
	}

	for (int i = 0; i < curvesX.size(); i++) {
		vector<vec4> points4D = curvesX[i].getQuadPoints(50);
		vector<vec3> pointsX;
		for (int j = 0; j < points4D.size(); j++) {
			const vec4& p = points4D[j];
			pointsX.push_back(vec3(p.x / p.w, p.y / p.w, p.z / p.w));
		}

		points4D = curvesZ[i].getQuadPoints(50);
		vector<vec3>pointsZ;
		for (int j = 0; j < points4D.size(); j++) {
			const vec4& p = points4D[j];
			pointsZ.push_back(vec3(p.x / p.w, p.y / p.w, p.z / p.w));
		}

		vector<vec3> points;
		float s = 0;				//s parameterizes line from pointY to center
		float sStep = 0.75 / 50.f;
		vec3 a = curvesX[i].control[0] / curvesX[i].control[0].w;
		vec3 bX = curvesX[i].control[1] / curvesX[i].control[1].w;
		vec3 bZ = curvesZ[i].control[1] / curvesZ[i].control[1].w;

		float xRatio = 1.f/length(a - bX);
		float zRatio = 1.f / length(a - bZ);

		while(s < 0.75) {
			vec3 p = a + s*vec3(0, 1, 0);
			float uX = 1.f - sqrt(1.f - s*xRatio);
			float uZ = 1.f - sqrt(1.f - s*zRatio);
//			vec3 pX = cu

			s += sStep;
		}


		vector<vec3> controlX;
		vector<vec3> controlZ;
		for (int j = 0; j < curvesX[i].control.size(); j++) {
			const vec4& pX = curvesX[i].control[j];
			const vec4& pZ = curvesZ[i].control[j];
			controlX.push_back(vec3(pX.x / pX.w, pX.y / pX.w, pX.z / pX.w) + vec3(0.003, 0.003, 0.003));
			controlZ.push_back(vec3(pZ.x / pZ.w, pZ.y / pZ.w, pZ.z / pZ.w) + vec3(0.003, 0.003, 0.003));
		}
		drawables.push_back(Drawable(
			new ColorMat(vec3(0.f, 0.f, 1.f)),
			new SimpleGeometry(controlX.data(), controlX.size(), GL_LINE_STRIP)));
		drawables.push_back(Drawable(
			new ColorMat(vec3(0.f, 0.f, 1.f)),
			new SimpleGeometry(controlZ.data(), controlZ.size(), GL_LINE_STRIP)));
		
		/*drawables.push_back(Drawable(
			new ColorMat(vec3(1.f, 0.f, 0.f)),
			new SimpleGeometry(pointsX.data(), pointsX.size(), GL_LINE_STRIP))); 
		drawables.push_back(Drawable(
			new ColorMat(vec3(0.f, 1.f, 0.f)),
			new SimpleGeometry(pointsZ.data(), pointsZ.size(), GL_LINE_STRIP)));*/
	}
}
	
void WindowManager::mainLoop() {

	glfwSetKeyCallback(window, keyCallback);
	glfwSetWindowSizeCallback(window, windowResizeCallback);
	//glfwSetCursorPosCallback(window, cursorPositionCallback);

	vector<Drawable> drawables;

	generateBlendedCurves(drawables, 40, 50, 0.1f);

	SimpleShader shader;

	vec3 lightPos(10.f, 10.f, 10.f);

	while (!glfwWindowShouldClose(window)) {

		if (windowResized) {
			window_width = windowWidth;
			window_height = windowHeight;
		}

		glClearColor(0.f, 0.f, 0.f, 1.f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		for (int i = 0; i < drawables.size(); i++) {
			shader.draw(cam, drawables[i]);
		}

		glfwSwapBuffers(window);
		glfwWaitEvents();
	}

	glfwTerminate();
}

void initGLExtensions() {
#ifndef USING_GLEW
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
		std::cout << "GLAD initialization failed" << std::endl;
	}
#else
	glewExperimental = true;
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		std::cout << "GLEW initialization failed" << std::endl;
	}
#endif
}


GLFWwindow *createWindow(int width, int height, std::string name)
{
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	
	GLFWwindow *window = glfwCreateWindow(
		width, height, name.c_str(), nullptr, nullptr);
	
	if (window == NULL) {
		glfwTerminate();
		return nullptr;
	}
	else {
		glfwMakeContextCurrent(window);
		return window;
	}
}

