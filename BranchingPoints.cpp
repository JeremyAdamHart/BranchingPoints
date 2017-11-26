// RenderingEngine.cpp : Defines the entry point for the console application.
//
#include "GLWindow.h"
#include <cstring>

int main(int argc, char **argv)
{
	WindowManager wm (800, 800, "Branching Points");

	wm.mainLoop();
}

