CC=g++


CFLAGS=-c -std=c++11 -O3 -Wall -g
LINKFLAGS=-O3

#debug = true
ifdef debug
	CFLAGS +=-g
	LINKFLAGS += -flto
endif

INCDIR= -I./include

LIBDIR=-L/usr/X11R6 -L/usr/local/lib -L../RenderingEngineLibrary

LIBS=

OS_NAME:=$(shell uname -s)

ifeq ($(OS_NAME),Darwin)
	LIBS += `pkg-config --static --libs glfw3 gl`
	LIBS += \
	 -lglfw3 \
	 -lGLEW \
	 -framework Cocoa \
	 -framework OpenGL \
	 -framework IOKit \
	-framework CoreVideo
endif
ifeq ($(OS_NAME),Linux)
	LIBS += `pkg-config --static --libs glfw3 gl`
endif

ENGINE_DIR=.

ENGINE_SOURCES=$(wildcard $(ENGINE_DIR)/*cpp)

ENGINE_HEADERS=-I$(ENGINE_DIR) 

HEADERS=-Iheaders -Iinclude

ENGINE_OBJS=$(addprefix $(ENGINE_DIR)/obj/,$(notdir $(ENGINE_SOURCES:.cpp=.o)))

.PHONY: default
default:  fetch buildDirectories $(ENGINE_DIR)/BranchingPoints.out
	echo "Build completed"

$(ENGINE_DIR)/BranchingPoints.out: $(ENGINE_OBJS)
	$(CC) $(LINK_FLAGS) $(ENGINE_OBJS) -o $@ -L. -lrenderingextensions -lrenderingcore $(LIBS) $(LIBDIR)

all: $(ENGINE_DIR)/BranchingPoints.out 

$(ENGINE_DIR)/obj/%.o: $(ENGINE_DIR)/%.cpp
	$(CC) -c $(CFLAGS) $(ENGINE_HEADERS) $(HEADERS) $(INCDIR) $(LIBDIR) $< -o $@

.PHONY: fetch
fetch: 
	if [ ! -d "./shaders" ]; then \
		svn export https://github.com/JeremyAdamHart/RenderingEngineLibrary/trunk/RenderingEngine/shaders; \
	fi
	if [ ! -d "./headers" ]; then \
		svn export https://github.com/JeremyAdamHart/RenderingEngineLibrary/trunk/RenderingLibrary/headers && \
		svn export https://github.com/JeremyAdamHart/RenderingEngineLibrary/trunk/RenderingExtensions/headers/ --force ; \
	fi


.PHONY: buildDirectories
buildDirectories:
	mkdir -p $(ENGINE_DIR)/obj

.PHONY: clean
clean:
	rm $(ENGINE_DIR)/*.out $(ENGINE_DIR)/obj/*.o
