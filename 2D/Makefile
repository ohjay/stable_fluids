CFLAGS = -Wno-deprecated-declarations
UNAME_S := $(shell uname -s)

ifeq ($(shell uname -s), Darwin)
LDFLAGS= -framework GLUT -framework OpenGL -framework Cocoa
else
LDFLAGS= -lGL -lGLU -lglut
endif

sim: Fluid.h solver.h
	g++ $(CFLAGS) -o sim main.cpp Fluid.cpp solver.cpp $(LDFLAGS)

clean:
	$(RM) sim
