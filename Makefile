all: fluid display

fluid: fluid.cpp fluid.h
	g++-14 -std=c++17 -fopenmp -I/usr/local/Cellar/libomp/19.1.7/include -L/usr/local/Cellar/libomp/19.1.7/lib -I/usr/local/include -L/usr/local/lib -I/usr/local/Cellar/sdl2/2.30.12/include -L/usr/local/Cellar/sdl2/2.30.12/lib -lSDL2 -lomp -o fluid fluid.cpp
	# g++ -std=c++17 -o fluid fluid.cpp -lSDL2

display: simDisplayer.cpp
	g++-14 -std=c++17 -I/usr/local/include -L/usr/local/lib -I/usr/local/Cellar/sdl2/2.30.12/include -L/usr/local/Cellar/sdl2/2.30.12/lib -o display simDisplayer.cpp -lSDL2