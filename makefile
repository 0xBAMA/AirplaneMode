FLAGS = -O3 -std=c++17 -lpthread
all: render

render: src/main.cc src/lib/*.h
		g++ -o render src/main.cc ${FLAGS}

run:
		./render out.png
