CXX = g++
CXXFLAGS = -std=c++17

analyzer: main.cpp fasta.cpp fasta.h
	$(CXX) $(CXXFLAGS) -o analyzer main.cpp fasta.cpp

clean:
	rm -f analyzer