all: poisson

# -g outputs debugging information
# -Wall enables all warnings
# -pthread configures threading
CCX=g++
CXXFLAGS = -std=c++20 -g -Wall -pthread -O3

poisson: main.cpp poisson_mt.cpp calc_scalar.cpp
	$(CCX) $^ -o $@ $(CXXFLAGS)

.PHONY: disassembly
disassembly: poisson.s

poisson.s: poisson
	objdump -S --disassemble $<> $@

.PHONY: test
test: poisson
	./test.sh

.PHONY: clean
clean:
	rm -f poisson *.o *.s
