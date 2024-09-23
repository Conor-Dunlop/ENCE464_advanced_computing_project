all: poisson

# -g outputs debugging information
# -Wall enables all warnings
# -pthread configures threading
CFLAGS = -pg -g -Wall -pthread
LDFLAGS += $(CFLAGS)
CXXFLAGS = $(CFLAGS)
VPATH += implementations

poisson: main.o poisson_r2.o

.PHONY: disassembly
disassembly: poisson.s

poisson.s: poisson
	objdump -S --disassemble $< > $@

.PHONY: test
test: poisson
	./test.sh

.PHONY: clean
clean:
	rm -f poisson *.o *.s



