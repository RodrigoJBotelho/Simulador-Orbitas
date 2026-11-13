CXX ?= c++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra

all: bin/orbitas

bin/orbitas: orbitas/orbitas.cpp
	@mkdir -p bin
	$(CXX) $(CXXFLAGS) -o $@ $<

run: all
	@mkdir -p saida
	cd saida && ../bin/orbitas

teste: all
	@mkdir -p saida
	cd saida && ../bin/orbitas --verificar

exemplos: run
	cp saida/*.csv docs/exemplos/

clean:
	rm -rf bin saida

.PHONY: all run teste exemplos clean
