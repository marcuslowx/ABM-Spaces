BASE=abm_tests
EXE=$(BASE)
CORES=0
DefBrowser=explorer.exe
ifeq ($(OS),Windows_NT)
	EXE=$(BASE)+.exe
	DefBrowser=explorer.exe
endif

release:
	g++ -O3 -std=c++11 -DNDEBUG $(BASE).cc -o $(EXE) -lpthread

debug:
	g++ -g -std=c++11 -Wall -pedantic abm_tests.cc -o $(EXE) -lpthread

lookcsv: release
	./$(EXE) --num_runs=1 --reg_tests=1 --cores=$(CORES) >tests.csv
	R -f lookcsv.R

tests: release
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=4 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a4.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=8 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=12 --tracing=1 --reg_tests_sensitivity=0.7 --cores=$(CORES) > a12.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=16 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a16.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=20 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a20.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=24 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a24.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=28 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a28.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=32 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a32.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=36 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a36.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=40 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a40.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=44 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a44.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=48 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a48.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=52 --reg_tests_sensitivity=0.7 --tracing=1 --cores=$(CORES) > a52.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=4 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a4s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=8 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a8s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=12 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a12s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=16 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a16s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=20 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a20s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=24 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a24s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=28 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a28s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=32 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a32s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=36 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a36s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=40 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a40s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=44 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a44s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=48 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a48s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=52 --reg_tests_sensitivity=0.8 --tracing=1 --cores=$(CORES) > a52s8.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=4 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a4s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=8 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a8s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=12 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a12s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=16 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a16s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=20 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a20s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=24 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a24s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=28 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a28s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=32 --tracing=1 --cores=$(CORES) > a32s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=36 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a36s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=40 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a40s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=44 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a44s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=48 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a48s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=52 --reg_tests_sensitivity=0.9 --tracing=1 --cores=$(CORES) > a52s9.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=4 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a4s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=8 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a8s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=12 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a12s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=16 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a16s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=20 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a20s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=24 --tracing=1 --cores=$(CORES) > a24s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=28 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a28s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=32 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a32s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=36 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a36s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=40 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a40s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=44 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a44s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=48 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a48s1.csv
	./$(EXE) --num_runs=1000 --reg_tests=1 --reg_tests_freq=52 --reg_tests_sensitivity=1 --tracing=1 --cores=$(CORES) > a52s1.csv
	R -f tests.R

