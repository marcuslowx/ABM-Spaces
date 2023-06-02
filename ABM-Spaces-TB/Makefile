BASE=abm_spaces_tb7
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
	g++ -g -std=c++11 -Wall -pedantic abm_spaces_tb7.cc -o $(EXE) -lpthread

analyze: release
	./$(EXE) --num_runs=5 --cores=$(CORES) >a.csv
	./$(EXE) --num_runs=5 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --num_runs=5 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --num_runs=5 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --num_runs=5 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

analyze100: release
	./$(EXE) --num_runs=100 --cores=$(CORES) >a.csv
	./$(EXE) --num_runs=100 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --num_runs=100 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --num_runs=100 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --num_runs=100 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

analyze1000: release
	./$(EXE) --num_runs=1000 --cores=$(CORES) >a.csv
	./$(EXE) --num_runs=1000 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --num_runs=1000 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --num_runs=1000 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --num_runs=1000 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

analyze1000alt_activation: release
	./$(EXE) --alt_activation=1 --num_runs=1000 --cores=$(CORES) >a.csv
	./$(EXE) --alt_activation=1 --num_runs=1000 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --alt_activation=1 --num_runs=1000 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --alt_activation=1 --num_runs=1000 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --alt_activation=1 --num_runs=1000 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

analyze1000alt_non_infectious: release
	./$(EXE) --alt_non_infectious=1 --num_runs=1000 --cores=$(CORES) >a.csv
	./$(EXE) --alt_non_infectious=1 --num_runs=1000 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --alt_non_infectious=1 --num_runs=1000 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --alt_non_infectious=1 --num_runs=1000 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --alt_non_infectious=1 --num_runs=1000 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

csv: release
	./$(EXE) --num_runs=10 --mass_x_ray=0 --cores=$(CORES) >tb_test.csv
	R -f lookcsv.R
	$(DefBrowser) TB_C.html

analyze20years: release
	./$(EXE) --num_time=1040 --num_runs=10 --cores=$(CORES) >a.csv
	./$(EXE) --num_time=1040 --num_runs=10 --tracing=1 --cores=$(CORES) > b.csv
	./$(EXE) --num_time=1040 --num_runs=10 --mass_x_ray=1 --tracing=1 --cores=$(CORES) > c.csv
	./$(EXE) --num_time=1040 --num_runs=10 --tut=1 --tracing=1 --cores=$(CORES) > d.csv
	./$(EXE) --num_time=1040 --num_runs=10 --tut=1 --tracing=1 --mass_x_ray=1 --cores=$(CORES) > e.csv
	R -f look.R
	$(DefBrowser) five.html

fit: release
	./$(EXE) --num_runs=1000 --cores=$(CORES) --testing_rate=0.647 --house_risk=0.0127 --class_risk=0.0076 --block_risk=0.00115 --taxi_risk=0.022 --work_risk=0.1851 >fit.csv
	R -f fit.R
	$(DefBrowser) fit.html



