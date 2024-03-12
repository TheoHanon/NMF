CC = gcc
CFLAGS = -Wall -Wextra

# Base directory for data files
DATA_DIR = post_processing/data
# Configuration for running the program
N = 128
SCHEME_TYPE = ED # Options: E2, E4, I4, ED
GRID_TYPE = UNIF # Options: UNIF, NONUNIF
INITIAL_TYPE = GAUSSIAN # Options: GAUSSIAN, WAVEPACKET

# Filenames for output (these need to be aligned with your program's expectations)
FILE_NAME1 = data_$(N).txt
FILE_NAME2 = diagnose_$(N).txt
FILE_NAME3 = grid_$(N).txt
FILE_NAME4 = param_$(N).txt

# Compile the source code into an executable	
hw1: src/hw1.c
	$(CC) $(CFLAGS) -o hw1 src/*.c

# Run the program with a specific configuration
run: hw1
	./hw1 $(N) $(SCHEME_TYPE) $(GRID_TYPE) $(INITIAL_TYPE) $(DATA_DIR)/$(FILE_NAME1) $(DATA_DIR)/$(FILE_NAME2) $(DATA_DIR)/$(FILE_NAME3) $(DATA_DIR)/$(FILE_NAME4)
	python ./post_processing/animation.py -N $(N) 


plot: 
	python ./post_processing/PARTI.py
	python ./post_processing/PARTII.py
	python ./post_processing/PARTIII.py
# Print usage instructions
help:
	@echo "Usage instructions for running the program:"
	@echo "  make run N=<grid size> SCHEME_TYPE=<scheme type> GRID_TYPE=<grid type> INITIAL_TYPE=<initial condition>"
	@echo "  Available SCHEME_TYPE options: E2, E4, I4, ED"
	@echo "  Available GRID_TYPE options: UNIF, NONUNIF"
	@echo "  Available INITIAL_TYPE options: GAUSSIAN, WAVEPACKET"
	@echo "Example:"
	@echo "  make run N=128 SCHEME_TYPE=ED GRID_TYPE=NONUNIF INITIAL_TYPE=GAUSSIAN"

# Clean up generated files
clean:
	rm -f hw1
	rm -f $(DATA_DIR)/*.txt




