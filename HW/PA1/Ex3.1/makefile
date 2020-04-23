CC = mpicc
EXEC = srun
OPT = -g -Wall
PARRALLEL_SRC = prog3.1_histo_dist.c
PARRALLEL_TARGET = histogram

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o $(PARRALLEL_TARGET)
run: $(PARRALLEL_TARGET)
	$(EXEC) -n 4 $(PARRALLEL_TARGET) < input
clean:
	rm $(PARRALLEL_TARGET)
	