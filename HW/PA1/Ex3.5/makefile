CC = mpicc
EXEC = srun
OPT = -g -Wall
PARRALLEL_SRC = mpi_mat_vect3-5.c
PARRALLEL_TARGET = parallel

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o  $(PARRALLEL_TARGET) -lm 

run: $(PARRALLEL_TARGET)
	$(EXEC) -n 4 $(PARRALLEL_TARGET) 

clean:
	rm $(PARRALLEL_TARGET)
