CC = mpicc
EXEC = mpiexec
OPT = -g -Wall
PARRALLEL_SRC = mpi_mat_vect_col.c
PARRALLEL_TARGET = parallel

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o  $(PARRALLEL_TARGET) -lm

run: $(PARRALLEL_TARGET)
	$(EXEC) -n 2 $(PARRALLEL_TARGET) 

clean:
	rm $(PARRALLEL_TARGET)
