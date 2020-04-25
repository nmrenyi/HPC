CC = mpicc
EXEC = mpiexec
OPT = -g -Wall
PARRALLEL_SRC = mpi_mat_vect3-6.c
PARRALLEL_TARGET = parallel
THREAD_COUNT = 10
$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o  $(PARRALLEL_TARGET) -lm -fopenmp -lpthread 

run: $(PARRALLEL_TARGET)
	$(EXEC) -n 4 $(PARRALLEL_TARGET) $(THREAD_COUNT)

clean:
	rm $(PARRALLEL_TARGET)
