CC = mpicc
STDCC = gcc
EXEC = srun
OPT = -g -Wall
PARRALLEL_SRC = mpi_vector_add_v.c
PARRALLEL_TARGET = parallel
SERIAL_SRC = vector_add.c
SERIAL_TARGET = serial

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o $(PARRALLEL_TARGET)

$(SERIAL_TARGET): $(SERIAL_SRC)
	$(STDCC) $(OPT) $(SERIAL_SRC) -o $(SERIAL_TARGET)

run:$(PARRALLEL_TARGET)
	$(EXEC) -n 4 $(PARRALLEL_TARGET)

run_serial: $(SERIAL_TARGET)
	./$(SERIAL_TARGET)

check: $(PARRALLEL_TARGET) $(SERIAL_TARGET)
	python check.py ./$(SERIAL_TARGET) "$(EXEC) -n 4 $(PARRALLEL_TARGET)" "python datamaker.py"

clean:
	rm $(PARRALLEL_TARGET) $(SERIAL_TARGET)
