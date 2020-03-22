CC = mpicc
STDCC = gcc
EXEC = mpiexec
OPT = -g -Wall
PARRALLEL_SRC = prefix_sum.c
PARRALLEL_TARGET = parallel
SERIAL_SRC = serial_prefix_sum.c
SERIAL_TARGET = serial

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o $(PARRALLEL_TARGET)

$(SERIAL_TARGET): $(SERIAL_SRC)
	$(STDCC) $(OPT) $(SERIAL_SRC) -o $(SERIAL_TARGET)

run_parallel: $(PARRALLEL_TARGET)
	$(EXEC) -n 4 $(PARRALLEL_TARGET) -g
run_serial: $(SERIAL_TARGET)
	./$(SERIAL_TARGET)

check: $(PARRALLEL_TARGET) $(SERIAL_TARGET)
	python check.py ./$(SERIAL_TARGET) "$(EXEC) -n 4 $(PARRALLEL_TARGET) -i" "python datamaker.py"

clean:
	rm $(PARRALLEL_TARGET) $(SERIAL_TARGET)
