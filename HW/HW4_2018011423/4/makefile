CC = gcc
STDCC = gcc
EXEC = ./
OPT = -g -Wall 
PARRALLEL_SRC = parallel_input.c
PARRALLEL_TARGET = parallel
SERIAL_SRC = serial.c
SERIAL_TARGET = serial
TASKS_COUNT = 60
THREADS = 10

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o $(PARRALLEL_TARGET) -lpthread

$(SERIAL_TARGET): $(SERIAL_SRC)
	$(STDCC) $(OPT) $(SERIAL_SRC) -o $(SERIAL_TARGET) -lpthread

run: $(PARRALLEL_TARGET)
	$(EXEC) $(THREADS) $(TASKS_COUNT)
run_serial: $(SERIAL_TARGET)
	./$(SERIAL_TARGET) $(TASKS_COUNT)

check: $(PARRALLEL_TARGET) $(SERIAL_TARGET)
	python check.py "./$(SERIAL_TARGET) $(TASKS_COUNT)" "$(EXEC)$(PARRALLEL_TARGET) $(THREADS) $(TASKS_COUNT)" "python datamaker.py"

clean:
	rm $(PARRALLEL_TARGET) $(SERIAL_TARGET)
