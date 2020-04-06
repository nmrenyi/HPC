CC = gcc
STDCC = gcc
EXEC = ./
OPT = -g -Wall 
PARRALLEL_SRC = main.c
PARRALLEL_TARGET = main

$(PARRALLEL_TARGET): $(PARRALLEL_SRC)
	$(CC) $(OPT) $(PARRALLEL_SRC) -o $(PARRALLEL_TARGET) -lpthread

run: $(PARRALLEL_TARGET)
	$(EXEC)$(PARRALLEL_TARGET)

check: $(PARRALLEL_TARGET) $(SERIAL_TARGET)
	python check.py "$(EXEC)$(PARRALLEL_TARGET)" "python datamaker.py"

clean:
	rm $(PARRALLEL_TARGET)
