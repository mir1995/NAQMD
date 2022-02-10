QUAD_DIR = ../src/Quadrature
HAG_DIR = ../src/Hagedorn
OBJ_DIR = ./Object

EXE = test_hermite_hagedorn
SRCS = $(wildcard test_hermite_hagedorn.c $(HAG_DIR)/*.c $(QUAD_DIR)/*.c) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS))) 

vpath %.c $(HAG_DIR) $(QUAD_DIR) 

CC = gcc
CFLAGS = -pedantic -Wall -g 
LDLIBS = -lm

#HEADERS = $(SH_DIR)/%.h $(POT_DIR)/%.h $(ODE_DIR)/%.h ?

.PHONY: all run clean


all: $(EXE)

$(EXE):  $(OBJS) # LINKING STAGE
	$(CC) -g $^ -o $@ $(LDLIBS) 

$(OBJ_DIR)/%.o: %.c  
	$(CC) $(CFLAGS) -c $< -o $@ 

run: $(EXE)
	./$<

clean:
	rm $(OBJ_DIR)/*.o
