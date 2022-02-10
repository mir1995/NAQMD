LINEAR_DIR = ../Linear
METRICS_DIR = ../../src/Auxiliary
GRID_DIR = ../../src/Grid
QUAD_DIR = ../../src/Quadrature
HAG_DIR = ../../src/Hagedorn
POT_DIR = ../../src/Potential
ODE_DIR = ../../src/Odeint
OBJ_DIR = ./Object

EXE = hagedorn_linear# hagedorn_linear_projection_comparison   
# should not be a problem # it would be good to keep a separate file with the dependencies
SRCS = $(wildcard hagedorn_linear.c $(HAG_DIR)/*.c $(GRID_DIR)/*.c $(POT_DIR)/potential.c $(POT_DIR)/pot_linear.c $(POT_DIR)/potential_io.c $(QUAD_DIR)/*.c $(METRICS_DIR)/*.c ) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS))) 

vpath %.c $(LINEAR_DIR) $(HAG_DIR) $(POT_DIR) $(QUAD_DIR) $(GRID_DIR) $(METRICS_DIR) 

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
