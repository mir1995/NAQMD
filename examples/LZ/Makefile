LZ_DIR = ../LZ
SH_DIR = ../../src/SurfaceHopping
POT_DIR = ../../src/Potential
ODE_DIR = ../../src/Odeint
AUX_DIR = ../../src/Auxiliary
OBJ_DIR = ./Object

EXE = lz_mass
# should not be a problem # it would be good to keep a separate file with the dependencies
SRCS = $(wildcard $(EXE).c $(SH_DIR)/*.c $(POT_DIR)/potential.c $(POT_DIR)/pot_lz.c $(ODE_DIR)/*.c $(AUX_DIR)/*.c) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS)) lz_mass.c) 

vpath %.c $(LZ_DIR) $(SH_DIR) $(POT_DIR) $(ODE_DIR) $(AUX_DIR)

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
