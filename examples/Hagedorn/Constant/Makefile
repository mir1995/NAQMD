CONSTANT_DIR = ../Constant
HAG_DIR = ../../../src/Hagedorn
POT_DIR = ../../../src/Potential
ODE_DIR = ../../../src/Odeint
OBJ_DIR = ./Object

EXE = hag_constant   
# should not be a problem # it would be good to keep a separate file with the dependencies
SRCS = $(wildcard hag_constant.c $(HAG_DIR)/hag_dynamics.c $(POT_DIR)/potential.c $(POT_DIR)/pot_constant.c $(ODE_DIR)/odeint.c) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS))) 

vpath %.c $(CONSTANT_DIR) $(HAG_DIR) $(POT_DIR) $(ODE_DIR) 

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
