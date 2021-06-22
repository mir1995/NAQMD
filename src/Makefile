SH_DIR = ./SurfaceHopping
POT_DIR = ./Potential
ODE_DIR = ./Odeint
OBJ_DIR = ./Object

SRCS = $(wildcard $(SH_DIR)/*.c $(POT_DIR)/*.c $(ODE_DIR)/*.c) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS))) 

vpath %.c $(SH_DIR) $(POT_DIR) $(ODE_DIR)

CC = gcc
CFLAGS = -Wall -g 
LDLIBS = -lm

EXE = test_evolution
#HEADERS = $(SH_DIR)/%.h $(POT_DIR)/%.h $(ODE_DIR)/%.h ?

.PHONY: all run clean


all: $(EXE)

$(EXE):  $(OBJS)
	$(CC) -g $^ -o $@ $(LDLIBS) 

$(OBJ_DIR)/%.o: %.c  
	$(CC) $(CFLAGS) -c $< -o $@ 

run: $(EXE)
	./$<

clean:
	rm $(OBJ_DIR)/*.o
