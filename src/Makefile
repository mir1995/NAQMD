NAI_DIR = ./NaI
SH_DIR = ./SurfaceHopping
POT_DIR = ./Potential
ODE_DIR = ./Odeint
OBJ_DIR = ./Object

# you want to have a file of target-dependencies where each gets called 
# depending on the EXE
SRCS = $(wildcard $(NAI_DIR)/*.c $(SH_DIR)/*.c $(POT_DIR)/pot_nai.c $(ODE_DIR)/*.c) 
OBJS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(notdir $(SRCS))) 

vpath %.c $(NAI_DIR) $(SH_DIR) $(POT_DIR) $(ODE_DIR) 

CC = gcc
CFLAGS = -Wall -g 
LDLIBS = -lm

EXE = main
#HEADERS = $(SH_DIR)/%.h $(POT_DIR)/%.h $(ODE_DIR)/%.h ?

.PHONY: all run clean


all: $(EXE)

$(EXE):  $(OBJS) # LINKING STAGE
	$(CC) -g $^ -o $@ $(LDLIBS) 

$(OBJ_DIR)/%.o: %.c  # generates the .o files - at this stage it does not really matter if there is no implr
	$(CC) $(CFLAGS) -c $< -o $@ 

run: $(EXE)
	./$<

clean:
	rm $(OBJ_DIR)/*.o
