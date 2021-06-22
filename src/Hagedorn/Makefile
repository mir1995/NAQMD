		
CC=gcc
CFLAGS += -std=c99 -pedantic -Wall -Wextra

.PHONY : all 
all : $() particle 

## 
test: test.o 
		$(CC) $(CFLAGS) $< -o $@ -lm

##
test.o : test.c hag_wavepackets.c hagedorn.h 
		$(CC) $(CFLAGS) -c $< -lm 
		
.PHONY : clean 
clean :
		rm -f *.o # replace this with variable 
