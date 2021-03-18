CC = gcc
H = types.h io.h
C = sam2ain.c io.c
EXEC = sam2main

build: ${C} ${H}
	${CC} ${C} -lm -o ${EXEC} 

debug: ${C} ${H}
	${CC} ${C} -lm -g -o ${EXEC} 

clean:
	rm -f ${EXEC} *.csv