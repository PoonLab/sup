CC = gcc
H = types.h io.h
C = sam2aln.c io.c
EXEC = sam2aln
CFLAGS = -Wall -Werror -g -std=c99

build: ${C} ${H}
	${CC} $(CFLAGS) `pkg-config --cflags glib-2.0` ${C} -lm -o ${EXEC} `pkg-config --libs glib-2.0` -pthread

debug: ${C} ${H}
	${CC} `pkg-config --cflags glib-2.0` ${C} -lm -g -o ${EXEC} `pkg-config --libs glib-2.0` -pthread

clean:
	rm -f ${EXEC} *.csv
