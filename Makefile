CC = mpicc
CFLAGS = -O2 -Wall -Werror -Wformat-security -Wignored-qualifiers -Winit-self -Wswitch-default -Wfloat-equal -Wpointer-arith -Wtype-limits -Wempty-body -Wstrict-prototypes -Wold-style-definition -Wmissing-field-initializers -Wnested-externs -Wno-pointer-sign -std=gnu99


kohonen_learn : kohonen_learn.o util.o
	$(CC) $(CFLAGS) kohonen_learn.o util.o -o kohonen_learn

util.o : util.c
	$(CC) $(CFLAGS) -c util.c

kohonen_learn.o : kohonen_learn.c
	$(CC) $(CFLAGS) -c kohonen_learn.c

tests.o : tests.c
	$(CC) $(CFLAGS) -c tests.c

test_binary : tests.o util.o
	$(CC) $(CFLAGS) tests.o util.o -o test_binary 

test : test_binary
	./test_binary
