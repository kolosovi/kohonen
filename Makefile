CC = mpicc
CFLAGS = -O2 -Wall -Werror -Wformat-security -Wignored-qualifiers -Winit-self -Wswitch-default -Wfloat-equal -Wpointer-arith -Wtype-limits -Wempty-body -Wstrict-prototypes -Wold-style-definition -Wmissing-field-initializers -Wnested-externs -Wno-pointer-sign -std=gnu99

SOURCEDIR = src
BUILDDIR = build


build_dir:
	mkdir -p $(BUILDDIR)


kohonen_learn : build_dir $(BUILDDIR)/kohonen_learn.o $(BUILDDIR)/util.o
	$(CC) $(CFLAGS) $(BUILDDIR)/kohonen_learn.o $(BUILDDIR)/util.o -o kohonen_learn

$(BUILDDIR)/util.o : build_dir $(SOURCEDIR)/util.c
	$(CC) $(CFLAGS) -c $(SOURCEDIR)/util.c -o $(BUILDDIR)/util.o

$(BUILDDIR)/kohonen_learn.o : build_dir $(SOURCEDIR)/kohonen_learn.c
	$(CC) $(CFLAGS) -c $(SOURCEDIR)/kohonen_learn.c -o $(BUILDDIR)/kohonen_learn.o

$(BUILDDIR)/tests.o : build_dir $(SOURCEDIR)/tests.c
	$(CC) $(CFLAGS) -c $(SOURCEDIR)/tests.c -o $(BUILDDIR)/tests.o

$(BUILDDIR)/test_binary : build_dir $(BUILDDIR)/tests.o $(BUILDDIR)/util.o
	$(CC) $(CFLAGS) $(BUILDDIR)/tests.o $(BUILDDIR)/util.o -o $(BUILDDIR)/test_binary

test : build_dir $(BUILDDIR)/test_binary
	./$(BUILDDIR)/test_binary
