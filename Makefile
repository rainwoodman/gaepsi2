CC=cc
LDSHARED=$(CC) -shared
all:
	LDSHARED="$(LDSHARED)" CC="$(CC)" python setup.py build_ext --inplace
