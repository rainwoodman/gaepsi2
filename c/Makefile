OBJS=quadcurve.o kernel.o svremap.o transform.o sphrasterize.o

libgaepsi.a: $(OBJS)
	ar r $@ $^
	ranlib $@
.c.o:
	$(CC) -c -o $@ $<

clean:
	rm *.o *.a
