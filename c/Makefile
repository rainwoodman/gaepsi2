OBJS=quadcurve.o spline.o svremap.o transform.o sphrasterize.o

libgaepsi.a: $(OBJS)
	ar r $@ $^
	ranlib $@
.c.o:
	$(CC) -c -o $@ $<

clean:
	rm *.o *.a
