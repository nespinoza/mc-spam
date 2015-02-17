INC = /usr/local/include
LIB = /usr/local/lib
CC  = gcc
OBJS = rv.o true_anomaly.o orbit.o transit_times.o convlib.o
AR = ar
LIB_EXT = .a
RANLIB = ranlib

test: $(OBJS) planet.h test.c
	$(CC) -o test test.c $(OBJS) -I$(INC) -L$(LIB) -lgsl -lgslcblas -lm

convlib.o: planet.h convlib.c
	$(CC) -c convlib.c -I$(INC)

rv.o: true_anomaly.o planet.h rv.c
	$(CC) -c rv.c -I$(INC)

true_anomaly.o: planet.h true_anomaly.c
	$(CC) -c true_anomaly.c -I$(INC)

orbit.o: planet.h orbit.c
	$(CC) -c orbit.c -I$(INC)

transit_times.o: planet.h transit_times.c
	$(CC) -c transit_times.c -I$(INC)

lib: $(OBJS)
	$(AR) rv libplanet$(LIB_EXT) $(OBJS)
	$(RANLIB) libplanet$(LIB_EXT)


