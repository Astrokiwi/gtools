FORT = gfortran
#OPTIONS = -O3 -Wall -mcmodel=medium
OPTIONS = -g -Wall -mcmodel=medium -fbacktrace

FFLAGS = ${OPTIONS}

%.o: %.mod

.SUFFIXES : .o .f90 .F

DISC_OBJS=gadget_ics_IO.o disc_ICs.o
DISC_EXEC=disc_ICs

disc : ${DISC_OBJS}
	${FORT} ${OPTIONS} -o ${DISC_EXEC} ${DISC_OBJS} ${LIBS}
	mv ${DISC_EXEC} ..


DROP_OBJS=gadget_ics_IO.o blobdrop_ICs.o
DROP_EXEC=blobdrop_ICs

drop : ${DROP_OBJS}
	${FORT} ${OPTIONS} -o ${DROP_EXEC} ${DROP_OBJS} ${LIBS}
	mv ${DROP_EXEC} ..


BLOB_OBJS=gadget_ics_IO.o blob_ICs.o
BLOB_EXEC=blob_ICs

blob : ${BLOB_OBJS}
	${FORT} ${OPTIONS} -o ${BLOB_EXEC} ${BLOB_OBJS} ${LIBS}
	mv ${BLOB_EXEC} ..

PERIODIC_OBJS=gadget_ics_IO.o periodic_ICs.o
PERIODIC_EXEC=periodic_ICs

periodic : ${PERIODIC_OBJS}
	${FORT} ${OPTIONS} -o ${PERIODIC_EXEC} ${PERIODIC_OBJS} ${LIBS}
	mv ${PERIODIC_EXEC} ..


TEST_OBJS=gadget_ics_IO.o opacity_test.o
TEST_EXEC=opacity_test_ICs

test : ${TEST_OBJS}
	${FORT} ${OPTIONS} -o ${TEST_EXEC} ${TEST_OBJS} ${LIBS}
	mv ${TEST_EXEC} ..

.f90.o :
	${FORT} ${FFLAGS} -c $<

.F.o :
	${FORT} ${FFLAGS} -c $<


clean :
	rm *.o *.mod
