TARGET = a.out

OBJECTS = read_input_nnp2.o \
          read_io_nnp2.o \
          read_io2.o \
          main_mean_disp.o 

MOD_FILES = 

FC = ifort -O1 -shared-intel -mcmodel=large

FFLAGS = 
LDFLAGS = 

.SUFFIXES: .o .f90
.f90.o:
	${FC} -c $<

${TARGET}: ${OBJECTS}
	${FC} -o $@ ${OBJECTS} ${LDFLAGS} ${FFLAGS}

	cp ${TARGET} ../

.PHONY: clean
clean:
	${RM} ${TARGET} ${OBJECTS} ${MOD_FILES}
