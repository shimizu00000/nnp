TARGET = a.out

OBJECTS = allocarray.o pre_random.o \
          read_input_nnp2.o read_io2.o read_layer_node2.o read_layer_node_frex2.o \
          energy_force_nnp2ab.o f_ef_deriv2ab.o \
          train2ab.o valid2ab.o \
          read_data_1.o read_data_2.o \
          timer.o blas.o linpack.o lbfgsb.o \
          main_training.o

MOD_FILES = allocarray.mod

FC = mpif90 -O2 -shared-intel -mcmodel=large -heap-arrays -traceback -fp-model strict

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
