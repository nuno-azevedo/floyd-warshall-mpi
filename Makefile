CC = mpicc
CCFLAGS = -ffast-math -Ofast -lm

OUT = floyd
SRC = src/floyd.c
HOSTFILE = hosts.txt

default: build

build:
	${CC} ${CCFLAGS} ${SRC} -o ${OUT}

clean:
	rm -f ${OUT} floyd-warshall-mpi.zip

archive: build
	zip -r floyd-warshall-mpi.zip inputs report scripts src times hosts.txt Makefile README.md

tests: build
	mpirun -np 1  -hostfile ${HOSTFILE} ${OUT} < inputs/input5x5.txt | diff - inputs/output5x5.txt

	mpirun -np 1  -hostfile ${HOSTFILE} ${OUT} < inputs/input6x6.txt | diff - inputs/output6x6.txt
	mpirun -np 4  -hostfile ${HOSTFILE} ${OUT} < inputs/input6x6.txt | diff - inputs/output6x6.txt
	@# mpirun -np 9  -hostfile ${HOSTFILE} ${OUT} < inputs/input6x6.txt | diff - inputs/output6x6.txt

	mpirun -np 1  -hostfile ${HOSTFILE} ${OUT} < inputs/input12x12.txt | diff - inputs/output12x12.txt
	mpirun -np 4  -hostfile ${HOSTFILE} ${OUT} < inputs/input12x12.txt | diff - inputs/output12x12.txt
	@# mpirun -np 9  -hostfile ${HOSTFILE} ${OUT} < inputs/input12x12.txt | diff - inputs/output12x12.txt
	@# mpirun -np 16 -hostfile ${HOSTFILE} ${OUT} < inputs/input12x12.txt | diff - inputs/output12x12.txt

	mpirun -np 1  -hostfile ${HOSTFILE} ${OUT} < inputs/input60x60.txt | diff - inputs/output60x60.txt
	mpirun -np 4  -hostfile ${HOSTFILE} ${OUT} < inputs/input60x60.txt | diff - inputs/output60x60.txt
	@# mpirun -np 9  -hostfile ${HOSTFILE} ${OUT} < inputs/input60x60.txt | diff - inputs/output60x60.txt
	@# mpirun -np 16 -hostfile ${HOSTFILE} ${OUT} < inputs/input60x60.txt | diff - inputs/output60x60.txt

	mpirun -np 1  -hostfile ${HOSTFILE} ${OUT} < inputs/input300x300.txt | diff - inputs/output300x300.txt
	mpirun -np 4  -hostfile ${HOSTFILE} ${OUT} < inputs/input300x300.txt | diff - inputs/output300x300.txt
	@# mpirun -np 9  -hostfile ${HOSTFILE} ${OUT} < inputs/input300x300.txt | diff - inputs/output300x300.txt
	@# mpirun -np 16 -hostfile ${HOSTFILE} ${OUT} < inputs/input300x300.txt | diff - inputs/output300x300.txt

	@# mpirun -np 16 -hostfile ${HOSTFILE} ${OUT} < inputs/input600x600.txt | diff - inputs/output600x600.txt
	@# mpirun -np 25 -hostfile ${HOSTFILE} ${OUT} < inputs/input600x600.txt | diff - inputs/output600x600.txt
	@# mpirun -np 36 -hostfile ${HOSTFILE} ${OUT} < inputs/input600x600.txt | diff - inputs/output600x600.txt

	@# mpirun -np 64 -hostfile ${HOSTFILE} ${OUT} < inputs/input1000x1000.txt | diff - inputs/output1000x1000.txt

	@# mpirun -np 64 -hostfile ${HOSTFILE} ${OUT} < inputs/input1200x1200.txt | diff - inputs/output1200x1200.txt

	@# mpirun -np 64 -hostfile ${HOSTFILE} ${OUT} < inputs/input1400x1400.txt | diff - inputs/output1400x1400.txt

	@echo "TESTS PASSED. ALL DONE."
