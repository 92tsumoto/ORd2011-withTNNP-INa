#INCLUDES = syspara.h

bin_PROGRAMS = ORd_TP 
#ORd_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
ORd_TP_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -vec-report2 -pg

ORd_TP_INCLUDE = -I./include
ORd_TP_LDADD = -L/HOME/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o
#ORd_LDADD = /home/tsumoto/lib/xhplot.o

ORd_TP_SOURCES = ./src/ord.c ./src/val_constants.c ./src/mem.c\
	./src/win_out.c ./src/eular.c ./src/input_para.c ./src/system.c ./src/currents.c\
	./src/eventloop.c ./src/exptable.c ./src/syspara.h

