# Makefile to compile new FLIP
# DATE 2010 09 28
#
FORTRAN_COMPILER = gfortran
PACKAGE	= test_flip.exe
#
SRCS	= CTIPE-int.FOR \
	CMINOR.mod.FOR \
	ELECXS.FOR \
	INIT-PROFILES.FOR \
	KEMPRN.FOR \
	MINORA.FOR \
	Neut_Heating.FOR \
	Photoel-Freqs.mod.FOR \
	RSDENA_EXB.FOR \
	RSJACA.FOR \
	RSLPSD.FOR \
	RSLPST.FOR \
	RSPE2B.FOR \
	RSPRIM.mod.FOR \
	RSTEMC_EXB.FOR \
	Rates.mod.FOR \
	Rsmnsd.FOR
#
#HEADS	= $(PACKAGE).h
OBJS	= $(SRCS:.FOR=.o)

#
FILES	= README Makefile $(HEADS) $(SRCS)
VER	= `date +%Y%m%d`


### command and flags ###
# uncomment when debugging
#DEBUG	= -ggdb -pg # -lefence

# common (*.o)
LD	= $(FORTRAN_COMPILER)
LDFLAGS	= $(DEBUG)
LDLIBS	= -lm

# C (*.c)
CC	= gcc
CFLAGS	= -g -O2 -Wall $(DEBUG)
CPPFLAGS= -I.

# C++ (*.cc)
CXX	= g++
CXXFLAGS= -g -O2 -Wall $(DEBUG)

# Fortran77 (*.f .for .FOR)
FC	= $(FORTRAN_COMPILER)
FFLAGS	= $(DEBUG)

# Fortran90 (*.f90)
FC90	= $(FORTRAN_COMPILER)
FFLAGS90= $(DEBUG)

# Pascal (*.p)
PC	= pc
PFLAGS	= -Wall $(DEBUG)

# etc
SHELL	= /bin/sh
RM	= rm -f
PROF	= gprof


### rules ###

.SUFFIXES:
.SUFFIXES: .o .c .cc .f .FOR .for .f90 .mod .p

all: $(PACKAGE)

$(PACKAGE): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $@ $(LDLIBS)

$(OBJS): $(HEADS) Makefile

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@
.cc.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
.f.o:
	$(FC) $(FFLAGS) -c $< -o $@
.for.o:
	$(FC) $(FFLAGS) -c $< -o $@
.FOR.o:
	$(FC) $(FFLAGS) -c $< -o $@
.f90.mod:
	$(FC90) $(FFLAGS90) -c $< -o $@
.f90.o:
	$(FC90) $(FFLAGS90) -c $< -o $@
.p.o:
	$(PC) $(PFLAGS) $(CPPFLAGS) -c $< -o $@


### useful commands ###

clean:
	$(RM) $(PACKAGE) $(OBJS)
	$(RM) core *.mod

tar:
	@echo $(PACKAGE)-$(VER) > .package
	@$(RM) -r `cat .package`
	@mkdir `cat .package`
	@ln $(FILES) `cat .package`
	tar cvf - `cat .package` | gzip -9 > `cat .package`.tar.gz
	@$(RM) -r `cat .package` .package

zip:
	zip -9 $(PACKAGE)-$(VER).zip $(FILES)


prof: run
	$(PROF) $(PACKAGE) | less

run: all
	./$(PACKAGE) < sample-data | less
