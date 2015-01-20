# Makefile
#
CC      = gcc#
CFLAGS	=-Wall -pedantic#

rm	= rm -f#				delete file

all:	_
clean::		;@ $(MAKE) T='$T' _clean
_clean:	_	;  $(rm) *.o $T a.out core *.tmp *.ps *.bak
run::	_
_:		;@ echo -------------------- $D --------------------


D =	Aggregated Hidden Markov Models

MCMC= mp_parameter.o mp_mcmc.o likelihood.o mp_proposal.o

HMM = matrixIO.o matrixExp.o matrixArith.o mctools.o mp_HMMlikelihoods.o mp_HMMutils.o nw_data.o

T = $(HMM) $(MCMC) mp_OneOpen$x

all:	$T

INCLUDE=.
GSLINC = `gsl-config --cflags`

GSL = `gsl-config --cflags` `gsl-config --libs` -lm


mp_parameter.o:	mp_parameter.c;
		$(CC) $(CFLAGS) -c -DPROP -o $@ mp_parameter.c -I$(INCLUDE) $(GSLINC) 

mp_proposal.o:	mp_proposal.c;
		$(CC) $(CFLAGS) -c -o $@ mp_proposal.c -I$(INCLUDE) $(GSLINC) 

mp_HMMlikelihoods.o:	mp_HMMlikelihoods.c;
		$(CC) $(CFLAGS) -c -o $@ mp_HMMlikelihoods.c -I$(INCLUDE) $(GSLINC)

mp_HMMutils.o:	mp_HMMutils.c;
		$(CC) $(CFLAGS) -c -o $@ mp_HMMutils.c -I$(INCLUDE) $(GSLINC)

mctools.o:	mctools.c;
	$(CC) $(CFLAGS) -c -o $@ mctools.c -I$(INCLUDE) $(GSLINC)

matrixArith.o:	matrixArith.c;
	$(CC) $(CFLAGS) -c -o $@ matrixArith.c -I$(INCLUDE) $(GSLINC)

matrixExp.o:	matrixExp.c;
		$(CC) $(CFLAGS) -c -o $@ matrixExp.c -I$(INCLUDE) $(GSLINC)

mp_OneOpen$x:	$(MCMC) $(HMM) mp_OneOpen.c
		$(CC) $(CFLAGS) -o $@ $(MCMC) $(HMM) mp_OneOpen.c -I$(INCLUDE) $(GSL) 

mp_mcmc.o:	mp_mcmc.c;
		$(CC) $(CFLAGS) -c -o $@ mp_mcmc.c -I$(INCLUDE) $(GSLINC)

likelihood.o:	likelihood.c
		$(CC) $(CFLAGS) -c -o $@ likelihood.c -I$(INCLUDE) $(GSLINC)

matrixIO.o:	matrixIO.c;
		$(CC) $(CFLAGS) -c -o $@ matrixIO.c -I$(INCLUDE) $(GSLINC)