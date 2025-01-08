#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lcblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native
