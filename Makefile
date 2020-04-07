CC = gcc
#C compiler
CFLAGS = -O3 -pg -std=c99 -o
#C compiler options
SRC = ./src
#Directory of sources
INCLUDE = -I $(SRC)/include
#Inclued folder
LIBS = -lm
RM = rm -vf
#Delete command


HEAD = finite_volume meshing file_io Riemann_solver tools
#Name of header files or subdirectories
SOURCE = hydrocode
#Name of the main source 

all : modules
	@echo "**********Generate executable file***********"
	$(CC) $(CFLAGS) $(SRC)/$(SOURCE).out $(SRC)/$(SOURCE).c $(addsuffix /*.a,$(addprefix $(SRC)/,$(HEAD))) $(INCLUDE) $(LIBS)
.PHONYP:all

modules:
#Enter each subdirectory
#Call the Makefile in the subdirectory
	@for n in $(HEAD); do \
	( make -f Makefile.sub --directory=$(SRC)/$$n )  \
	done;
.PHONYP:modules

clean:
#Clean in the subdirectory
	@$(RM) $(SRC)/$(SOURCE).o
	@for n in $(HEAD); do \
	( make -f Makefile.sub --directory=$(SRC)/$$n clean ) \
	done;
.PHONYP:clean
