CC = mpic++
FLAGS = -O3 -std=c++11 -pedantic -g -DUSE_MPI=1 \
        -Werror=return-type -Werror=uninitialized -Wall

LYRA_DIR = /home/stefan/projects/Lyra/include/
MPI_DIR = /usr/lib/openmpi/
MPI_INC = /usr/include/openmpi/

LNK = -L./ -L$(MPI_DIR) -lmpi -lpng
INC = -I./ -I$(MPI_INC) -I$(XOSHIRO_DIR) -I$(LYRA_DIR)

COMP = $(CC) $(FLAGS) $(INC)
LINK = $(CC) $(FLAGS) $(INC) $(LNK)

EXE = mpi_xy
EXT = cpp
SRC = $(wildcard *.$(EXT))

# For windows:
#MAKE_DIR = $(if exist $(1),,mkdir $(1))
#S=\\
# Linux and Unix-like:
MAKE_DIR = mkdir -p $(1)
S=/



OBJ_DIR = obj
OBJ = $(SRC:%.$(EXT)=$(OBJ_DIR)$(S)%.o)
OBJ_DIRS = $(dir $(OBJ))
DEPS = $(OBJ:%.o=%.d)

.PHONY: dirs all help clean

all : dirs $(EXE)

dirs : $(OBJ_DIR)

$(OBJ_DIR) :
	$(call $(MAKE_DIR),$@)

help :
	@echo "SRC is $(SRC)"
	@echo "OBJ is $(OBJ)"
	@echo "DEPS is $(DEPS)"

$(EXE) : $(OBJ)
	$(LINK) $(OBJ) -o $@

$(OBJ_DIR)$(S)%.o : %.$(EXT)
	$(call MAKE_DIR,$(dir $@))
	$(COMP) -c $< -o $@
	$(COMP) -M -MT '$@' $< -MF $(@:%.o=%.d)

clean:
	rm -r $(OBJ_DIR)
	rm -f $(EXE)

-include $(DEPS)
