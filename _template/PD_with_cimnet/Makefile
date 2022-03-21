# Project
BIN      = async

# Specify project files except BIN.
OBJS     =
HEADERS  =
FLAGS    = --std=c++11 -Wall -Wextra -O3
INC      =

# Config _lib path
LIB      = ../../_lib

# Compiler
CPP      = g++
RM       = rm -f
MKDIR    = mkdir -p

_INC     = $(INC) -I $(LIB)
_OBJS    = $(BIN).o $(OBJS)
_FLAGS   = $(_INC) $(FLAGS)

.PHONY: all clean clean-data clean-all

all: $(BIN).out

clean-all: clean clean-data

clean:
	$(RM) $(BIN).out $(_OBJS)

clean-data:
	$(RM) -r output

$(BIN).out: $(_OBJS) $(HEADERS)
	$(CPP) $(_OBJS) -o $(BIN).out

$(_OBJS): %.o: %.cc
	$(CPP) -c $(_FLAGS) $< -o $@

