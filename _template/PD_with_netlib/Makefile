# Project
BIN      = async

# Specify project files except BIN.
OBJS     =
HEADERS  =
FLAGS    = --std=c++11 -Wall -Wextra
INC      =

# Config _lib path
LIB      = ../../_lib
LIBA     = $(LIB)/netlib/netlib.a

# Compiler
CPP      = g++
RM       = rm -f
MKDIR    = mkdir -p

_LIB     = $(LIBA)
_INC     = $(INC) -I $(LIB)
_OBJS    = $(BIN).o $(OBJS)
_FLAGS   = $(_INC) $(FLAGS)

.PHONY: all clean clean-data clean-all

all: $(BIN)

clean-all: clean clean-data

clean:
	$(RM) $(BIN) $(_OBJS)

clean-data:
	$(RM) -r output

$(BIN): $(_OBJS) $(HEADERS)
	$(CPP) $(_OBJS) -o $(BIN) $(_LIB)

$(_OBJS): %.o: %.cc
	$(CPP) -c $(_FLAGS) $< -o $@

