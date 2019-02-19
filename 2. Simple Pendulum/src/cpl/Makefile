##################################################
# Author: József Stéger
# Date: 02 13 2017
# Usage (tested in linux):
#   make
# Summary:
#   1. compiles objects in ./build/
#   2. creates lib archive in .
##################################################

UTILS=datafit fft interpolate matrix odeint optimize random vector
BUILD=build

OBJ_LIBS=$(patsubst %,$(BUILD)/%.o,$(UTILS))
A_LIB=libcpl.a 

.PHONY: all distclean clean

all: $(A_LIB)

clean:
	rm -rf $(A_LIB) $(BUILD)

distclean: clean
	rm -rf $(PDIR)

$(OBJ_LIBS):
	@make $(BUILD)
	g++ -Wall -c -I. -O2 $(patsubst $(BUILD)/%.o,./%.cpp,$@) -o $@

$(A_LIB): $(OBJ_LIBS)
	ar -srcl $@ $?

$(BUILD):
	@mkdir -p $@

