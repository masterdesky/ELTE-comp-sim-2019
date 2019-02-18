##################################################
# Author: József Stéger
# Date: 14 02 2017
# Usage (tested in linux):
#   make
# Summary:
#   1. compiles objects in build/
#   2. links binary in ./
##################################################

APPS=pendulum pendulum-gl
LIB=-L../cpl -lcpl
BUILD=build
INCLUDE=-I../cpl
LIB_GL=-lGL -lGLU -lglut

OBJ_APPS=$(patsubst %,$(BUILD)/%.o,$(APPS))
A_APPS=$(patsubst %,%.bin,$(APPS))

.PHONY: all clean

all: $(OBJ_APPS) $(A_APPS)

clean:
	rm -rf $(A_APPS) $(BUILD) pendulum.data

$(BUILD)/%.o: %.cpp
	@make $(BUILD)
	g++ -Wall -c $(INCLUDE) -O2 $< -o $@

$(BUILD):
	@mkdir -p $@

%.bin: $(BUILD)/%.o
	g++ -Wall $< $(LIB) -o $@

pendulum-gl.bin: $(BUILD)/pendulum-gl.o
	g++ -Wall $< $(LIB) $(LIB_GL) -o $@

