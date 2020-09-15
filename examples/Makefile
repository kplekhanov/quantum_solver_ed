#--- main directory of quantum_solver_ed

MAIN_DIR = ..

#--- binaries to make

BIN = test test_sym test_spin test_read

#--- geting all the important variables

include $(MAIN_DIR)/options.mk

SRC = $(BIN:%=%.cpp)
OBJ = $(BIN:%=%.o)
DEP = $(BIN:%=%.d)

OBJ_MAIN = $(wildcard $(MAIN_DIR)/$(BUILD_DIR)/*.o)
DEP_MAIN = $(OBJ_MAIN:%.o=%.d)

#--- rules

$(BIN): %: %.o $(OBJ_MAIN)
	$(CCCOM) -o $@ $^ $(CCFLAGS) $(LDFLAGS)

%.o: %.cpp
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.PHONY: all clean

all: $(BIN)

clean:
	rm -rf $(BIN) $(OBJ) $(DEP)

#--- dependency files

-include $(DEP)
-include $(DEP_MAIN)
