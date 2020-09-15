#--- flags and options; to change if needed, especially the armadillo one

include options.mk

#--- file variables

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))
DEP = $(OBK:%.o=%.d)

#--- rules

$(OBJ): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(dir $@)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.PHONY: all clean test

all: $(OBJ)

clean:
	rm -rf build

test:
	echo $(SRC_DIR)
	echo $(SRC)

#--- dependency files

-include $(DEP)
