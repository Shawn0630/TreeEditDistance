#############################################################################
#
#    Tree Edit Distance
#
#    (C) Copyright 2010 Qichan Ma
#    Modified by Shaofeng Jiang
#

# This Makefile uses powerful features of GNU Make.
# On Linux platforms, just follow the usage to compile the program.
# On other platforms which have installed GNU Make, replace "make" with "gmake" in the following usage to compile the program.
#
# Usage: 
# ------
#   $ make            (=make all) compile and link
#   $ make objs       compile only (no linking)
#   $ make clean      clean objects and dependencies
#   $ make distclean  clean objects, dependencies and the executable
#   $ make rebuild    rebuild the program. The same as make distclean & make all.
#   $ make help       get the usage of the makefile
#===========================================================================


## Customizable Section
##==========================================================================

# The executable file name.
# If not specified, current directory name or `a.out' will be used.
PROGRAM      = TreeEditDistance

# The source files directories (do not include last '/').
# If not specified, only the current directory will be searched.
SRC_DIRS     = src test

# Header file directories (do not include last '/').
INCLUDE_DIRS = 

# The library directories (do not include last '/').
LIB_DIRS     = 

# The library names.
LIB_NAMES    = 

# The static 

# Object files directory (do not include last '/').
OBJ_DIR      = obj

# Output file directory (do not include last '/').
OUTPUT_DIR   = out

# The source file types (headers excluded).
# .c indicates C source files, and others indicate C++ ones.
SRCEXTS      = .c .C .cc .cpp .CPP .c++ .cxx .cp

vpath %.c   $(SRC_DIRS)
vpath %.C   $(SRC_DIRS)       
vpath %.cc  $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)
vpath %.CPP $(SRC_DIRS)
vpath %.c++ $(SRC_DIRS)
vpath %.cxx $(SRC_DIRS)
vpath %.cp  $(SRC_DIRS)
vpath %.o   $(OBJ_DIR)
vpath %.d   $(OBJ_DIR)

# The C program compiler.
CC           = gcc

# The C++ program compiler.
CXX          = g++

# The pre-processor and compiler options for C programs.
CFLAGS       = -g3 -O2 -c -g -std=c++11 #-pg 

# The pre-processor and compiler options for C++ programs.
CXXFLAGS     = $(CFLAGS) 

# The pre-processor and compiler options for C++ Gtest programs.
GTESTLIBS = -lgtest -lgtest_main -lpthread

# The linker options for C programs.
CLDFLAGS     = -g3 -O2 #-pg

# The linker options for C++ programs.
CXXLDFLAGS   = $(CLDFLAGS)

# The command used to delete file.
RM           = rm -f

## Stable Section
##==========================================================================

SHELL   = /bin/sh
EMPTY   =
SPACE   = $(EMPTY) $(EMPTY)
ifeq ($(PROGRAM),)
  CUR_PATH_NAMES = $(subst /,$(SPACE),$(subst $(SPACE),_,$(CURDIR)))
  PROGRAM = $(word $(words $(CUR_PATH_NAMES)),$(CUR_PATH_NAMES))
  ifeq ($(PROGRAM),)
    PROGRAM = a.out
  endif
endif
ifeq ($(SRC_DIRS),)
  SRC_DIRS = .
endif
ifeq ($(INCLUDE_DIRS),)
  INCLUDE_DIRS = .
endif
ifeq ($(OBJ_DIR),)
  OBJ_DIR = .
endif
ifeq ($(OUTPUT_DIR),)
  OUTPUT_DIR = .
endif


SOURCES   = $(foreach d,$(SRC_DIRS),$(wildcard $(addprefix $(d)/*,$(SRCEXTS))))
SRC_CXX   = $(filter-out %.c,$(SOURCES))
OBJS      = $(addprefix $(OBJ_DIR)/,$(notdir $(foreach x,$(SRCEXTS),$(patsubst %$(x),%.o,$(filter %$(x),$(SOURCES))))))
DEPS      = $(OBJS:.o=.d)

## Define some useful variables.
INCLUDES    = $(addprefix -I ,$(INCLUDE_DIRS)) $(addprefix -I ,$(SRC_DIRS))
LIBS        = $(addprefix -L ,$(LIB_DIRS)) $(addprefix -l ,$(LIB_NAMES))
OUT         = $(OUTPUT_DIR)/$(PROGRAM)
DEP_OPT     = -MMD -MP -MT"$(@:%.o=%.d)"
COMPILE.c   = $(CC)  $(DEP_OPT) $(CFLAGS)
COMPILE.cxx = $(CXX) $(DEP_OPT) $(CXXFLAGS)
LINK.c      = $(CC)  $(CLDFLAGS)
LINK.cxx    = $(CXX) $(CXXLDFLAGS)

.PHONY: all objs clean cleantest distclean rebuild help check

# Delete the default suffixes
.SUFFIXES:

all: $(OUT)

# Rules for generating object files (.o) and dependency files (.d).
#----------------------------------------
objs:$(OBJS)

$(OBJ_DIR)/%.o:%.c
	@echo 'Building target: $@'
	$(COMPILE.c) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.C
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.cc
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.cpp
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.CPP
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.c++
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.cp
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

$(OBJ_DIR)/%.o:%.cxx
	@echo 'Building target: $@'
	$(COMPILE.cxx) $< $(INCLUDES) -o $@
	@echo 'Finished building: $@'
	@echo ''

# Rules for generating the executable.
#-------------------------------------
$(OUT):$(OBJS)
ifeq ($(SRC_CXX),)              # C program
	@echo 'Building target: $@'
	$(LINK.c)   $(OBJS) $(LIBS) -o $@ -lgtest -lpthread
	@echo 'Finished building: $@'
	@echo ''
	@echo 'Type "cd $(OUTPUT_DIR)" to change to output file directory.'
	@echo 'Then type "./$(PROGRAM)" to execute the program.'
	@echo ''
#	@echo 'Type ./$@ to execute the program.'
else                            # C++ program
	@echo 'Building target: $@'
	$(LINK.cxx) $(OBJS) $(LIBS) -o $@ -lgtest -lpthread
	@echo 'Finished building: $@'
	@echo ''
	@echo 'Type "cd $(OUTPUT_DIR)" to change to output file directory.'
	@echo 'Then type "./$(PROGRAM)" to execute the program.'
	@echo ''
#	@echo 'Type ./$@ to execute the program.'
endif

ifneq ($(MAKECMDGOALS),clean)
  ifneq ($(MAKECMDGOALS),distclean)
    ifneq ($(MAKECMDGOALS),rebuild)
      ifneq ($(strip $(DEPS)),)
        -include $(DEPS)
      endif
    endif
  endif
endif


# clean objects and dependencies
clean:
	$(RM) $(OBJS) $(DEPS)

cleantest:
	$(RM) $(OBJS) $(DEPS)

# clean objects, dependencies and the executable
distclean: clean
	$(RM) $(OUT)

rebuild: distclean all
	
# Show help.
help:
	@echo 'This Makefile uses powerful features of GNU Make.'
	@echo 'On Linux platforms, just follow the usage to compile the program.'
	@echo 'On other platforms which have installed GNU Make, replace "make" with "gmake" in the following usage to compile the program.'
	@echo 'Usage: make [TARGET]'
	@echo 'TARGETS:'
	@echo '  all        (=make) compile and link.'
	@echo '  objs       compile only (no linking).'
	@echo '  clean      clean objects and dependencies.'
	@echo '  distclean  clean objects, dependencies and the executable.'
	@echo '  rebuild    rebuild the program. The same as make distclean & make all.'
	@echo '  help       print this message.'


.PHONY: show
# Show variables (for debug use only.)
show:
	@echo 'PROGRAM     :' $(PROGRAM)
	@echo 'SRC_DIRS    :' $(SRC_DIRS)
	@echo 'SOURCES     :' $(SOURCES)
	@echo 'TESTFILES   :' $(TESTFILES)
	@echo 'MAINFILES   :' $(MAINFILES)
	@echo 'SRC_CXX     :' $(SRC_CXX)
	@echo 'OBJS        :' $(OBJS)
	@echo 'TESTOBJS    :' $(TESTOBJS)
	@echo 'MAINOBJS    :' $(MAINOBJS)
	@echo 'DEPS        :' $(DEPS)
	@echo 'INCLUDES    :' $(INCLUDES)
	@echo 'LIBS        :' $(LIBS)
	@echo 'OUT         :' $(OUT)
	@echo 'DEP_OPT     :' $(DEP_OPT)
	@echo 'COMPILE.c   :' $(COMPILE.c)
	@echo 'COMPILE.cxx :' $(COMPILE.cxx)
	@echo 'link.c      :' $(LINK.c)
	@echo 'link.cxx    :' $(LINK.cxx)
	
## End of the Makefile ##
#############################################################################
