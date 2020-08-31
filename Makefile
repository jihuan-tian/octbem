# List of Octave scripts
OCTAVE_FILES := $(wildcard *.m)
# List of Matlab scripts
MATLAB_FILES := $(subst .m,.mlb,$(OCTAVE_FILES))

# List of source oct files written in C++
SOURCE_OCT_FILES := GaussLegendreRule.cpp GaussJacobiRule.cpp
# List of object files from which oct files will be generated
OCT_OBJ_FILES := $(subst .cpp,.o,$(SOURCE_OCT_FILES))
# List of oct files compiled from C++ source code
OCT_FILES := $(subst .cpp,.oct,$(SOURCE_OCT_FILES))

# CPP compiler options
DEFS :=
CXXFLAGS := -std=gnu++11 -O2 -g $(DEFS)

# ##################################
# Definition of Makefile targets
# ##################################
all: matlab octs doc
matlab: $(MATLAB_FILES)

# ##################################
# Generate GNU Octave oct files
# ##################################
octs: $(OCT_FILES)
GaussLegendreRule.oct : GaussLegendreRule.cpp quadrule.o
	mkoctfile $(DEFS) $^

GaussJacobiRule.oct : GaussJacobiRule.cpp quadrule.o
	mkoctfile $(DEFS) $^

# Generate documentation.
doc:
	cd doc
	doxygen Doxyfile

# Clean object files from compiling C++ source code
clean_objs:
	rm -f $(OCT_OBJ_FILES) quadrule.o

# Clean GNU Octave oct files
clean_octs:
	rm -f $(OCT_FILES)

# Clean automatically generated GNU Octave files.
clean_mats:
	rm $(MATLAB_FILES)

# Clean all generated files
clean: clean_octs clean_objs clean_checks
	rm -f *~

%.o : %.cpp
	g++ $(CXXFLAGS) -c $< -fPIC

%.mlb: %.m
	oct2mlb.pl $< $@
