CXX = g++
CFLAGS = -O2

LODEPNGSRC = lodepng/lodepng.cpp
LODEPNGOBJ = lodepng/lodepng.o

DIR_1 = 1-matrixless/
SRCF_1 = $(DIR_1)Fluid.cpp
OBJF_1 = $(subst .cpp,.o,$(SRCF_1))
DIR_2 = 2-better-advection/
SRCF_2 = $(DIR_2)Fluid.cpp
OBJF_2 = $(subst .cpp,.o,$(SRCF_2))

all: matrixless better-advection

matrixless: $(OBJF_1) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_1)$@ $^ $(LDFLAGS) 
    
better-advection: $(OBJF_2) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_2)$@ $^ $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^

lodepng:
	$(CXX) $(CFLAGS) -c -o $(LODEPNGOBJ) $(LODEPNGSRC)

clean:
	rm -f $(DIR_1)matrixless
	rm -f $(OBJF_1)
	rm -f $(DIR_2)better-advection
	rm -f $(OBJF_2)
	rm -f lodepng/lodepng.o

.PHONY: clean all 
