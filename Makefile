CXX = g++
CFLAGS = -O2

LODEPNGSRC = lodepng/lodepng.cpp
LODEPNGOBJ = lodepng/lodepng.o

DIR_1 = 1-matrixless/
DIR_2 = 2-better-advection/
DIR_3 = 3-conjugate-gradients/
SRCF_1 = $(DIR_1)Fluid.cpp
SRCF_2 = $(DIR_2)Fluid.cpp
SRCF_3 = $(DIR_3)Fluid.cpp
OBJF_1 = $(subst .cpp,.o,$(SRCF_1))
OBJF_2 = $(subst .cpp,.o,$(SRCF_2))
OBJF_3 = $(subst .cpp,.o,$(SRCF_3))

all: matrixless better-advection conjugate-gradients

matrixless: $(OBJF_1) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_1)$@ $^ $(LDFLAGS) 
    
better-advection: $(OBJF_2) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_2)$@ $^ $(LDFLAGS) 
    
conjugate-gradients: $(OBJF_3) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_3)$@ $^ $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $^

lodepng:
	$(CXX) $(CFLAGS) -c -o $(LODEPNGOBJ) $(LODEPNGSRC)

clean:
	rm -f $(DIR_1)matrixless
	rm -f $(OBJF_1)
	rm -f $(DIR_2)better-advection
	rm -f $(OBJF_2)
	rm -f $(DIR_3)conjugate-gradients
	rm -f $(OBJF_3)
	rm -f lodepng/lodepng.o

.PHONY: clean all 
