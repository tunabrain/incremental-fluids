CXX = g++
CFLAGS = -O2

LODEPNGSRC = lodepng/lodepng.cpp
LODEPNGOBJ = lodepng/lodepng.o

DIR_1 = 1-matrixless/
DIR_2 = 2-better-advection/
DIR_3 = 3-conjugate-gradients/
DIR_4 = 4-solid-boundaries/
DIR_5 = 5-curved-boundaries/
DIR_6 = 6-heat/
DIR_7 = 7-variable-density/
DIR_8 = 8-flip/
SRCF_1 = $(DIR_1)Fluid.cpp
SRCF_2 = $(DIR_2)Fluid.cpp
SRCF_3 = $(DIR_3)Fluid.cpp
SRCF_4 = $(DIR_4)Fluid.cpp
SRCF_5 = $(DIR_5)Fluid.cpp
SRCF_6 = $(DIR_6)Fluid.cpp
SRCF_7 = $(DIR_7)Fluid.cpp
SRCF_8 = $(DIR_8)Fluid.cpp
OBJF_1 = $(subst .cpp,.o,$(SRCF_1))
OBJF_2 = $(subst .cpp,.o,$(SRCF_2))
OBJF_3 = $(subst .cpp,.o,$(SRCF_3))
OBJF_4 = $(subst .cpp,.o,$(SRCF_4))
OBJF_5 = $(subst .cpp,.o,$(SRCF_5))
OBJF_6 = $(subst .cpp,.o,$(SRCF_6))
OBJF_7 = $(subst .cpp,.o,$(SRCF_7))
OBJF_8 = $(subst .cpp,.o,$(SRCF_8))

all: matrixless better-advection conjugate-gradients solid-boundaries curved-boundaries heat variable-density flip

matrixless: $(OBJF_1) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_1)$@ $^ $(LDFLAGS) 
    
better-advection: $(OBJF_2) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_2)$@ $^ $(LDFLAGS) 
    
conjugate-gradients: $(OBJF_3) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_3)$@ $^ $(LDFLAGS) 
    
solid-boundaries: $(OBJF_4) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_4)$@ $^ $(LDFLAGS) 
    
curved-boundaries: $(OBJF_5) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_5)$@ $^ $(LDFLAGS) 
    
heat: $(OBJF_6) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_6)$@ $^ $(LDFLAGS) 
    
variable-density: $(OBJF_7) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_7)$@ $^ $(LDFLAGS) 
    
flip: $(OBJF_8) $(LODEPNGOBJ)
	$(CXX) -o $(DIR_8)$@ $^ $(LDFLAGS) 

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
	rm -f $(DIR_4)solid-boundaries
	rm -f $(OBJF_4)
	rm -f $(DIR_5)curved-boundaries
	rm -f $(OBJF_5)
	rm -f $(DIR_6)heat
	rm -f $(OBJF_6)
	rm -f $(DIR_7)variable-density
	rm -f $(OBJF_7)
	rm -f $(DIR_8)flip
	rm -f $(OBJF_8)
	rm -f lodepng/lodepng.o

.PHONY: clean all 
