CXXFLAGS := -O3  -std=c++11 -mavx #-g 

.PHONY: all

SOURCES = $(wildcard source/*.cpp)
OBJS = $(patsubst %.cpp,%.o,$(SOURCES))

SOURCES1 = $(wildcard source/*.c)
OBJS1 = $(patsubst %.c,%.o,$(SOURCES1))

all: $(OBJS) $(OBJS1)
	ar -ruv lib/libGraph.a $(OBJS) $(OBJS1)
	rm -f $(OBJS) $(OBJS1)
	g++ $(CXXFLAGS) -o bin/hb_network_hydration.x hb_network_hydration.cpp -I"include" -L"lib" -lGraph -I"../libgmxcpp-5.0.2/include/gmxcpp" -L"../libgmxcpp-5.0.2/build/lib" -lgmxcpp -I"../libxdrfile-master/include" -L"../libxdrfile-master/build/lib" -lxdrfile -L"../gromacs-2016.4/build/lib" -lgromacs
$(OBJS): %.o: %.cpp
	g++ $(CXXFLAGS) -c $< -o $@ -I"include"
	
$(OBJS1): %.o: %.c
	gcc $(CXXFLAGS) -c $< -o $@ -I"include"

clean:
	rm -f $(OBJS) $(OBJS1) *.o lib/*.a bin/ksp
	

	
