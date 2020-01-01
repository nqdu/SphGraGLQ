makemat=bin/makemat

# compiler and options
cc = gcc
CFLAGS= -lm -fopenmp
includepath=include/

# c scripts
deps = $(shell find include/*.h)
src = $(shell find src/*.c)

# objective
obj = $(src:%.c=%.o) 

# compile
$(makemat): $(obj)
	$(cc) -o $(makemat) $(obj) $(CFLAGS)
	
%.o: %.c $(deps)
	$(cc) -g -c -I$(includepath) $< -o $@ -fopenmp

clean:
	rm  src/*.o 
	
