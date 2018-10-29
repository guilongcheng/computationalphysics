CC = g++
Cflag =

obj=#lib_linalg.o

all:$(obj) main.cpp
	$(CC) main.cpp $(obj) -o test


$(obj):%.o:%.cpp
	$(CC) -c $(Cflag) $< -o $@

clean:
	rm -rf $(obj) test




