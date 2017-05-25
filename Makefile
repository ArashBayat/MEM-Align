

all:
	rm -f *.o
	gcc 	-c -DTICKPROF 	-m64 -msse4.2 -O3 ssw.c -o ssw.o
	g++ 	-c -DTICKPROF 	-m64 -msse4.2 -O3 main.cpp -o main.o
	g++ 	   -DTICKPROF 	-m64 -msse4.2 -O3 main.o ssw.o -o MA.o
	gcc 	-c 				-m64 -msse4.2 -O3 ssw.c -o ssw_f.o
	g++ 	-c 				-m64 -msse4.2 -O3 main.cpp -o main_f.o
	g++ 	   				-m64 -msse4.2 -O3 main_f.o ssw_f.o -o MA_f.o
	gcc -g3 -c 				-m64 -msse4.2 ssw.c -o ssw_g.o
	g++ -g3 -c 				-m64 -msse4.2 main.cpp -o main_g.o
	g++ -g3    				-m64 -msse4.2 main_g.o ssw_g.o -o MA_g.o			
	gcc -pg -c 				-m64 -msse4.2 ssw.c -o ssw_pg.o
	g++ -pg -c 				-m64 -msse4.2 main.cpp -o main_pg.o
	g++ -pg    				-m64 -msse4.2 main_pg.o ssw_pg.o -o MA_pg.o

clean:
	rm -f *.o
