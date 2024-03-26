# ------------------------------------------------------------
# Makefile for ECE420 W24 Lab4
# Edited by:
# Brandon Hoynick, CCID: hoynick
# Justin Javier, CCID: jjavier, ID: 1573754
# Saif Husnain, CCID: sahusnai, ID: 1573497
# ------------------------------------------------------------
# Makefile commands and Usage:
#
#   make main  			// Compiles main program to executable (using gcc and openmp).
#   make datatrim	  	// Compiles file parser program to executable (using gcc).
#   make zip      		// Create a '(target).zip' archive of '(allFiles)'.
#   make clean     		// Removes unneeded files produced in compilation (i.e. intermediate files), and main exe, zip created.
# ------------------------------------------------------------
# -std=c99 used to remove warnings about for-loop vars declared inside loop
# ------------------------------------------------------------

target = 	user28_lab4
allFiles = 	Makefile main.c Lab4_IO.c Lab4_IO.h timer.h


main: main.c Lab4_IO.c Lab4_IO.h timer.h
	mpicc -g -Wall -o main main.c Lab4_IO.c -lm -std=c99

datatrim: datatrim.c Lab4_IO.c Lab4_IO.h timer.h
	gcc -g -Wall -o datatrim datatrim.c Lab4_IO.c -lm

zip:
	zip $(target).zip $(allFiles)

clean:
	rm -f main
	rm -f $(target).zip
