// example from
//https://www.geeksforgeeks.org/named-pipe-fifo-example-c-program/

// C program to implement one side of FIFO 
// This side writes first, then reads 
#include <stdio.h> 
#include <fcntl.h> 
#include <sys/stat.h> 
#include <sys/types.h> 
#include <unistd.h> 
  
int main() 
{ 
	int fd; 
	int fd1;

	double data[2];
	data[0] = 5.0;
	data[1] = 5.0;

	int OK[1];
	OK[0] = 0;

	printf("%g %g\n", data[0], data[1]);

	// FIFO file path 
	const char *myfifo = "myfifo"; 
	const char *fifoCheck = "fifoCheck"; 

	// Creating the named file(FIFO) 
	// mkfifo(<pathname>, <permission>) 
	mkfifo(myfifo, 0666); //path, permission mode 
	mkfifo(fifoCheck, 0666); //path, permission mode 
  
        // Open FIFO for write only 
        fd = open(myfifo, O_WRONLY); 
  
        write(fd, &data, sizeof(double) * 2); 
        close(fd); 
  
  
        // Read from FIFO 
        // Open FIFO for Read only 
        fd1 = open(fifoCheck, O_RDONLY); 
        read(fd1, &OK, sizeof(int)); 
  
        // Print the read message 
        printf("User2: %d\n", OK[0]); 
        close(fd1); 
 
   return 0; 
} 

