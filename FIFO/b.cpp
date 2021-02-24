//example from
//https://www.geeksforgeeks.org/named-pipe-fifo-example-c-program/

// C program to implement one side of FIFO 
// This side reads first, then reads 
#include <stdio.h> 
#include <fcntl.h> 
#include <sys/stat.h> 
#include <sys/types.h> 
#include <unistd.h> 
  
int main(){
 
	int fd1; 
 
	double data[2];
 
	// FIFO file path 
	const char *myfifo = "myfifo"; 
	const char *fifoCheck = "fifoCheck"; 
  
	// ****************************************
	// Creating the named file(FIFO) 
	// mkfifo(<pathname>,<permission>) 
	mkfifo(myfifo, 0666); //path, permission mode  
	mkfifo(fifoCheck, 0666); //path, permission mode  
  
        // First open in read only and read 
        fd1 = open(myfifo,O_RDONLY); 
        read(fd1, &data, sizeof(double)*2); 
  
        // Print the read string and close 
        printf("User1: %g %g\n", data[0], data[1]); 
        close(fd1);
	
	// **************************************** 
	printf("wait\n");
	sleep(0.5);
	printf("Send OK\n");
  

	int OK = 100;
        // Now open in write mode and write 
        // string taken from user. 
        fd1 = open(fifoCheck,O_WRONLY); 
        write(fd1, &OK, sizeof(int)); 
        close(fd1); 
    return 0; 
} 

