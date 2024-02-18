#include <stdio.h>
#include "mpi.h"
#include <time.h>


#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255

typedef struct C{
double real;
double imag;

}complex;

void save_pgm(const char *filename, int image[HEIGHT][WIDTH]) {
    FILE* pgmimg; 
    int temp;
    pgmimg = fopen(filename, "wb"); 
    fprintf(pgmimg, "P2\n"); // Writing Magic Number to the File   
    fprintf(pgmimg, "%d %d\n", WIDTH, HEIGHT);  // Writing Width and Height
    fprintf(pgmimg, "255\n");  // Writing the maximum gray value 
    int count = 0; 
    
    for (int i = 0; i < HEIGHT; i++) { 
        for (int j = 0; j < WIDTH; j++) { 
            temp = image[i][j]; 
            fprintf(pgmimg, "%d ", temp); // Writing the gray values in the 2D array to the file 
        } 
        fprintf(pgmimg, "\n"); 
    } 
    fclose(pgmimg); 
} 

int cal_pixel(complex c) {

            double z_real = 0;
            double z_imag = 0;
            
            double z_real2, z_imag2, lengthsq;

            int iter = 0;
            do {
                z_real2 = z_real * z_real;
                z_imag2 = z_imag * z_imag;

                z_imag = 2 * z_real * z_imag + c.imag;
                z_real = z_real2 - z_imag2 + c.real;
                lengthsq =  z_real2 + z_imag2;
                iter++;
            }
            while ((iter < MAX_ITER) && (lengthsq < 4.0));

            return iter;

}




int main(int argc, char *argv[]){
double AVG = 0;

int IMAGE[HEIGHT][WIDTH];
double real_max = 2.0, real_min = -2.0, imag_max = 2.0, imag_min = -2.0;
int myRank, worldSize;
MPI_Init(NULL, NULL);

MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

int processRows = HEIGHT/worldSize;

int myStart = myRank * processRows;
int myEnd = myStart + processRows;

int localImage[processRows][WIDTH];
complex c;
double scale_real = (real_max - real_min) / WIDTH;
double scale_imag = (imag_max - imag_min) / HEIGHT;

for (int k = 0;k < 10; k++){
clock_t start_time = clock();

for (int i = myStart; i < myEnd; i++){
for (int j = 0; j < WIDTH; j++){
c.real = real_min + ((double) j * scale_real);
c.imag = imag_min + ((double) i * scale_imag);
localImage[i-myStart][j] = cal_pixel(c);

}
}
MPI_Gather(localImage, processRows * WIDTH, MPI_INT, IMAGE, processRows * WIDTH, MPI_INT, 0, MPI_COMM_WORLD);
if (myRank==0){
clock_t end_time = clock();
double total_time;
total_time = ((double) (end_time - start_time))/CLOCKS_PER_SEC;
printf("Elapsed time of trial %d: %f\n",k+1, total_time);
AVG += total_time;

}



}

if (myRank==0){
save_pgm("staticMandel.pgm", IMAGE);
printf("The average execution time of 10 trials is: %f ms", AVG/10*1000);
}


MPI_Finalize();
return 0;

}








