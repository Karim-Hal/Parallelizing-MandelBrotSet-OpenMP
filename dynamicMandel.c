#include <stdio.h>
#include <time.h>
#include "mpi.h"


#define WIDTH 640
#define HEIGHT 480
#define MAX_ITER 255
#define DATA_TAG 0
#define TERMINATION_TAG 2

typedef struct c{
  double real;
  double imag;
}complex ;


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


int main(int argc, char *argv[]){

int IMAGE[HEIGHT][WIDTH];
double AVG = 0;
double real_max = 2.0, real_min = -2.0, imag_max = 2.0, imag_min = -2.0;
int myRank, worldSize;
MPI_Init(NULL, NULL);

MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
double scale_real = (real_max - real_min) / WIDTH;
double scale_imag = (imag_max - imag_min) / HEIGHT;

for (int k = 0; k < 10; k++){

clock_t start_time = clock();
//Master
if (myRank==0){
int row =0;
int count = 0;

int result[WIDTH];
MPI_Status status;


for (int i = 1; i < worldSize; i++){

MPI_Send(&row, 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
count++;
row++;
}

do{
MPI_Recv(result, WIDTH, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
count--;
int source = status.MPI_SOURCE;
if (row < HEIGHT){

MPI_Send(&row, 1, MPI_INT,source , DATA_TAG, MPI_COMM_WORLD);
row++;
count++;

}

int rowDone = status.MPI_TAG;
for (int r =0; r < WIDTH; r++){
IMAGE[rowDone][r] = result[r];
}
}while(count>0);

clock_t end_time = clock();
double total_time;

for (int i = 1; i < worldSize; i++) {
        MPI_Send(&row, 1, MPI_INT, i, TERMINATION_TAG, MPI_COMM_WORLD);
    }
    
total_time = ((double) (end_time - start_time))/CLOCKS_PER_SEC;
printf("Elapsed time of trial %d: %f\n",k+1, total_time);
AVG += total_time;


}


//Workers
else{
int row;
int result[WIDTH];
MPI_Status status;
complex c;

MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
int tag = status.MPI_TAG;
while(tag==DATA_TAG){

c.imag = imag_min + ((double) row * scale_imag);
for (int x = 0; x < WIDTH; x++){
c.real = real_min + ((double) x * scale_real);

result[x] = cal_pixel(c);
}

MPI_Send(result, WIDTH, MPI_INT, 0, row, MPI_COMM_WORLD);
MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
tag = status.MPI_TAG;

}

}

}

if (myRank==0){
save_pgm("dynamicMandel.pgm", IMAGE);
printf("The average execution time of 10 trials is: %f ms", AVG/10*1000);
}


MPI_Finalize();
return 0;




}
