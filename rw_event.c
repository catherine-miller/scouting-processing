//reads events into buffer and writes into file one by one

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#define BUF_SIZE 8*1000
int main() {
    clock_t start, end;
    double time_used;
    start = clock ();
    char buf[BUF_SIZE];
    //const char *fromfile = "files/puppi_SingleNeutrino_PU200.dump";
    const char *fromfile = "files/puppi_SingleNeutrino_PU200_out_merge.dump";
    //const char *fromfile = "files/puppi_SingleNeutrino_PU200_out_merge.dump";
    const char *tofile = "files/puppi_SingleNeutrino_PU200_out_byevent.dump";
    int fromfd = open(fromfile, O_RDONLY);
    int tofd = open(tofile, O_WRONLY | O_CREAT, 00644);


    size_t bytes_to_read = 8;
    int k = 0; //the bx
    unsigned long long particles; //the number of particles from the header

    //read all the particles from one event plus the next event's header
    //the first time just read the first header
    int n;
    while (n = read(fromfd, &buf, bytes_to_read) > 0) {
        memcpy(&particles, &buf[bytes_to_read-8], 8);
        particles &= ((1 << 12) - 1); //should clear any bits above 12 (ie parts of the header not the particle)
        k += 1;
        write(tofd, &buf, bytes_to_read);
        bytes_to_read = (particles + 1) * 8; //next event read number of particles
    }

    printf("events written: %i \n",k);
    k = 0;

    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n",time_used);
    printf("n = %d, bytes_to_read = %d, particles = %d \n", n, bytes_to_read, particles);
    return 0;
}