//reads events into buffer one by one and writes out whole orbits

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#define BUF_SIZE 8*50*3564
int main() {
    clock_t start, end;
    double time_used;
    start = clock ();
    char buf[BUF_SIZE];
    //const char *fromfile = "files/puppi_SingleNeutrino_PU200.dump";
    const char *fromfile = "files/puppi_SingleNeutrino_PU200_out_merge.dump";
    const char *tofile = "files/puppi_SingleNeutrino_PU200_out_byorbit.dump";
    int fromfd = open(fromfile, O_RDONLY);
    int tofd = open(tofile, O_WRONLY | O_CREAT, 00644);

    int n = 8; //the result of read
    int j; //the position in the buffer to read into
    size_t bytes_to_read = 8;
    int k; //the bx
    unsigned long long particles; //the number of particles from the header
    int orb = 0;
    int w;
    int n2;


    while (n > 0 && orb < 100) {
        //k % 6 keeps track of file to read from (luckily 3564 is divisible by 6)
        for (k = 0; k < 3564; ++k) {
            n = read(fromfd, &buf[j], 8);
            j += n;
          /*  if (n != 8) {
                printf("ERROR reading header \n");
                printf("Event: %d, Orbit: %d, Bytes read: %d \n", k, orb, n);
                break;
            } */
            memcpy(&particles, &buf[j], 8);
            particles &= ((1 << 12) - 1); //should clear any bits above 12 (ie parts of the header not the particle)
            n2 = read(fromfd, &buf[j], particles * 8);
            j += n2;
          /* if (n2 != particles * 8) {
                printf("ERROR reading particles \n");
                printf("Result of read: %d Event: %d, Orbit: %d, Particles * 8: %llu \n", n2, k, orb, particles * 8);
                break;
            } */
        }
        w = write(tofd, &buf, j); //at the end of the orbit, write to file
      /*  if (w != j) {
            printf("ERROR writing: attempted to write %d bytes, wrote %d \n", j, w);
            printf("Event: %d, Orbit: %d, Particles: %d\n", k, orb, particles);
            return 1;
        } */
        j = 0;
        ++orb;
    }

    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n", time_used);
    return 0;
}
