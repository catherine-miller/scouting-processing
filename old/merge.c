//reads files time-multiplexed by 6 and writes into one file
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
    const char *fromfile1 = "files/out.dump.bx_0_of_6.dump";
    const char *fromfile2 = "files/out.dump.bx_1_of_6.dump";
    const char *fromfile3 = "files/out.dump.bx_2_of_6.dump";
    const char *fromfile4 = "files/out.dump.bx_3_of_6.dump";
    const char *fromfile5 = "files/out.dump.bx_4_of_6.dump";
    const char *fromfile6 = "files/out.dump.bx_5_of_6.dump";
    const char *tofile = "files/puppi_SingleNeutrino_PU200_out_merge.dump";

    //open the files
    int fromfd[6];
    fromfd[0] = open(fromfile1, O_RDONLY);
    fromfd[1] = open(fromfile2, O_RDONLY);
    fromfd[2] = open(fromfile3, O_RDONLY);
    fromfd[3] = open(fromfile4, O_RDONLY);
    fromfd[4] = open(fromfile5, O_RDONLY);
    fromfd[5] = open(fromfile6, O_RDONLY);





    //fstat(fromfd, &stat_buf);
    int tofd = open(tofile, O_WRONLY | O_CREAT, 00644);
    if (tofd < 0) {
        printf("ERROR opening file \n");
        return 1;
    }
    //read the first header

    int n = 8; //the result of read
    size_t j = 0; //the position in the buffer to read into
    unsigned long long particles; //the number of particles from the header
    int orb = 0;
    int part;
    int n2 = 0;
    int k;
    int w;
    //read all the particles from one event plus the next event's header
    //the first time just read the first header
    while (n > 0 && orb < 100) {
        //k % 6 keeps track of file to read from (luckily 3564 is divisible by 6)
        for (k = 0; k < 3564; ++k) {
            n = read(fromfd[k % 6], &buf[j], 8);
            j += n;
            if (n != 8) {
                printf("ERROR reading header \n");
                printf("Event: %d, Orbit: %d \n", k, orb);
                break;
            }
            memcpy(&particles, &buf[j], 8);
            particles &= ((1 << 12) - 1); //should clear any bits above 12 (ie parts of the header not the particle)
            n2 = read(fromfd[k % 6], &buf[j], particles * 8);
            j += n2;
            if (n2 != particles * 8) {
                printf("ERROR reading particles \n");
                printf("Event: %d, Orbit: %d, Particles: %d", k, orb, particles);
                break;
            } 
        }
        w = write(tofd, &buf, j); //at the end of the orbit, write to file
        if (w != j) {
            printf("ERROR writing: attempted to write %d bytes, wrote %d \n", j, w);
            printf("Event: %d, Orbit: %d, Particles: %d\n", k, orb, particles);
            return 1;
        }
        j = 0;
        ++orb;
    }
    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n", time_used);
    printf("Orbit: %d, n: %d, last call to write: %d \n", orb, n, w);
    return 0;
}
