//most recent
//reads files time-multiplexed by 6
//distributes orbits into 6 files

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#define BUF_SIZE 8*50*3564 //this is the maximum size of an orbit
int main() {
    clock_t start, end;
    double time_used;
    start = clock ();
    char buf[BUF_SIZE];

    //files to read from
    const char *fromfile1 = "files/out.dump.bx_0_of_6.dump";
    const char *fromfile2 = "files/out.dump.bx_1_of_6.dump";
    const char *fromfile3 = "files/out.dump.bx_2_of_6.dump";
    const char *fromfile4 = "files/out.dump.bx_3_of_6.dump";
    const char *fromfile5 = "files/out.dump.bx_4_of_6.dump";
    const char *fromfile6 = "files/out.dump.bx_5_of_6.dump";

    //open the files
    int fromfd[6];
    fromfd[0] = open(fromfile1, O_RDONLY);
    fromfd[1] = open(fromfile2, O_RDONLY);
    fromfd[2] = open(fromfile3, O_RDONLY);
    fromfd[3] = open(fromfile4, O_RDONLY);
    fromfd[4] = open(fromfile5, O_RDONLY);
    fromfd[5] = open(fromfile6, O_RDONLY);

    //output files
    const char *tofile1 = "files/mergesplit1.dump";
    const char *tofile2 = "files/mergesplit2.dump";
    const char *tofile3 = "files/mergesplit3.dump";
    const char *tofile4 = "files/mergesplit4.dump";
    const char *tofile5 = "files/mergesplit5.dump";
    const char *tofile6 = "files/mergesplit6.dump";

    int tofd[6];
    tofd[0] = open(tofile1, O_WRONLY | O_CREAT, 00700);
    tofd[1] = open(tofile2, O_WRONLY | O_CREAT, 00700);
    tofd[2] = open(tofile3, O_WRONLY | O_CREAT, 00700);
    tofd[3] = open(tofile4, O_WRONLY | O_CREAT, 00700);
    tofd[4] = open(tofile5, O_WRONLY | O_CREAT, 00700);
    tofd[5] = open(tofile6, O_WRONLY | O_CREAT, 00700);

    int n = 8; //the result of read
    size_t j = 0; //the position in the buffer to read into
    unsigned long long header; //the number of particles from the header
    int orb = 0;
    int part;
    int n2 = 0;
    int w; //result of write
    int k;
    //read all the particles from one event plus the next event's header
    //the first time just read the first header
    while (n > 0 && orb < 100) {
        //k % 6 keeps track of file to read from (luckily 3564 is divisible by 6)
        for (k = 0; k < 3564; ++k) { //k is the bunch crossing
            n = read(fromfd[k % 6], &buf[j], 8); //read header into buffer
            j += n;
           /* if (n != 8) {
                printf("ERROR reading header \n");
                printf("Event: %d, Orbit: %d \n", k, orb);
                break;
            } */
            memcpy(&header, &buf[j], 8);
            part = header & ((1 << 12) - 1); //should clear any bits above 12 (ie parts of the header not the particle)
            n2 = read(fromfd[k % 6], &buf[j], part * 8); //read appropriate number of particles
            j += n2;
          /*  if (n2 != part * 8) {
                printf("ERROR reading particles \n");
                printf("Event: %d, Orbit: %d, Particles: %d", k, orb, part);
                break;
            } */
            /*if (part > 256) {
                printf("ERROR: Asked to read out %llu particles, but the maximum number of particles is 256 \n", part);
                printf("Event = %d, Orbit = %d \n", k, orb + 1);
                //return 1;
            }*/
        }
    w = write(tofd[orb++ % 6], &buf, j); //at the end of the orbit, write to file
   /* if (w != j) {
        printf("ERROR writing: attempted to write %d bytes, wrote %d \n", j, w);
        printf("Event: %d, Orbit: %d, Particles: %d\n", k, orb, part);
        //printf("SSIZE_MAX: %d", SSIZE_MAX)
        return 1;
    } */
    j = 0;
    }
    //printf("events in the last orbit: %i \n",k);
    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n", time_used);
    printf("Orbit: %d, n: %d, last call to write: %d \n", orb, n, w);
    return 0;
}
