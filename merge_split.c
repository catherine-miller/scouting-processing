//OLD
//reads files time-multiplexed by 6 (by event)
//distributes into n files by orbit
//reads events into buffer one by one and writes out whole orbits

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#include <string.h>
#define BUF_SIZE 8*256*3564
int main(int argc, char *argv[]) {
    clock_t start, end;
    double time_used;
    start = clock ();
    char buf[BUF_SIZE];

    //files to read from and write to
    const char *fromfile1 = "files/out.dump.bx_0_of_6_.dump";
    const char *fromfile2 = "files/out.dump.bx_1_of_6_.dump";
    const char *fromfile3 = "files/out.dump.bx_2_of_6_.dump";
    const char *fromfile4 = "files/out.dump.bx_3_of_6_.dump";
    const char *fromfile5 = "files/out.dump.bx_4_of_6_.dump";
    const char *fromfile6 = "files/out.dump.bx_5_of_6_.dump";

    int fromfd[6];
    fromfd[0] = open(fromfile1, O_RDONLY);
    fromfd[1] = open(fromfile2, O_RDONLY);
    fromfd[2] = open(fromfile3, O_RDONLY);
    fromfd[3] = open(fromfile4, O_RDONLY);
    fromfd[4] = open(fromfile5, O_RDONLY);
    fromfd[5] = open(fromfile6, O_RDONLY);

    const char *tofile1 = "files/mergesplit1.dump";
    const char *tofile2 = "files/mergesplit2.dump";
    const char *tofile3 = "files/mergesplit3.dump";
    const char *tofile4 = "files/mergesplit4.dump";
    const char *tofile5 = "files/mergesplit5.dump";
    const char *tofile6 = "files/mergesplit6.dump";

    int tofd[6];
    tofd[0] = open(tofile1, O_WRONLY | O_CREAT);
    tofd[1] = open(tofile2, O_WRONLY | O_CREAT);
    tofd[2] = open(tofile3, O_WRONLY | O_CREAT);
    tofd[3] = open(tofile4, O_WRONLY | O_CREAT);
    tofd[4] = open(tofile5, O_WRONLY | O_CREAT);
    tofd[5] = open(tofile6, O_WRONLY | O_CREAT);

/*
    int nout = argv[1] - '0';
    int tofd[nout];
    for (int i = 0; i < nout; ++i) {
        char tofile[50];
        //char prefix = "files/merge_split";
        //char extension = ".dump";
        sprintf(tofile,"%s","files/merge_split");
        sprintf(tofile, "%d", i+1);
        sprintf(tofile, "%s", ".dump");
        const char *tofile2 = tofile;
        tofd[i] = open(tofile2, O_WRONLY | O_CREAT);
    }
*/



    //read the first header

    int n = 8; //the result of read
    int j = 0; //the position in the buffer to read into
    int k = 0; //the bx
    int particles; //the number of particles from the header
    int orb = 0;

    //read all the particles from one event plus the next event's header
    //the first time just read the first header
    while (n > 0) {
        //k % 6 keeps track of file to read from (luckily 3564 is divisible by 6)
        while ((n = read(fromfd[k % 6], &buf[j], 8) == 8) && (k < 3564)) { //read header
            memcpy(&particles, &buf[j+7], 1);
            j += n;
            n = read(fromfd[k % 6], &buf[j], particles * 8);
            if (particles > 256) {
                printf("ERROR: Asked to read out %llu particles, but the maximum number of particles is 256 \n", particles);
                return 1;
            }
            k += 1;
            j += n;
        }
    if (n != particles * 8 && n != -1) j += n;
    write(tofd[orb++ % 6], &buf, j); //at the end of the orbit, write to file, then switch the file you're writing to for next time
    //printf("events in the orbit: %i \n",k);
    k = 0;
    j = 0;
    }
    //printf("events in the last orbit: %i \n", k);
    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n", time_used);
    return 0;
}