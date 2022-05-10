//reads from dump file of simulated events in the L1 trigger
//into a buffer w the average size of an orbit
//writes out "orbit" by "orbit"

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <fcntl.h>
#define BUF_SIZE 8*31*3564 //this is the average size of an orbit
int main() {
    clock_t start, end;
    double time_used;
    start = clock ();
    char buf[BUF_SIZE];
    const char *fromfile = "files/puppi_SingleNeutrino_PU200.dump";
    const char *tofile = "files/puppi_SingleNeutrino_PU200_out_byorbit.dump";
    //struct stat stat_buf;
    int fromfd = open(fromfile, O_RDONLY);
    //fstat(fromfd, &stat_buf);
    int tofd = open(tofile, O_WRONLY | O_CREAT);
    int n;
    while ((n = read(fromfd, &buf, sizeof(buf))) > 0) {
	    write(tofd, &buf, n);
    }
    end = clock();
    time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Reading and writing took %f s \n",time_used);
    return 0;
}
