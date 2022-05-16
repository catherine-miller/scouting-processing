//reads part of the dummy files I am using for these exercises

#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    
    if (argc < 1) {
        printf("ERROR: no file name to read \n");
        return 1;
    }
    const char *filetoread = argv[1];
    FILE *fptr;

    fptr = fopen(filetoread, "r");
    if (fptr == NULL) {
        printf("Error opening the file to read :( \n");
        return 1;
    }
    unsigned long long nparticlesout;
    //fread(&nparticlesout, 8, 1, fptr); //read out the number of particles

    char pout[256*8];
    int bx = 0;
    while (fread(&nparticlesout, 8, 1, fptr) > 0 && bx < 10 ) {
        nparticlesout &= 0xFFF;
        printf("Event %i: %llu particles \n", bx, nparticlesout);
        fread(&pout, 8, (size_t) nparticlesout, fptr); //read out puppi candidates
        bx += 1;
    }
    int orbits = (int) (bx / 3564);
    printf("Orbits: %i \n", orbits);
    return 0;
}
