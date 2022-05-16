#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <cassert>
#include <thread>
#include <atomic>
#include <vector>
#include <oneapi/tbb.h>
#include <oneapi/tbb/parallel_pipeline.h>

const unsigned int BX_PER_ORBIT = 3564;

typedef std::chrono::time_point<std::chrono::steady_clock> tp; 
void report(unsigned long long totorb, unsigned long long totev, unsigned long long totpuppi, 
            tp t0, tp t1, tp t2) {
    unsigned long long totbytes = (totpuppi + totev)*8;
    std::chrono::duration<double> dt1 = t2 - t1, dt0 = t1-t0;
    printf("Read %llu orbits, %llu events, %llu candidates, %llu bytes\n", totorb, totev, totpuppi, totbytes);
    printf("Time to init %.2f ms, to read file %.2f ms\n", dt0*1000, dt1*1000);
    printf("Event rate %.2f kHz, data rate %.3f GB/s\n",
                totev/1000./dt1.count(), totbytes/1024.0/1024.0/1024.0/dt1.count());
}

int sequential_1toNscatter(const char *fname, int nouts, const char *foutname) {
    assert(nouts >= 1);
    unsigned long long int totorb = 0, totev = 0, totpuppi = 0;
    auto t0 = std::chrono::steady_clock::now();
    std::fstream fin(fname, std::ios_base::in | std::ios_base::binary); //open the file (?)
    bool write = true;
    std::vector<std::fstream> fouts;
    if (std::string(foutname) == "null") {
        write = false;
        printf("Will pretend to distribute data from %s to %d sinks\n", fname, nouts);
    } else {
        printf("Will distribute data from %s to %d sinks:\n", fname, nouts);
        char fnamebuff[1024];
        fouts.resize(nouts); //resize vector of output files to # of output files
        for (int i = 0; i < nouts; ++i) {
            snprintf(fnamebuff,1023,foutname,i); // print file names to a buffer
            fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
            printf(" - %s\n", fnamebuff); //print file names to terminal
        }
        printf("\n");
    }
    uint64_t header, data[255], npuppi;
    uint32_t orbit, lastorbit = 1<<30; //ig the last orbit is 2^30?
    int iout = -1;
    auto t1 = std::chrono::steady_clock::now();    
    while(fin.good()) {
        fin.read(reinterpret_cast<char*>(&header), sizeof(uint64_t)); //I think read the header into header
        if (fin.gcount() == 0) break;
        orbit  = (header >> 24) & 0x3FFFFFFF; //takes the bits of the header between 3 and 30? or something
        if (orbit != lastorbit) { //why would you repeat orbits anyway?
            totorb++; 
iout = (iout + 1) % nouts;
            lastorbit = orbit;
        }
        if (write) fouts[iout].write(reinterpret_cast<char*>(&header), sizeof(uint64_t));
        npuppi = header & 0xFFF; //the last 12? bytes
        totev++;
        if (npuppi) {
            fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
            assert(fin.gcount() == npuppi*sizeof(uint64_t)); //program will fail if last input byte # not same as Puppi*8
            if (write) fouts[iout].write(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
            totpuppi += npuppi;
        }
    }
    if (write) {
        for (auto & f : fouts) f.close();
    }
    auto t2 = std::chrono::steady_clock::now();
    report(totorb, totev, totpuppi, t0, t1, t2);
    return 0;
}

class OrbitBuffer {
    public:
        OrbitBuffer(int bufsizeMb) { //makes a buffer using assigned_alloc
        //also contains the methods that read
            assert(bufsizeMb >= 1);
            size_t size = bufsizeMb*1024*1024;
            mem_ = reinterpret_cast<char *>(std::aligned_alloc(4096u, size)); //this is assigning the memory for the buffer
            //I guess we want to do things so explicitly in this case b/c we're dealing with a lot of memory
            //and later will deallocate?
            assert(mem_ != nullptr);
            ptr_ = mem_;
            writeptr = mem_;
            end_ = mem_ + size;
        }
        ~OrbitBuffer() { //deallocates the memory
            std::free(mem_);
        }
        void clear() {
            ptr_ = mem_;
        }
        bool writeptrclear() {
            if (writeptr != ptr_) { //make sure we wrote out everything
                printf("ERROR: writeptr != ptr_ \n");
                printf("writeptr - mem_ = %d, ptr - mem_ = %d \n", writeptr-mem_, ptr_-mem_);
                //printf("orbit = %d, bx = %d \n", ctx->totorb, ctx->totevent);
                return false;
            }
            //initialize both pointers to mem_, the start of buffer
            ptr_ = mem_;
            writeptr = mem_;
            //printf("initialized pointers \n"); fflush(stdout);
            return true;
        }
        bool empty() const { return ptr_ == mem_; } //ptr_ keeps track of where we are w/in allocated memory
        bool full() const { return ptr_ == end_; }
        void write64(uint64_t word) {
            uint64_t *ptr64 = reinterpret_cast<uint64_t *>(ptr_);
            *ptr64 = word; //write 64-bit word to the address indicated by ptr_
            ptr_ += 8; 
            assert(ptr_ < end_);
        }
        void readChunk(std::fstream & in, size_t length) {
            in.read(ptr_, length); //read from in file
            assert(in.gcount() == length); //check that we did read length bytes
            ptr_ += length; //move the address to read into next time
            assert(ptr_ < end_); //check we're not at the end of allocated mem
        }
        size_t dataSize() const { return ptr_ - mem_; } //how much we have in
        void writeToFile(std::fstream & out) {
            out.write(mem_, dataSize());
        }
        void writeevent(std::fstream & out) { //writes one event from the buffer into a file
            uint64_t header;
            memcpy(&header, writeptr, 8);
            uint32_t npuppi = header & 0xFFF;
            out.write(writeptr, 8*(npuppi+1));
            writeptr += 8*(npuppi+1);
        }
       /* void write_demux(std::vector<std::fstream> &outs, unsigned int nin, unsigned int nout, Merge_Context::Pair p) {
            uint64_t npuppi;
            for (unsigned int i = 0; i < 3564; ++i) {
                unsigned int ibuff = (i + p.first) % nin; //orbit doesn't always start at 0
                memcpy(&npuppi, &)
            }
        } */ //simply not possible within one orbitbuffer lol
        bool readNEvents(std::fstream & fin, 
                         unsigned int n, 
                         unsigned long long int & totev, //why would you want to use reference to your args?
                         unsigned long long int & totpuppi,
                         unsigned long long int myorbit) {
            uint64_t header, npuppi;
            uint32_t orbit;
            for (unsigned int i = 0; fin.good() && i < n; ++i) {
                fin.read(reinterpret_cast<char*>(&header), 8); //read header from file into "header" 
                if (fin.gcount() == 0) break; //break if we read nothing
                //getting orbit (&checking), number of puppi, etc from file
                orbit = (header >> 24) & 0x3FFFFFFF;
                npuppi = header & 0xFFF;
                uint32_t bx = (header >> 12) & 0xFFF;
                //printf("orbit: %lu \n", orbit);
                //printf("bx: %lu \n", bx);
                //printf("npuppi: %llu \n", npuppi);
                /*if (i == 0) myorbit = orbit;
                else */ if (orbit != myorbit) {
                    printf("ERROR: orbit from header: %lu, expected orbit: %d \n", orbit, myorbit);
                    printf("bunch crossings: %d, i = %d \n", totev, i);
                    exit(1);
                }
                write64(header);
                totev++;
                if (npuppi) {
                    readChunk(fin, npuppi*sizeof(uint64_t)); //read
                    totpuppi += npuppi;
                }
            }
            return fin.good();
        }
        //what orbit header??
        //looks like there is supposed to be some header with the size of the entire orbit
        bool readOrbitWithHeader(std::fstream & fin, 
                         unsigned long long int & totpuppi) {
            uint64_t orbitHeader, orbit;
            fin.read(reinterpret_cast<char*>(&orbitHeader), sizeof(uint64_t)); //read orbit header into orbitHeader
            if (fin.gcount() == 0) return false;
            orbit  = (orbitHeader >> 24) & 0x3FFFFFFF;
            size_t ndata = orbitHeader & ((1 << 24)-1), length = ndata * sizeof(uint64_t);
            //printf("orbit %llu of length %llu\n", orbit, ndata);
            write64(orbitHeader); //write orbit header into buffer
            fin.read(ptr_, length);
            assert(fin.gcount() == length);
            ptr_ += length;
            totpuppi += ndata;
            return fin.good();
        }

    private:
        char *mem_, *ptr_, *end_, *writeptr;
};

int sequential_fullorbit_1toNscatter(const char *fname, int nouts, int bufsizeMb, const char *foutname) {
    assert(nouts >= 1);
    unsigned long long int totorb = 0, totev = 0, totpuppi = 0;
    auto t0 = std::chrono::steady_clock::now();
    std::fstream fin(fname, std::ios_base::in | std::ios_base::binary);
    bool write = true;
    printf("Will distribute data from %s to %d sinks with one %dMB orbit buffer:\n", fname, nouts, bufsizeMb);
    std::vector<std::fstream> fouts;
    if (std::string(foutname) == "null") {
        write = false;
    } else {
        char fnamebuff[1024];
        fouts.resize(nouts);
        for (int i = 0; i < nouts; ++i) {
            snprintf(fnamebuff,1023,foutname,i);
            fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
            printf(" - %s\n", fnamebuff);
        }
        printf("\n");
    }
    OrbitBuffer buffer(bufsizeMb);
    unsigned int nEvts = BX_PER_ORBIT;
    auto t1 = std::chrono::steady_clock::now();    
    bool fileok;
    int iout = 0;
    do {
        buffer.clear();
        fileok = buffer.readNEvents(fin, nEvts, totev, totpuppi, totorb);
        totorb += 1;
        if (write) buffer.writeToFile(fouts[iout]);
        iout = (iout+1)%nouts;
    } while (fileok);
    if (write) {
        for (auto & f : fouts) f.close();
    }
    auto t2 = std::chrono::steady_clock::now();
    report(totorb, totev, totpuppi, t0, t1, t2);
    return 0;
}
int sequential_fullorbitOH_1toNscatter(const char *fname, int nouts, int bufsizeMb, const char *foutname) {
    assert(nouts >= 1);
    unsigned long long int totorb = 0, totev = 0, totpuppi = 0;
    auto t0 = std::chrono::steady_clock::now();
    std::fstream fin(fname, std::ios_base::in | std::ios_base::binary);
    bool write = true;
    printf("Will distribute data from %s to %d sinks with one %dMB orbit buffer:\n", fname, nouts, bufsizeMb);

    //making and naming out files
    std::vector<std::fstream> fouts;
    if (std::string(foutname) == "null") {
        write = false;
    } else {
        char fnamebuff[1024];
        fouts.resize(nouts);
        for (int i = 0; i < nouts; ++i) {
            snprintf(fnamebuff,1023,foutname,i);
            fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
            printf(" - %s\n", fnamebuff);
        }
        printf("\n");
    }
    OrbitBuffer buffer(bufsizeMb);
    auto t1 = std::chrono::steady_clock::now();    
    bool fileok;
    int iout = 0;
    do {
        buffer.clear();
        fileok = buffer.readOrbitWithHeader(fin,totpuppi); //reads whole orbit at once
        totorb += 1; totev += BX_PER_ORBIT;
        if (write) buffer.writeToFile(fouts[iout]);
        iout = (iout+1)%nouts;
    } while (fileok);
    if (write) {
        for (auto & f : fouts) f.close();
    }
    auto t2 = std::chrono::steady_clock::now();
    report(totorb, totev, totpuppi, t0, t1, t2);
    return 0;
}

class ManySinks_Context { 
    //opens files and names them, makes nouts + 2 buffers
    //function orbitbuffer gives address of next buffer
    //nextpair returns address of next buffer and next output file number
    public:
        typedef std::pair<OrbitBuffer *, int> Pair;
        ManySinks_Context(const char *fname, int nout, int bufsizeMb, const char *foutname) :
            fin(fname, std::ios_base::in | std::ios_base::binary),
            nouts(nout),
            fouts(nouts),
            write(true),
            totorb(0), totev(0), totpuppi(0),
            nbuffs(nouts+2),
            ibuff(0), iout(0)
        {
            printf("Will distribute data from %s to %d sinks with one %dMB orbit buffer:\n", fname, nouts, bufsizeMb);
            if (std::string(foutname) == "null") {
                write = false;
            } else {
                char fnamebuff[1024];
                for (int i = 0; i < nouts; ++i) {
                    snprintf(fnamebuff,1023,foutname,i);
                    fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                }
                printf("\n");
            }
            buffers.resize(nbuffs);
            for (int i = 0; i < nbuffs; ++i) {
                buffers[i] = new OrbitBuffer(bufsizeMb);
            }
            ibuff = 0;
            iout = 0;
        }

        OrbitBuffer *nextBuffer() { //pointer to the next buffer
            int curr = ibuff;
            while (!(ibuff.compare_exchange_weak(curr,(curr+1)%nbuffs))) ; //add 1 or go back to 0
            OrbitBuffer *ret = buffers[curr];
            //assert(ret->empty());
            return ret;
        }
        int nextOut() { //don't understand what's supposed to change between loops here?
            int curr, next;
            while(true) {
                curr = iout;
                next = (iout+1)%nouts;
                if (iout.compare_exchange_weak(curr,next)) {
                    return curr;
                }
            } 
        }

        Pair nextPair() { return Pair(nextBuffer(), nextOut()); }
    
    public:
        std::fstream fin;
        unsigned int nouts;
        std::vector<std::fstream> fouts;
        bool write;
        unsigned long long int totorb, totev, totpuppi;
    private:
        unsigned int nbuffs;
        std::vector<OrbitBuffer *> buffers;
        std::atomic<int> ibuff, iout; //atomic: access guaranteed not to cause "data races," can use to synchronize data access across threads
        uint32_t lastorbit;

};

class MergeSplit_Context {
    public:
       //typedef std::tuple<uint, uint, uint> Triplet; //contains a buffer for each input file and an index corresponding to the output file
        MergeSplit_Context(int opt, const char* finname, int nin, int nout, int bufsizeMb, const char *foutname) :
            //fin[nin](fname, std::ios_base::in | std::ios_base::binary),
            nouts(nout),
            nins(nin),
            //fouts(nouts),
            write(true),
            totorb(0), totev(0), totpuppi(0),
            nbuffs(nin * (nout + nin + 1)),
            iout(0),
            start(0),
            bufdex(0)
        {
            printf("Will distribute data from %d to %d files with %d %dMB orbit buffers.\n", nin, nouts, nbuffs, bufsizeMb);
            if (std::string(foutname) == "null" or std::string(finname) == "null") {
                write = false;
            } else {
                char fnamebuff[1024];
                fins.resize(nins);
                if (opt == 0) {
                //make files
                //construct names of input files (must follow the template *n.dump where n is an index 0 to N-1)
                printf("Input files:\n");
                for (int i = 0; i < nins; ++i) {
                    snprintf(fnamebuff,1023, "%s%d.dump",finname,i);
                    fins[i].open(fnamebuff, std::ios_base::in | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                } } else if (opt == 1) {
                    if (nins != 1) {
                        printf("ERROR: Method 1 takes only one file \n");
                        exit(1);
                    }
                    fins[0].open(finname, std::ios_base::in | std::ios_base::binary);
                } else if (opt == 2) {
                    printf("Input files:\n");
                for (int i = 0; i < nins; ++i) {
                    snprintf(fnamebuff,1023, "%s.bx_%d_of_%d.dump",finname,i,nout);
                    fins[i].open(fnamebuff, std::ios_base::in | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                } } else {
                    printf("Invalid file name format \n");
                    exit(1);
                }
                //make output files
                fouts.resize(nout);
                printf("Output files:\n");
                for (int i = 0; i < nouts; ++i) {
                    snprintf(fnamebuff,1023, "%s%d.dump",foutname,i);
                    fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                }
                }
                printf("\n");
            buffers.resize(nbuffs);
            for (int i = 0; i < nbuffs; ++i) {
                buffers[i] = new OrbitBuffer(bufsizeMb);
            }
            iout = 0;
        }

       /* OrbitBuffer *nextBuffer() { //pointer to the next buffer
            int curr = ibuff;
            while (!(ibuff.compare_exchange_weak(curr,(curr+1)%nbuffs))) ;
            //add 1 to ibuff or go back to 0
            OrbitBuffer *ret = buffers[curr];
            //assert(ret->empty());
            return ret;
        }*/

        uint nextOut() { //the next file to write out to
            uint curr, next;
            while(true) {
                curr = iout;
                next = (iout+1)%nouts;
                if (iout.compare_exchange_weak(curr,next)) {
                    return curr;
                }
            } 
        }

        uint OrbitStart() {
            uint c, n;
            while(true) {
                c = start;
                n = start + ((3564 % nins) * totorb) % nins;
                if (start.compare_exchange_weak(c, n)) return c;
            }
        }

        uint nextBuf() { 
            uint curbuf, nextbuf;
            while (true) {
                curbuf = bufdex;
                nextbuf = (curbuf + 1) % (nouts + nins + 1);
                if (bufdex.compare_exchange_weak(curbuf, nextbuf)) return curbuf;
            }
        }
        uint nextOrbit() {
            uint curorbit, nextorbit;
            while (true) {
                curorbit = totorb;
                nextorbit = curorbit + 1;
                if (totorb.compare_exchange_weak(curorbit, nextorbit)) return curorbit;
            }
        }

        std::vector<uint> nextQuad() {
            std::vector<uint> v{nextOut(), OrbitStart(), nextBuf(), nextOrbit()}; 
            return v;
            }

    public:
        std::vector<std::fstream> fins;
        unsigned int nouts;
        unsigned int nins;
        std::vector<std::fstream> fouts;
        bool write;
        unsigned long long int totev, totpuppi;
        std::atomic<uint> totorb;
        std::vector<OrbitBuffer *> buffers; //is it ok if this is public?
    private:
        unsigned int nbuffs;
        std::atomic<uint> iout, start, bufdex; //atomic: access guaranteed not to cause "data races," can use to synchronize data access across threads
        //std::vector<std::atomic<int>> eventstoread[fins];
        uint32_t lastorbit;

};

struct ManySinks_Reader {
    ManySinks_Reader(ManySinks_Context &context) : ctx(&context) {} //takes a pointer to context, which is now called ctx (?)
    ManySinks_Context::Pair operator()(oneapi::tbb::flow_control& fc) const {
        ManySinks_Context::Pair ret(nullptr, 0);
        if (!ctx->fin.good()) {
            fc.stop();
            return ret;
        }
        ret = ctx->nextPair();
        assert(ret.first->empty());
        bool fileok = ret.first->readNEvents(ctx->fin, BX_PER_ORBIT, ctx->totev, ctx->totpuppi, ctx->totorb);
        ctx->totorb++;
        if (!fileok) fc.stop(); 
        return ret;
    }
    ManySinks_Context *ctx;
};
struct ManySinks_ReaderOH { //uses readOrbitWithHeader
    ManySinks_ReaderOH(ManySinks_Context &context) : ctx(&context) {} //takes a pointer to context
    ManySinks_Context::Pair operator()(oneapi::tbb::flow_control& fc) const { //overloading "pair"
        ManySinks_Context::Pair ret(nullptr, 0);
        if (!ctx->fin.good()) { //stop everything if there is an error (can't do io operations)
            fc.stop();
            return ret;
        }
        ret = ctx->nextPair();
        assert(ret.first->empty());
        bool fileok = ret.first->readOrbitWithHeader(ctx->fin, ctx->totpuppi); //reading step
        ctx->totorb++; ctx->totev += BX_PER_ORBIT;
        if (!fileok) fc.stop(); 
        return ret;
    }
    ManySinks_Context *ctx;
};

//Merging readers
struct First_Merge_Reader {
    First_Merge_Reader(MergeSplit_Context &context, unsigned int nin) : ctx(&context), nin(nin) {}
    std::vector<uint> operator()(oneapi::tbb::flow_control& fc) const {
        std::vector<uint> p{0, 0, 0}; //return index of next file to write out to
        if (!ctx->fins[0].good()) {
            printf("there is some kind of error with the in file");
            fc.stop();
            return p;
        }
        //switch file to write out to and find new start of orbit
        p = ctx->nextQuad();
        int eventstoread = 3564 / nin + ((3564 % nin) > 0 ? 1 : 0);
        bool fileok = ctx->buffers[nin*p[2]]->readNEvents(ctx->fins[p[1]], eventstoread, ctx->totev, ctx->totpuppi, p[3]);
        if (!fileok) {
            fc.stop();
            //printf("file not ok \n"); fflush(stdout);
        }
        //printf("output file: %d, orbit start: %d, buffer batch: %d, orbit: %d \n", p[0], p[1], p[2], p[3]);
        return p;
    }
    MergeSplit_Context *ctx;
    unsigned int nin;
};

struct Merge_Reader {
    Merge_Reader() {};
    Merge_Reader(MergeSplit_Context &context, unsigned int index, unsigned int nin) : ctx(&context), i(index), nin(nin) {}
    std::vector<uint> operator()(std::vector<uint> p) const {
        /*if (!ctx->fins[i].good()) {
            fc.stop();
            return ret;
        }*/
        //fix this:
        //printf("got here \n"); fflush(stdout);
        int eventstoread = 3564 / nin + ((3564 % nin) > i ? 1 : 0);
        int fileidx = (p[1] + i) % nin; 
        //printf("will try to read from file %d \n", fileidx); fflush(stdout);
        bool fileok = ctx->buffers[nin*p[2]+i]->readNEvents(ctx->fins[fileidx], eventstoread, ctx->totev, ctx->totpuppi, p[3]);
        //printf("got here 2 \n"); fflush(stdout);
        if (!fileok) {
            printf("file not ok \n");
            fflush(stdout);
            exit(1);
        }
        //printf("read into buffer \n"); fflush(stdout);
        return p;
    }
    MergeSplit_Context *ctx;
    unsigned int i; //the index of the reader (and therefore the buffer)
    unsigned int nin;
    //unsigned int out;
};

//NEXT: make demux writer take nin
struct Demux_Writer {
    Demux_Writer() {};
    Demux_Writer(MergeSplit_Context &context, int index, int nin, std::fstream &fout) : 
        ctx(&context), self(index), fout(&fout), nin(nin) {}
    std::vector<uint> operator()(std::vector<uint> p) const {
        if (p[0] == self) {
            int nevents = 0;
            for (uint i = 0; i < 3564; ++i) {
                unsigned int ibuff = i % nin + nin * p[2]; //buffer to read from
                ctx->buffers[ibuff]->writeevent(*fout);
                //printf("event written: %d \n", nevents++, nin); fflush(stdout);
            }
            //printf("got here \n");
            for (uint i = 0; i < nin; ++i) { //after we're done writing out initialize buffers
                //printf("clearing buffer # %d \n", i);
                if (!ctx->buffers[i + nin*p[2]]->writeptrclear()) {
                    printf("failed to clear buffer # %d \n orbit = %d, bx = %d \n", i + nin*p[2], p[3], ctx->totev);
                    exit(1);
                } //else printf("succeeded in clearing buffer # %d", i); fflush(stdout);
            }
        }
    return p;
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};

struct Last_Demux_Writer { //the only difference is the constructor is VOID
    Last_Demux_Writer() {};
    Last_Demux_Writer(MergeSplit_Context &context, int index, std::fstream &fout, uint nin) : 
        ctx(&context), self(index), fout(&fout), nin(nin) {}
    void operator()(std::vector<uint> p) const {
        if (p[0] == self) {
            for (uint i = 0; i < 3564; ++i) {
                unsigned int ibuff = i % nin + nin*p[2]; //buffer to read from
                ctx->buffers[ibuff]->writeevent(*fout);
            }
            for (uint i = 0; i < nin; ++i) { //after we're done writing out initialize buffers
                ctx->buffers[i + nin*p[2]]->writeptrclear();
            }
        }
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};


struct ManySinks_Writer {
    ManySinks_Writer() {}
    ManySinks_Writer(ManySinks_Context &context, int index, std::fstream &out) : 
        ctx(&context), self(index), fout(&out) {} 
    ManySinks_Context::Pair operator()(ManySinks_Context::Pair data) const {
        if (data.second == self && data.first != nullptr) {
            if (ctx->write) data.first->writeToFile(*fout);
            data.first->clear();
        }
        return data;
    }
    
    ManySinks_Context *ctx;
    unsigned int self;
    std::fstream *fout;
};
struct ManySinks_LastWriter {
    ManySinks_LastWriter() {}
    ManySinks_LastWriter(ManySinks_Context &context, int index, std::fstream &out) : 
        ctx(&context), self(index), fout(&out) {}
    void operator()(ManySinks_Context::Pair data) const {
        if (data.second == self && data.first != nullptr) {
            if (ctx->write) data.first->writeToFile(*fout);
            data.first->clear();
        }
    }
    
    ManySinks_Context *ctx;
    unsigned int self;
    std::fstream *fout;
};

int manysinks_fullorbit_1toNscatter(const char *fname, int nouts, int bufsizeMb, const char *foutname) {
    assert(nouts >= 1);
    auto t0 = std::chrono::steady_clock::now();
    ManySinks_Context ctx(fname, nouts, bufsizeMb, foutname);
    ManySinks_Reader reader(ctx);
    std::vector<ManySinks_Writer> writers(nouts-1);
    for (unsigned int i = 0; i < nouts-1; ++i) {
        writers[i] = ManySinks_Writer(ctx, i, ctx.fouts[i]);
    }
    ManySinks_LastWriter lastwriter(ctx, nouts-1, ctx.fouts.back());
    oneapi::tbb::filter<void,ManySinks_Context::Pair> queue;
    queue = oneapi::tbb::make_filter<void,ManySinks_Context::Pair>(oneapi::tbb::filter_mode::serial_in_order, reader);
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<ManySinks_Context::Pair,ManySinks_Context::Pair>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    auto full_queue = queue & oneapi::tbb::make_filter<ManySinks_Context::Pair,void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(nouts+1,full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}
int manysinks_fullorbitOH_1toNscatter(const char *fname, int nouts, int bufsizeMb, const char *foutname) {
    assert(nouts >= 1);
    auto t0 = std::chrono::steady_clock::now();
    ManySinks_Context ctx(fname, nouts, bufsizeMb, foutname);
    ManySinks_ReaderOH reader(ctx);
    std::vector<ManySinks_Writer> writers(nouts-1);
    for (unsigned int i = 0; i < nouts-1; ++i) {
        writers[i] = ManySinks_Writer(ctx, i, ctx.fouts[i]);
    }
    ManySinks_LastWriter lastwriter(ctx, nouts-1, ctx.fouts.back());
    oneapi::tbb::filter<void,ManySinks_Context::Pair> queue;
    queue = oneapi::tbb::make_filter<void,ManySinks_Context::Pair>(oneapi::tbb::filter_mode::serial_in_order, reader);
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<ManySinks_Context::Pair,ManySinks_Context::Pair>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    auto full_queue = queue & oneapi::tbb::make_filter<ManySinks_Context::Pair,void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(nouts+1,full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}

int mergesplit(int opt, const char *finname, unsigned int nin, unsigned int nout, unsigned int bufsizeMb, const char *foutname) {
    assert (nin >= 1 && nout >= 1);
    auto t0 = std::chrono::steady_clock::now();
    MergeSplit_Context ctx(opt, finname, nin, nout, bufsizeMb, foutname);
    First_Merge_Reader firstreader(ctx, nin);
    oneapi::tbb::filter<void,std::vector<uint>> queue = oneapi::tbb::make_filter<void, std::vector<uint>>(oneapi::tbb::filter_mode::serial_in_order, firstreader);
    if (nin > 1) {
        std::vector<Merge_Reader> readers(nin-1);
        for (unsigned int i = 0; i < nin-1; ++i) {
            readers[i] = Merge_Reader(ctx, i+1, nin);
        }
    for (auto & r : readers) {
        queue = queue & oneapi::tbb::make_filter<std::vector<uint>, std::vector<uint>>(oneapi::tbb::filter_mode::serial_in_order, r);
    }
    }
    std::vector<Demux_Writer> writers(nout-1);
    for (unsigned int i = 0; i < nout-1; ++i) {
        writers[i] = Demux_Writer(ctx, i, nin, ctx.fouts[i]);
    }
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<std::vector<uint>, std::vector<uint>>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    Last_Demux_Writer lastwriter(ctx, nout-1, ctx.fouts[nout-1], nin);
    oneapi::tbb::filter<void,void> full_queue = queue & oneapi::tbb::make_filter<std::vector<uint>, void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(nin + nout + 1, full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}


int main(int argc, char**argv) {
    if (argc <= 2) {
        printf("Usage: %s test args\n", argv[0]);
        return 1;
    }
    std::string method(argv[1]);
    if (method == "sequential_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts outnames\n", argv[0], argv[1]);
        }
        return sequential_1toNscatter(argv[2], std::atoi(argv[3]), argc >= 5 ? argv[4] : "null");
    /*} else if (method == "sequential_orbit_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts bufsize_Mb outnames\n", argv[0], argv[1]);
        }
        return sequential_orbit_1toNscatter(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), argc >= 6 ? argv[5] : "null");*/
    } else if (method == "sequential_fullorbit_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts bufsize_Mb outnames\n", argv[0], argv[1]);
        }
        return sequential_fullorbit_1toNscatter(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), argc >= 6 ? argv[5] : "null");
    } else if (method == "sequential_fullorbitOH_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts bufsize_Mb outnames\n", argv[0], argv[1]);
        }
        return sequential_fullorbitOH_1toNscatter(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), argc >= 6 ? argv[5] : "null");
   } else if (method == "manysinks_fullorbit_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts bufsize_Mb outnames\n", argv[0], argv[1]);
        }
        return manysinks_fullorbit_1toNscatter(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), argc >= 6 ? argv[5] : "null");
    } else if (method == "manysinks_fullorbitOH_1toNscatter") {
        if (argc <= 2) {
            printf("Usage: %s %s nouts bufsize_Mb outnames\n", argv[0], argv[1]);
        }
        return manysinks_fullorbitOH_1toNscatter(argv[2], std::atoi(argv[3]), std::atoi(argv[4]), argc >= 6 ? argv[5] : "null");
     } else if (method == "mergesplit") {
         if (argc <= 2) {
            printf("Usage: %s %s nouts outnames\n", argv[0], argv[1]);
        }
        return mergesplit(std::atoi(argv[2]), argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]), argv[7]);
     } else {
        printf("Unknown test %s\n", method.c_str());
        return 2;
    }
    return 0;
}
