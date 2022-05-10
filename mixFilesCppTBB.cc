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
            snprintf(fnamebuff,1023, "%s%d.dump",foutname,i); // print file names to a buffer
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
        //function write to file writes to output file
            assert(bufsizeMb >= 1);
            size_t size = bufsizeMb*1024*1024;
            mem_ = reinterpret_cast<char *>(std::aligned_alloc(4096u, size)); //this is assigning the memory for the buffer
            //I guess we want to do things so explicitly in this case b/c we're dealing with a lot of memory
            //and later will deallocate?
            assert(mem_ != nullptr);
            ptr_ = mem_;
            end_ = mem_ + size;
        }
        ~OrbitBuffer() { //deallocates the memory (tilde defines a destructor)
            std::free(mem_);
        }
        void clear() { //just means return to the beginning of the buffer
            ptr_ = mem_;
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
            //reads a chunk of length and moves pointer
            in.read(ptr_, length); //read from in file
            assert(in.gcount() == length); //check that we did read length bytes
            ptr_ += length; //move the address to read into next time
            assert(ptr_ < end_); //check we're not at the end of allocated mem
        }
        size_t dataSize() const { return ptr_ - mem_; } //how much we have in
        void writeToFile(std::fstream & out) {
            out.write(mem_, dataSize());
        }
        //reads events one-by-one
        bool readNEvents(std::fstream & fin, 
                         unsigned int n, 
                         unsigned long long int & totev, //why would you want to use reference to your args?
                         unsigned long long int & totpuppi) {
            //I guess totev and totpuppi are just used for report
            uint64_t header, npuppi;
            uint32_t myorbit, orbit;
            for (unsigned int i = 0; fin.good() && i < n; ++i) {
                fin.read(reinterpret_cast<char*>(&header), sizeof(uint64_t)); //read header from file into "header"
                if (fin.gcount() == 0) break; //break if we read nothing
                //getting orbit (&checking), number of puppi, etc from file
                orbit = (header >> 24) & 0x3FFFFFFF;
                if (i == 0) myorbit = orbit;
                else assert(orbit == myorbit);
                write64(header);
                //what are totev and totpuppi used for?
                totev++;
                npuppi = header & 0xFFF;
                if (npuppi) {
                    readChunk(fin, npuppi*sizeof(uint64_t)); //read
                    totpuppi += npuppi;
                }
            }
            return fin.good();
        }
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
        char *mem_, *ptr_, *end_;
};

/*
int sequential_orbit_1toNscatter(const char *fname, int nouts, int bufsizeMb, const char *foutname) {
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
    uint64_t header, npuppi;
    uint32_t orbit, lastorbit = 1<<30;
    OrbitBuffer buffer(bufsizeMb);
    int iout = -1;
    auto t1 = std::chrono::steady_clock::now();    
    while(fin.good()) {
        fin.read(reinterpret_cast<char*>(&header), sizeof(uint64_t));
        if (fin.gcount() == 0) break;
        orbit  = (header >> 24) & 0x3FFFFFFF;
        if (orbit != lastorbit) {
            if (iout != -1) {
               if (write) buffer.writeToFile(fouts[iout]);
               buffer.clear();
            }
            totorb++; 
            iout = (iout + 1) % nouts;
            lastorbit = orbit;
        }
        buffer.write64(header);
        npuppi = header & 0xFFF;
        totev++;
        if (npuppi) {
            buffer.readChunk(fin, npuppi*sizeof(uint64_t));
            totpuppi += npuppi;
        }
    }
    printf("I still have to write out %llu data on buffer %d\n", buffer.dataSize(), iout);
    if (write) buffer.writeToFile(fouts[iout]);
    if (write) {
        for (auto & f : fouts) f.close();
    }
    auto t2 = std::chrono::steady_clock::now();
    report(totorb, totev, totpuppi, t0, t1, t2);
    return 0;
}
*/

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
            snprintf(fnamebuff,1023, "%s%d.dump",foutname,i);
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
        fileok = buffer.readNEvents(fin, nEvts, totev, totpuppi);
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
            snprintf(fnamebuff,1023, "%s%d.dump",foutname,i);
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
            nbuffs(nouts+2), //nouts+2
            ibuff(0), iout(0)
        {
            printf("Will distribute data from %s to %d sinks with one %dMB orbit buffer:\n", fname, nouts, bufsizeMb);
            if (std::string(foutname) == "null") {
                write = false;
            } else {
                //make files
                char fnamebuff[1024];
                for (int i = 0; i < nouts; ++i) {
                    snprintf(fnamebuff,1023, "%s%d.dump",foutname,i);
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
            while (!(ibuff.compare_exchange_weak(curr,(curr+1)%nbuffs))) ;
            //add 1 to ibuff or go back to 0
            OrbitBuffer *ret = buffers[curr];
            //assert(ret->empty());
            return ret;
        }
        int nextOut() {
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

struct ManySinks_Reader {
    //first moves to the next buffer
    //then reads an orbit into the buffer
    ManySinks_Reader(ManySinks_Context &context) : ctx(&context) {} //takes a pointer to context, which is now called ctx (?)
    ManySinks_Context::Pair operator()(oneapi::tbb::flow_control& fc) const {
        ManySinks_Context::Pair ret(nullptr, 0);
        if (!ctx->fin.good()) {
            fc.stop();
            return ret;
        }
        ret = ctx->nextPair();
        assert(ret.first->empty());
        bool fileok = ret.first->readNEvents(ctx->fin, BX_PER_ORBIT, ctx->totev, ctx->totpuppi);
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
struct ManySinks_Writer {
    //if the buffer matches the index of the file to write out to, writes out to file
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
     } else {
        printf("Unknown test %s\n", method.c_str());
        return 2;
    }
    return 0;
}
