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

class MuxBuffer { //group of nin buffers
    public:
        MuxBuffer(const uint nin) : n(nin) {
            mem.resize(n); ptr.resize(n); end.resize(n); writeptr.resize(n); events_per_buffer.resize(n);
            size_t size = 50 * 8 * 3564; //maximum size of events in buffer
            char* memtot = reinterpret_cast<char *>(std::aligned_alloc(4096u, size));
            char* creator = memtot;
            for (uint i = 0; i < n; ++i) { //make nin buffers
                events_per_buffer[i] = (BX_PER_ORBIT / nin) + ((i < (BX_PER_ORBIT % n)) ? 1 : 0);
                size_t size_i = events_per_buffer[i] * 50 * 8;
                mem[i] = creator;
                ptr[i] = mem[i];
                assert(mem[i] != nullptr);
                creator += size_i;
                end[i] = creator;
                writeptr[i] = mem[i];
            }
        }
        ~MuxBuffer() {
            for (uint i = 0; i < n; ++i) {
                std::free(mem[i]);
            }
        }
        void cleartest() {
            for (uint i = 0; i < n; ++i) {
                ptr[i] = mem[i];
                writeptr[i] = mem[i];
            }
        }
        bool clear() {
            for (uint i = 0; i < n; ++i) {
                if (ptr[i] != writeptr[i]) {
                    printf("ERROR: attempted to clear buffer before writing finished \n"); fflush(stdout);
                    printf("buffer index = %u, ptr - mem = %d, writeptr - mem = %d \n", i, ptr[i] - mem[i], writeptr[i] - mem[i]);
                    return false;
                }
                ptr[i] = mem[i];
                writeptr[i] = mem[i];
            }
            return true;
        }
        bool empty() const {
            for (uint i = 0; i < n; ++i) {
                if (ptr[i] != mem[i]) return false;
            }
            return true;
        }
        //the reading methods below work with 1/nin buffers at a time
        void write64(uint64_t word, uint i) {
            uint64_t *ptr64 = reinterpret_cast<uint64_t *>(ptr[i]);
            *ptr64 = word;
            ptr[i] += 8; 
            assert(ptr[i] < end[i]);
        }
        void readChunk(std::fstream & in, size_t length, uint i) {
            in.read(ptr[i], length); //read from in file
            assert(in.gcount() == length); //check that we did read length bytes
            ptr[i] += length; //move the address to read into next time
            assert(ptr[i] < end[i]); //check we're not at the end of allocated mem
        }
        size_t dataSize(uint i) const { return ptr[i] - mem[i]; } //how much we have in
        void writeToFile(std::fstream & out, uint i) {
            out.write(mem[i], dataSize(i));
        }
        void writeevent(std::fstream & out, uint i) { //writes one event from the buffer into a file
            uint64_t header;
            memcpy(&header, writeptr[i], 8);
            uint32_t npuppi = header & 0xFFF;
            size_t length = 8 * (npuppi + 1);
            out.write(writeptr[i], length);
            writeptr[i] += length;
        }
        void writedemux(std::fstream &out) {
            for (uint i = 0; i < BX_PER_ORBIT; ++i) {
                uint index = i % n;
                writeevent(out, index);
            }
        }
        void writedemux_OH(std::fstream &out) {
            for (uint i = 0; i < n; ++i) {
                out.write(writeptr[i], 8);
                writeptr[i] +=8;
            }
            for (uint i = 0; i < BX_PER_ORBIT; ++i) {
                uint index = i % n;
                writeevent(out, index);
            }
        }
        bool readOrbit(std::fstream & fin, 
                         unsigned int index,
                         unsigned long long int & totev, //reference lets you modify (this function can change totev)
                         unsigned long long int & totpuppi,
                         unsigned int myorbit) {
            uint64_t header, npuppi;
            uint32_t orbit;
            for (unsigned int i = 0; fin.good() && i < events_per_buffer[index]; ++i) {
                fin.read(reinterpret_cast<char*>(&header), 8); //read header from file into "header" 
                if (fin.gcount() == 0) break; //break if we read nothing
                orbit = (header >> 24) & 0x3FFFFFFF;
                npuppi = header & 0xFFF;
                uint32_t bx = (header >> 12) & 0xFFF;
                //printf("orbit: %lu \n", orbit);
                //printf("bx: %lu \n", bx);
                //printf("npuppi: %llu \n", npuppi);
                /*if (i == 0) myorbit = orbit;
                else */ if (orbit != myorbit) {
                    printf("ERROR: orbit from header: %lu, expected orbit: %u \n", orbit, myorbit);
                    printf("bunch crossings: %llu, i = %d \n", totev, i);
                    exit(1);
                }
                write64(header, index);
                totev++;
                if (npuppi) {
                    readChunk(fin, npuppi*sizeof(uint64_t), index); //read
                    totpuppi += npuppi;
                }
            }
            return fin.good();
        }
        bool readOrbitWithHeader(std::fstream & fin, uint i, unsigned long long int & totpuppi) {
            uint64_t orbitHeader, orbit;
            fin.read(reinterpret_cast<char*>(&orbitHeader), sizeof(uint64_t)); //read orbit header into orbitHeader
            if (fin.gcount() == 0) return false;
            orbit  = (orbitHeader >> 24) & 0x3FFFFFFF;
            size_t ndata = orbitHeader & ((1 << 24)-1);
            size_t length = ndata * sizeof(uint64_t);
            //printf("orbit %llu of length %llu\n", orbit, ndata);
            write64(orbitHeader, i); //write orbit header into buffer
            fin.read(ptr[i], length);
            if (fin.gcount() != length) {
                printf("ERROR: tried to write %d, but only wrote %d \n", length, fin.gcount());
                printf("%d \n", ndata);
                if (!fin.good()) printf("file not good \n"); else printf("file still good \n");
                printf("orbit = %d \n", orbit);
                return false;
            }
            ptr[i] += length;
            totpuppi += ndata;
            return fin.good();
        }
    public:
        std::vector<uint> events_per_buffer;
    private:
        const uint n; //number of in files = number of buffers in each muxbuffer object
        std::vector<char *> mem, ptr, end, writeptr;
};

class AssemblyBuffer { //like MuxBuffer, but also includes a buffer to assemble demultiplexed orbit in
    public:
        AssemblyBuffer(const uint nin) : n(nin) {
            mem.resize(n); ptr.resize(n); end.resize(n); writeptr.resize(n); events_per_buffer.resize(n);
            size_t size = 50 * 8 * 3564; //maximum size of events in buffer
            //make in buffers
            char* memtot = reinterpret_cast<char *>(std::aligned_alloc(4096u, size));
            char* creator = memtot;
            for (uint i = 0; i < n; ++i) { //make nin buffers
                events_per_buffer[i] = (BX_PER_ORBIT / nin) + ((i < (BX_PER_ORBIT % n)) ? 1 : 0);
                size_t size_i = events_per_buffer[i] * 50 * 8;
                mem[i] = creator;
                ptr[i] = mem[i];
                assert(mem[i] != nullptr);
                creator += size_i;
                end[i] = creator;
                writeptr[i] = mem[i];
            }
            //make out buffer
            omem = reinterpret_cast<char *>(std::aligned_alloc(4096u, size));
            optr = omem;
            owriteptr = omem;
            oend = omem + size;
        }
        ~AssemblyBuffer() {
            for (uint i = 0; i < n; ++i) {
                std::free(mem[i]);
            }
            std::free(omem);
        }
        void cleartest() {
            for (uint i = 0; i < n; ++i) {
                ptr[i] = mem[i];
                writeptr[i] = mem[i];
            }
        }
        bool clearin() {
            for (uint i = 0; i < n; ++i) {
                if (ptr[i] != writeptr[i]) {
                    printf("ERROR: attempted to clear in buffer before writing finished \n"); fflush(stdout);
                    printf("buffer index = %u, ptr - mem = %d, writeptr - mem = %d \n", i, ptr[i] - mem[i], writeptr[i] - mem[i]);
                    return false;
                }
                ptr[i] = mem[i];
                writeptr[i] = mem[i];
            }
            printf("cleared in buffer successfully \n");
            return true;
        }
        void clearout() {
            if (optr != owriteptr) {
                printf("ERROR: attempted to clear out buffer before writing finished \n");
                exit(1);
            }
            optr = omem;
            owriteptr = omem;
        }
        bool inempty() const {
            for (uint i = 0; i < n; ++i) {
                if (ptr[i] != mem[i]) return false;
            }
            return true;
        }
        bool outempty() const {
            return (optr == omem);
        }
        //the reading methods below work with 1/nin buffers at a time
        void write64(uint64_t word, uint i) {
            uint64_t *ptr64 = reinterpret_cast<uint64_t *>(ptr[i]);
            *ptr64 = word;
            ptr[i] += 8; 
            assert(ptr[i] < end[i]);
        }
        void readinChunk(std::fstream & in, size_t length, uint i) {
            in.read(ptr[i], length); //read from in file
            assert(in.gcount() == length); //check that we did read length bytes
            ptr[i] += length; //move the address to read into next time
            assert(ptr[i] < end[i]); //check we're not at the end of allocated mem
        }
        size_t dataSize(uint i) const { return ptr[i] - mem[i]; } //how much we have in
       /* void writeToFile(std::fstream & out, uint i) {
            out.write(mem[i], dataSize(i));
        } */
        void writeevent(uint i) { //writes one event from in buffer to out buffer
            uint64_t header;
            memcpy(&header, writeptr[i], 8);
            uint32_t npuppi = header & 0xFFF;
            size_t length = 8 * (npuppi + 1);
            memcpy(optr, writeptr[i], length);
            writeptr[i] += length;
            optr += length;
        }
        void bufferdemux() {
            for (uint i = 0; i < BX_PER_ORBIT; ++i) {
                uint index = i % n;
                writeevent(index);
            }
        }
        void bufferdemux_OH() { //for now we just copy in header
            for (uint i = 0; i < n; ++i) {
                memcpy(optr,writeptr[i], 8);
                writeptr[i] += 8;
                optr += 8;
            }
            for (uint i = 0; i < BX_PER_ORBIT; ++i) {
                uint index = i % n;
                writeevent(index);
            }
        }
        void writetofile(std::fstream & out) {
            size_t size = optr - omem;
            out.write(optr, size);
            owriteptr += size;
            assert(owriteptr == optr);
        }
        bool readOrbit(std::fstream & fin, 
                         unsigned int index,
                         unsigned long long int & totev, //reference lets you modify (this function can change totev)
                         unsigned long long int & totpuppi,
                         unsigned int myorbit) {
            uint64_t header, npuppi;
            uint32_t orbit;
            for (unsigned int i = 0; fin.good() && i < events_per_buffer[index]; ++i) {
                fin.read(reinterpret_cast<char*>(&header), 8); //read header from file into "header" 
                if (fin.gcount() == 0) break; //break if we read nothing
                orbit = (header >> 24) & 0x3FFFFFFF;
                npuppi = header & 0xFFF;
                uint32_t bx = (header >> 12) & 0xFFF;
                //printf("orbit: %lu \n", orbit);
                //printf("bx: %lu \n", bx);
                //printf("npuppi: %llu \n", npuppi);
                /*if (i == 0) myorbit = orbit;
                else */ if (orbit != myorbit) {
                    printf("ERROR: orbit from header: %lu, expected orbit: %u \n", orbit, myorbit);
                    printf("bunch crossings: %llu, i = %d \n", totev, i);
                    exit(1);
                }
                write64(header, index);
                totev++;
                if (npuppi) {
                    readinChunk(fin, npuppi*sizeof(uint64_t), index); //read
                    totpuppi += npuppi;
                }
            }
            return fin.good();
        }
        bool readOrbitWithHeader(std::fstream & fin, uint i, unsigned long long int & totpuppi) {
            uint64_t orbitHeader, orbit;
            fin.read(reinterpret_cast<char*>(&orbitHeader), sizeof(uint64_t)); //read orbit header into orbitHeader
            if (fin.gcount() == 0) return false;
            orbit  = (orbitHeader >> 24) & 0x3FFFFFFF;
            size_t ndata = orbitHeader & ((1 << 24)-1);
            size_t length = ndata * sizeof(uint64_t);
            //printf("orbit %llu of length %llu\n", orbit, ndata);
            write64(orbitHeader, i); //write orbit header into buffer
            fin.read(ptr[i], length);
            if (fin.gcount() != length) {
                printf("ERROR: tried to write %d, but only wrote %d \n", length, fin.gcount());
                printf("%d \n", ndata);
                if (!fin.good()) printf("file not good \n"); else printf("file still good \n");
                printf("orbit = %d \n", orbit);
                return false;
            }
            ptr[i] += length;
            totpuppi += ndata;
            return fin.good();
        }
    public:
        std::vector<uint> events_per_buffer;
    private:
        const uint n; //number of in files = number of buffers in each muxbuffer object
        std::vector<char *> mem, ptr, end, writeptr;
        char *omem, *optr, *oend, *owriteptr; //o for out
};

class MergeSplit_Context {
    public:
        MergeSplit_Context(const char* finname, const char* foutname, uint nameopt, uint nin, uint nout, bool assembly) :
        nout(nout),
        nin(nin),
        write(true),
        nbuff(nout + nin + 1),
        foutidx(0),
        bufidx(0),
        orbitstart(0),
        totorb(0), totev(0), totpuppi(0)
        {
            printf("Will distribute data from %d to %d files with %d %dMB buffers.\n", nin, nout, nbuff, 7.33);
            if (std::string(foutname) == "null" or std::string(finname) == "null") {
                write = false;
            } else {
                char fnamebuff[1024];
                fins.resize(nin);
                if (nameopt == 0) {
                //make files
                //option 0: format *n.dump where n is an index 0 to N-1)
                printf("Input files:\n");
                for (int i = 0; i < nin; ++i) {
                    snprintf(fnamebuff,1023, "%s%d.dump",finname,i);
                    fins[i].open(fnamebuff, std::ios_base::in | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                } } else if (nameopt == 1) {
                //option 1: just one input file
                    if (nin != 1) {
                        printf("ERROR: Method 1 takes only one file \n");
                        exit(1);
                    }
                    fins[0].open(finname, std::ios_base::in | std::ios_base::binary);
                } else if (nameopt == 2) {
                //option 2: format *.bx_n_of_m.dump
                    printf("Input files:\n");
                for (int i = 0; i < nin; ++i) {
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
                for (int i = 0; i < nout; ++i) {
                    snprintf(fnamebuff,1023, "%s%d.dump",foutname,i);
                    fouts[i].open(fnamebuff, std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
                    printf(" - %s\n", fnamebuff);
                }
                }
                printf("\n");
            //make the buffers
            if (assembly) {
                abuffers.resize(nbuff);
                for (uint i = 0; i < nbuff; ++i) {
                    abuffers[i] = new AssemblyBuffer(nin);
                }
            } else {
                buffers.resize(nbuff);
                for (int i = 0; i < nbuff; ++i) {
                    buffers[i] = new MuxBuffer(nin);
                    //printf("made one buffer \n"); fflush(stdout);
                }
            }
            //printf("finished making buffers \n"); fflush(stdout);
        }
        MuxBuffer* nextBuffer() {
            //returns buffer to write into this time
            //and increments the buffer for next time
            uint curbuf, nextbuf;
            while(true) {
                curbuf = bufidx;
                nextbuf = (curbuf + 1) % nbuff;
                if (bufidx.compare_exchange_weak(curbuf, nextbuf)) {
                    return buffers[curbuf];
                }
            }
        }
        uint nextABuffer() {
            uint curbuf, nextbuf;
            while(true) {
                curbuf = bufidx;
                nextbuf = (curbuf + 1) % nbuff;
                if (bufidx.compare_exchange_weak(curbuf, nextbuf)) {
                    return curbuf;
                }
            }
        }
        uint nextOut() {
            uint curout, nextout;
            while(true) {
                curout = foutidx;
                nextout = (curout + 1) % nout;
                if (foutidx.compare_exchange_weak(curout, nextout)) {
                    return curout;
                }
            }
        }
        uint nextOrbitStart() { //this only necessary in the case where 3564 % nin != 0 which will not happen in real life
            uint curstart, nextstart;
            while (true) {
                curstart = orbitstart;
                nextstart = (curstart + BX_PER_ORBIT % nin) % nin;
                if (orbitstart.compare_exchange_weak(curstart, nextstart)) {
                    return curstart;
                }
            }
        }
        uint nextOrbit() {
            uint curorbit, nextorbit;
            while (true) {
                curorbit = totorb;
                nextorbit = totorb + 1;
                if (totorb.compare_exchange_weak(curorbit, nextorbit)) {
                    return curorbit;
                }
            }
        }
        class TokInfo {
            public:
                TokInfo(MuxBuffer* buffer, uint file, uint orbitstart, uint orbit) :
                    buffer(buffer),
                    file(file),
                    orbitstart(orbitstart),
                    orbit(orbit)
                    {}
                TokInfo (std::vector<AssemblyBuffer *> abuff, uint ABufIdx, uint file, uint orbitstart, uint orbit) :
                    abuffer(abuff[ABufIdx]),
                    bufidx(ABufIdx),
                    file(file),
                    orbitstart(orbitstart),
                    orbit(orbit)
                    {}
            public:
                MuxBuffer* buffer;
                AssemblyBuffer* abuffer;
                uint bufidx;
                uint file;
                uint orbitstart;
                uint orbit;
        };
        TokInfo nextTok() {
            return TokInfo(nextBuffer(), nextOut(), nextOrbitStart(), nextOrbit());
        }
        TokInfo nextATok() {
            return TokInfo(abuffers, nextABuffer(), nextOut(), nextOrbitStart(), nextOrbit());
        }
    public:
        std::vector<std::fstream> fins, fouts;
        uint nin, nout;
        unsigned long long int totev, totpuppi;
        std::atomic<uint> totorb;
        std::vector<MuxBuffer *> buffers;
        std::vector<AssemblyBuffer *> abuffers;
        bool write;
    private:
        uint nbuff;
        std::atomic<uint> foutidx, bufidx, orbitstart;
};

struct First_Reader {
    First_Reader(MergeSplit_Context &context, unsigned int nin) : ctx(&context), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(oneapi::tbb::flow_control& fc) const {
        MergeSplit_Context::TokInfo t = ctx->nextTok();
        for (uint i = 0; i < nin; ++i) {
            if (!ctx->fins[0].good()) {
                printf("there is some kind of error with the in file");
                fc.stop();
                return t;
            }
        }
        if (!t.buffer->empty()) {
            printf("ERROR: attempted to read into buffer before empty \n");
            printf("totorb = %d, totev = %d \n", t.orbit, ctx->totev);
        }
        bool fileok = t.buffer->readOrbit(ctx->fins[t.orbitstart], 0, ctx->totev, ctx->totpuppi, t.orbit);
        //printf("got here 1 \n"); fflush(stdout);
        if (!fileok) {
            fc.stop();
            //printf("file not ok \n"); fflush(stdout);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int nin;
};

struct First_Assembly_Reader {
    First_Assembly_Reader(MergeSplit_Context &context, unsigned int nin) : ctx(&context), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(oneapi::tbb::flow_control& fc) const {
        MergeSplit_Context::TokInfo t = ctx->nextATok();
        for (uint i = 0; i < nin; ++i) {
            if (!ctx->fins[0].good()) {
                printf("there is some kind of error with the in file");
                fc.stop();
                return t;
            }
        }
        if (!t.abuffer->inempty()) {
            printf("ERROR: attempted to read into buffer before empty \n");
            printf("totorb = %d, totev = %d \n", t.orbit, ctx->totev);
        }
        bool fileok = t.abuffer->readOrbit(ctx->fins[t.orbitstart], 0, ctx->totev, ctx->totpuppi, t.orbit);
        //printf("got here 1 \n"); fflush(stdout);
        if (!fileok) {
            fc.stop();
            //printf("file not ok \n"); fflush(stdout);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int nin;
};

struct Merge_Reader {
    Merge_Reader() {};
    Merge_Reader(MergeSplit_Context &context, unsigned int index, unsigned int nin) : ctx(&context), i(index), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        bool fileok = t.buffer->readOrbit(ctx->fins[t.orbitstart + i], i, ctx->totev, ctx->totpuppi, t.orbit);
        if (!fileok) {
            printf("file not ok \n");
            fflush(stdout);
            exit(1);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int i; //the index of the reader (and therefore the buffer)
    unsigned int nin;
};

struct Assembly_Reader {
    Assembly_Reader() {};
    Assembly_Reader(MergeSplit_Context &context, unsigned int index, unsigned int nin) : ctx(&context), i(index), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        bool fileok = t.abuffer->readOrbit(ctx->fins[t.orbitstart + i], i, ctx->totev, ctx->totpuppi, t.orbit);
        if (!fileok) {
            printf("file not ok \n");
            fflush(stdout);
            exit(1);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int i; //the index of the reader (and therefore the buffer)
    unsigned int nin;
};

struct First_Reader_OH {
    First_Reader_OH(MergeSplit_Context &context, unsigned int nin) : ctx(&context), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(oneapi::tbb::flow_control& fc) const {
        MergeSplit_Context::TokInfo t = ctx->nextTok();
        for (uint i = 0; i < nin; ++i) {
            if (!ctx->fins[0].good()) {
                printf("there is some kind of error with the in file");
                fc.stop();
                return t;
            }
        }
        if (!t.buffer->empty()) {
            printf("ERROR: attempted to read into buffer before empty \n");
            printf("totorb = %d, totev = %d \n", t.orbit, ctx->totev);
        }
        bool fileok = t.buffer->readOrbitWithHeader(ctx->fins[t.orbitstart], 0, ctx->totpuppi);
        if (!fileok) {
            fc.stop();
        }
        ctx->totev += t.buffer->events_per_buffer[0];
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int nin;
};

struct First_Assembly_Reader_OH {
    First_Assembly_Reader_OH(MergeSplit_Context &context, unsigned int nin) : ctx(&context), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(oneapi::tbb::flow_control& fc) const {
        MergeSplit_Context::TokInfo t = ctx->nextATok();
        for (uint i = 0; i < nin; ++i) {
            if (!ctx->fins[0].good()) {
                printf("there is some kind of error with the in file");
                fc.stop();
                return t;
            }
        }
        if (!t.abuffer->inempty()) {
            printf("ERROR: attempted to read into buffer before empty \n");
            printf("totorb = %d, totev = %d \n", t.orbit, ctx->totev);
        }
        bool fileok = t.abuffer->readOrbitWithHeader(ctx->fins[t.orbitstart], 0, ctx->totpuppi);
        if (!fileok) {
            fc.stop();
        }
        ctx->totev += t.abuffer->events_per_buffer[0];
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int nin;
};

struct Merge_Reader_OH {
    Merge_Reader_OH() {};
    Merge_Reader_OH(MergeSplit_Context &context, unsigned int index, unsigned int nin) : ctx(&context), i(index), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        bool fileok = t.buffer->readOrbitWithHeader(ctx->fins[t.orbitstart + i], i, ctx->totpuppi);
        ctx->totev += t.buffer->events_per_buffer[i];
        if (!fileok) {
            printf("file not ok \n");
            fflush(stdout);
            exit(1);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int i; //the index of the reader (and therefore the buffer)
    unsigned int nin;
};

struct Assembly_Reader_OH {
    Assembly_Reader_OH() {};
    Assembly_Reader_OH(MergeSplit_Context &context, unsigned int index, unsigned int nin) : ctx(&context), i(index), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        bool fileok = t.abuffer->readOrbitWithHeader(ctx->fins[t.orbitstart + i], i, ctx->totpuppi);
        ctx->totev += t.abuffer->events_per_buffer[i];
        if (!fileok) {
            printf("file not ok \n");
            fflush(stdout);
            exit(1);
        }
        return t;
    }
    MergeSplit_Context *ctx;
    unsigned int i; //the index of the reader (and therefore the buffer)
    unsigned int nin;
};

struct Assembler {
    Assembler() {};
    Assembler(MergeSplit_Context &context, uint index) : self(index), ctx(&context) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        if (t.bufidx == self) {
            t.abuffer->bufferdemux();
            bool ok = t.abuffer->clearin();
            if (!ok) {
                printf("orbit: %d, total events: %d, buffer group: %d \n", t.orbit, ctx->totev, t.bufidx);
                exit(1);
            }
        }
        return t;
    }
    uint self;
    MergeSplit_Context *ctx;
};

struct Assembler_OH {
    Assembler_OH() {};
    Assembler_OH(MergeSplit_Context &context, uint index) : self(index), ctx(&context) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        if (t.bufidx == self) {
            t.abuffer->bufferdemux_OH();
            bool ok = t.abuffer->clearin();
            if (!ok) {
                printf("orbit: %d, total events: %d, buffer group: %d \n", t.orbit, ctx->totev, t.bufidx);
                exit(1);
            }
        }
        return t;
    }
    uint self;
    MergeSplit_Context *ctx;
};

struct Demux_Writer {
    Demux_Writer() {};
    Demux_Writer(MergeSplit_Context &context, int index, int nin, std::fstream &fout) : 
        ctx(&context), self(index), fout(&fout), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.buffer->writedemux(*fout);
            bool ok = t.buffer->clear();
            if (!ok) {
                printf("orbit: %d, totev: %d \n", t.orbit, ctx->totev);
                exit(1);
            }
        }
    return t;
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};

struct Demux_Writer_OH {
    Demux_Writer_OH() {};
    Demux_Writer_OH(MergeSplit_Context &context, int index, int nin, std::fstream &fout) : 
        ctx(&context), self(index), fout(&fout), nin(nin) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.buffer->writedemux_OH(*fout);
            bool ok = t.buffer->clear();
            if (!ok) {
                printf("orbit: %d, totev: %d \n", t.orbit, ctx->totev);
                exit(1);
            }
        }
    return t;
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};

struct Last_Demux_Writer {
    Last_Demux_Writer() {};
    Last_Demux_Writer(MergeSplit_Context &context, int nin, std::fstream &fout) : 
        ctx(&context), self(0), fout(&fout), nin(nin) {}
    void operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.buffer->writedemux(*fout);
            bool ok = t.buffer->clear();
            if (!ok) {
                printf("orbit: %d, totev: %d \n", t.orbit, ctx->totev);
                exit(1);
            }
        }
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};

struct Last_Demux_Writer_OH {
    Last_Demux_Writer_OH() {};
    Last_Demux_Writer_OH(MergeSplit_Context &context, int nin, std::fstream &fout) : 
        ctx(&context), self(nin-1), fout(&fout), nin(nin) {}
    void operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.buffer->writedemux_OH(*fout);
            bool ok = t.buffer->clear();
            if (!ok) {
                printf("orbit: %d, totev: %d \n", t.orbit, ctx->totev);
                exit(1);
            }
        }
    }
    MergeSplit_Context *ctx;
    unsigned int self;
    std::fstream *fout; //= ctx->fouts[self];
    unsigned int nin;
};

struct Assembled_Writer {
    Assembled_Writer() {};
    Assembled_Writer(uint index, uint nin, std::fstream &fout) :
        self(index), nin(nin), fout(&fout) {}
    MergeSplit_Context::TokInfo operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.abuffer->writetofile(*fout);
            t.abuffer->clearout();
        }
    return t;
    }
    uint nin, self;
    std::fstream *fout;
};

struct Last_Assembled_Writer {
    Last_Assembled_Writer() {};
    Last_Assembled_Writer(uint nin, std::fstream &fout) :
        self(nin), fout(&fout) {}
    void operator()(MergeSplit_Context::TokInfo t) const {
        if (t.file == self) {
            t.abuffer->writetofile(*fout);
            t.abuffer->clearout();
        }
    }
    uint self;
    std::fstream *fout;
};

int mergesplit(const char* finname, const char* foutname, uint nameopt, uint nin, uint nout) {
    assert(nin > 0 && nout > 0);
    auto t0 = std::chrono::steady_clock::now();
    MergeSplit_Context ctx(finname, foutname, nameopt, nin, nout, false);
    First_Reader firstreader(ctx, nin);
    oneapi::tbb::filter<void,MergeSplit_Context::TokInfo> queue = oneapi::tbb::make_filter<void, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, firstreader);
    if (nin > 1) {
        std::vector<Merge_Reader> readers(nin-1);
        for (unsigned int i = 0; i < nin-1; ++i) {
            readers[i] = Merge_Reader(ctx, i+1, nin);
        }
    for (auto & r : readers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, r);
    }
    }
    std::vector<Demux_Writer> writers(nout-1);
    for (unsigned int ix = 0; ix < nout-1; ++ix) {
        uint i = nout - ix;
        writers[i] = Demux_Writer(ctx, i, nin, ctx.fouts[i]);
    }
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    Last_Demux_Writer lastwriter(ctx, 0, ctx.fouts[0]);
    oneapi::tbb::filter<void,void> full_queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(nin + nout + 1, full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}

int mergesplit_OH(const char* finname, const char* foutname, uint nameopt, uint nin, uint nout) {
    assert(nin > 0 && nout > 0);
    auto t0 = std::chrono::steady_clock::now();
    MergeSplit_Context ctx(finname, foutname, nameopt, nin, nout, false);
    First_Reader_OH firstreader(ctx, nin);
    oneapi::tbb::filter<void,MergeSplit_Context::TokInfo> queue = oneapi::tbb::make_filter<void, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, firstreader);
    if (nin > 1) {
        std::vector<Merge_Reader_OH> readers(nin-1);
        for (unsigned int i = 0; i < nin-1; ++i) {
            readers[i] = Merge_Reader_OH(ctx, i+1, nin);
        }
    for (auto & r : readers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, r);
    }
    }
    std::vector<Demux_Writer_OH> writers(nout-1);
    for (unsigned int i = 0; i < nout-1; ++i) {
        writers[i] = Demux_Writer_OH(ctx, i, nin, ctx.fouts[i]);
    }
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    Last_Demux_Writer_OH lastwriter(ctx, nin, ctx.fouts[nout-1]);
    oneapi::tbb::filter<void,void> full_queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(1 + nin + nout, full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}

int assembleinbuffer(const char* finname, const char* foutname, uint nameopt, uint nin, uint nout) {
    assert(nin > 0 && nout > 0);
    auto t0 = std::chrono::steady_clock::now();
    MergeSplit_Context ctx(finname, foutname, nameopt, nin, nout, true);
    First_Assembly_Reader firstreader(ctx, nin);
    oneapi::tbb::filter<void,MergeSplit_Context::TokInfo> queue = oneapi::tbb::make_filter<void, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, firstreader);
    if (nin > 1) {
        std::vector<Assembly_Reader> readers(nin-1);
        for (unsigned int i = 0; i < nin-1; ++i) {
            readers[i] = Assembly_Reader(ctx, i+1, nin);
        }
    for (auto & r : readers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, r);
    }
    std::vector<Assembler> assemblers(nin + nout + 1);
    for (unsigned int i = 0; i < nin + nout + 1; ++i) {
        assemblers[i] = Assembler(ctx, i);
    }
    for (auto &a : assemblers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::parallel, a);
    }
    }
    std::vector<Assembled_Writer> writers(nout-1);
    for (unsigned int i = 0; i < nout-1; ++i) {
        writers[i] = Assembled_Writer(i, nin, ctx.fouts[i]);
    }
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    Last_Assembled_Writer lastwriter(nin, ctx.fouts[nout-1]);
    oneapi::tbb::filter<void,void> full_queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(1, full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}

int assembleinbuffer_OH(const char* finname, const char* foutname, uint nameopt, uint nin, uint nout) {
    assert(nin > 0 && nout > 0);
    auto t0 = std::chrono::steady_clock::now();
    MergeSplit_Context ctx(finname, foutname, nameopt, nin, nout, true);
    First_Assembly_Reader_OH firstreader(ctx, nin);
    oneapi::tbb::filter<void,MergeSplit_Context::TokInfo> queue = oneapi::tbb::make_filter<void, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, firstreader);
    if (nin > 1) {
        std::vector<Assembly_Reader_OH> readers(nin-1);
        for (unsigned int i = 0; i < nin-1; ++i) {
            readers[i] = Assembly_Reader_OH(ctx, i+1, nin);
        }
    for (auto & r : readers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, r);
    }
    std::vector<Assembler_OH> assemblers(nin + nout + 1);
    for (unsigned int i = 0; i < nin + nout + 1; ++i) {
        assemblers[i] = Assembler_OH(ctx, i);
    }
    for (auto &a : assemblers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::parallel, a);
    }
    }
    std::vector<Assembled_Writer> writers(nout-1);
    for (unsigned int i = 0; i < nout-1; ++i) {
        writers[i] = Assembled_Writer(i, nin, ctx.fouts[i]);
    }
    for (auto & w : writers) {
        queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, MergeSplit_Context::TokInfo>(oneapi::tbb::filter_mode::serial_in_order, w);
    }
    Last_Assembled_Writer lastwriter(nin, ctx.fouts[nout-1]);
    oneapi::tbb::filter<void,void> full_queue = queue & oneapi::tbb::make_filter<MergeSplit_Context::TokInfo, void>(oneapi::tbb::filter_mode::serial_in_order, lastwriter);
    auto t1 = std::chrono::steady_clock::now();    
    oneapi::tbb::parallel_pipeline(1, full_queue);
    auto t2 = std::chrono::steady_clock::now();
    report(ctx.totorb, ctx.totev, ctx.totpuppi, t0, t1, t2);
    return 0;
}

/*ARGUMENTS:
1: method (mergesplit)
2: in file name (/scratch/gpetrucc/puppi -- will fill in .bx_n_of_m.dump)
3: out file name (/run/user/$UID/out_)
4: in file name opt (2 gets the bx... of...  format)
5: n in (6)
6: n out (6)
*/
int main(int argc, char**argv) {
    if (argc <= 2) {
        printf("Usage: %s test args\n", argv[0]);
        return 1;
    }
    std::string method(argv[1]);
    if (method == "mergesplit") {
        if (argc <= 6) {
            printf("ERROR: too few arguments; 7 expected");
            return 1;
        }
        return mergesplit(argv[2], argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]));
    } else if (method == "mergesplit_OH") {
        if (argc <= 6) {
            printf("ERROR: too few arguments; 7 expected");
            return 1;
        }
        return mergesplit_OH(argv[2], argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]));
    } else if (method == "assembleinbuffer") {
        if (argc <= 6) {
            printf("ERROR: too few arguments; 7 expected");
            return 1;
        }
        return assembleinbuffer(argv[2], argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]));
    } else if (method == "assembleinbuffer_OH") {
        if (argc <= 6) {
            printf("ERROR: too few arguments; 7 expected");
            return 1;
        }
        return assembleinbuffer_OH(argv[2], argv[3], std::atoi(argv[4]), std::atoi(argv[5]), std::atoi(argv[6]));
    } else {
    printf("ERROR: method not recognized \n");
    return 1;
    }
}
