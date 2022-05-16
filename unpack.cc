/* ARGUMENTS:
method name: "combined" (branches contain both charged and neutral) 
    or "separate" (two sets of branches, one for charged and one for neutral)
type: "float" or "int"
OPTIONAL: for multithreading: "-j #threads")
OPTIONAL: compression name
OPTIONAL: compression level
in file name
out file name
*/

#include <cstdio>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <Compression.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include <string>
//#include <filesystem>
#include <sys/stat.h>
#define PI 3.14159265359

void readheader(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, UInt_t &npuppi) {
        fin.read(reinterpret_cast<char*>(&header), sizeof(uint64_t));
        npuppi = header          & 0xFFF;
        bx     = (header >> 12)  & 0xFFF;
        orbit  = (header >> 24)  & 0X3FFFFFFF;
        run    = (header >> 54);
}
void assignpdgid(uint pid, short int &pdgid) {
    if (pid == 0) pdgid = 130;
    else if (pid == 1) pdgid = 22;
    else if (pid == 2) pdgid = -211;
    else if (pid == 3) pdgid = 211;
    else if (pid == 4) pdgid = 11;
    else if (pid == 5) pdgid = -11;
    else if (pid == 6) pdgid = 13;
    else if (pid == 7) pdgid = -13;
}
bool readpid(uint64_t &data, short int &pdgid) {
    uint pid = (data >> 37) & 0x7;
    assignpdgid(pid, pdgid);
    return (pid > 1);
}
bool readpid(uint64_t &data, short int &pdgid_c, short int &pdgid_n) { //overload for "charged/neutral" version
    uint pid = (data >> 37) & 0x7;
    if (pid > 1) {
        assignpdgid(pid, pdgid_c);
        return true;
    } else {
        assignpdgid(pid, pdgid_n);
        return false;
    }
}
void readshared(uint64_t &data, uint16_t &pt, int16_t &eta, int16_t &phi) { //int
    pt = data & 0x3FFF;
    eta = ((data >> 25) & 1) ? ((data >> 14) | (- 0x800)) : ((data >> 14) & (0xFFF));
    phi = ((data >> 36) & 1) ? ((data >> 26) | (- 0x400)) : ((data >> 26) & (0x7FF));
}
void readshared(uint64_t &data, float &pt, float &eta, float &phi) { //float
    uint ptint = data & 0x3FFF;
    pt = ptint * 0.25;

    int etaint = ((data >> 25) & 1) ? ((data >> 14) | (- 0x800)) : ((data >> 14) & (0xFFF));
    eta = etaint * PI / 720. ;

    int phiint = ((data >> 36) & 1) ? ((data >> 26) | (- 0x400)) : ((data >> 26) & (0x7FF));
    phi = phiint * PI / 720. ;
}
void readcharged(uint64_t data, int16_t &z0, int8_t &dxy, UInt_t &quality) { //int
    z0 = ((data >> 49) & 1) ? ((data >> 40) | (- 0x200)) : ((data >> 40) & 0x3FF);

    dxy = ((data >> 57) & 1) ? ((data >> 50) | (- 0x100)) : ((data >> 50) & 0xFF);
    quality = (data >> 58) & 0x7; //3 bits
}
void readcharged(uint64_t data, float & z0, int8_t &dxy, UInt_t &quality) { //float
    int z0int = ((data >> 49) & 1) ? ((data >> 40) | (- 0x200)) : ((data >> 40) & 0x3FF);
    z0 = z0int * .05; //conver to centimeters

    dxy = ((data >> 57) & 1) ? ((data >> 50) | (- 0x100)) : ((data >> 50) & 0xFF);
    quality = (data >> 58) & 0x7; //3 bits
}
void readneutral(uint64_t data, uint16_t &wpuppi, UInt_t &id) {
    wpuppi = (data >> 23) & 0x3FF;
    id = (data >> 13) & 0x3F;
}
void readneutral(uint64_t data, float &wpuppi, UInt_t &id) {
    int wpuppiint = (data >> 23) & 0x3FF;
    wpuppi = wpuppiint / 256;
    id = (data >> 13) & 0x3F;
}
void readevent(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, 
                UInt_t &npuppi, uint64_t (&data)[255],
                uint16_t (&pt)[255], int16_t (&eta)[255], int16_t (&phi)[255], UInt_t (&pid)[255],
                int16_t (&z0)[255], int8_t (&dxy)[255], UInt_t (&quality)[255],
                uint16_t (&wpuppi)[255], UInt_t (&id)[255]
                ) { //int, combined
    readheader(fin, header, run, bx, orbit, npuppi);
    if (npuppi) fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
    for (uint i = 0; i < npuppi; ++i) {
        readshared(data[i], pt[i], eta[i], phi[i]);
        pid[i] = (data[i] >> 37) & 0x7;
        if (pid[i] > 1) {
            readcharged(data[i], z0[i], dxy[i], quality[i]);
            wpuppi[i] = 0; id[i] = 0;
        } else {
            readneutral(data[i], wpuppi[i], id[i]);
            z0[i] = 0; dxy[i] = 0; quality[i] = 0;
        }
    }
}
void readevent(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, 
                UInt_t &npuppi, UInt_t &npuppi_c, UInt_t &npuppi_n, uint64_t (&data)[255],
                uint16_t (&pt_c)[255], uint16_t (&pt_n)[255], int16_t (&eta_c)[255], int16_t (&eta_n)[255],
                int16_t (&phi_c)[255], int16_t (&phi_n)[255], UInt_t (&pid_c)[255], UInt_t (&pid_n)[255],
                int16_t (&z0)[255], int8_t (&dxy)[255], UInt_t (&quality)[255],
                uint16_t (&wpuppi)[255], UInt_t (&id)[255]
                ) { //int, separate
    npuppi_c = 0; npuppi_n = 0;
    readheader(fin, header, run, bx, orbit, npuppi);
    if (npuppi) fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
    for (uint i = 0; i < npuppi; ++i) {
        UInt_t pid = (data[i] >> 37) & 0x7;
        if (pid > 1) {
            pid_c[i] = pid;
            readshared(data[i], pt_c[npuppi_c], eta_c[npuppi_c], phi_c[npuppi_c]);
            readcharged(data[i], z0[npuppi_c], dxy[npuppi_c], quality[npuppi_c]);
            npuppi_c++;
        } else {
            pid_n[i] = pid;
            readshared(data[i], pt_n[npuppi_n], eta_n[npuppi_n], phi_n[npuppi_n]);
            readneutral(data[i], wpuppi[npuppi_n], id[npuppi_n]);
            npuppi_n++;
        }
    }
}
void readevent(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, 
                UInt_t &npuppi, uint64_t (&data)[255],
                float (&pt)[255], float (&eta)[255], float (&phi)[255], short int (&pdgid)[255],
                float (&z0)[255], int8_t (&dxy)[255], UInt_t (&quality)[255],
                float (&wpuppi)[255], UInt_t (&id)[255]
                ) { //float, combined
    readheader(fin, header, run, bx, orbit, npuppi);
    if (npuppi) fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
    for (uint i = 0; i < npuppi; ++i) {
        readshared(data[i], pt[i], eta[i], phi[i]);
        if (readpid(data[i], pdgid[i])) {
            readcharged(data[i], z0[i], dxy[i], quality[i]);
            wpuppi[i] = 0; id[i] = 0;
        } else {
            readneutral(data[i], wpuppi[i], id[i]);
            z0[i] = 0; dxy[i] = 0; quality[i] = 0;
        }
    }
}
void readevent(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, 
                UInt_t &npuppi, UInt_t &npuppi_c, UInt_t &npuppi_n, uint64_t (&data)[255],
                float (&pt_c)[255], float (&pt_n)[255], float (&eta_c)[255], float (&eta_n)[255],
                float (&phi_c)[255], float (&phi_n)[255], short int (&pdgid_c)[255], short int (&pdgid_n)[255],
                float (&z0)[255], int8_t (&dxy)[255], UInt_t (&quality)[255],
                float (&wpuppi)[255], UInt_t (&id)[255]
                ) { //float, separate
    npuppi_c = 0; npuppi_n = 0;
    readheader(fin, header, run, bx, orbit, npuppi);
    if (npuppi) fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
    for (uint i = 0; i < npuppi; ++i) {
        if (readpid(data[i], pdgid_c[npuppi_c], pdgid_n[npuppi_n])) {
            readshared(data[i], pt_c[npuppi_c], eta_c[npuppi_c], phi_c[npuppi_c]);
            readcharged(data[i], z0[npuppi_c], dxy[npuppi_c], quality[npuppi_c]);
            npuppi_c++;
        } else {
            readshared(data[i], pt_n[npuppi_n], eta_n[npuppi_n], phi_n[npuppi_n]);
            readneutral(data[i], wpuppi[npuppi_n], id[npuppi_n]);
            npuppi_n++;
        }
    }
}

void report(double tcpu, double real, std::string infile, std::string outfile, int entries) {
    printf("Done in %.2fs (cpu), %.2fs (real). Event rate: %.1f kHz\n", tcpu, real, entries/real/1000.);
    struct stat stat_buf;
    int rc = stat(infile.c_str(), &stat_buf);
    float insize = ((rc == 0) ? stat_buf.st_size : -1) / 1024. / 1024.;
    rc = stat(outfile.c_str(), &stat_buf);
    float outsize = ((rc == 0) ? stat_buf.st_size : -1) / 1024. / 1024.;
    float ratio = outsize / insize;
    printf("Input file size: %f MB, Output file size: %f MB, File size ratio: %f \n", insize, outsize, ratio);
}

int main(int argc, char**argv) {
    std::string method = std::string(argv[1]);
    std::string type = std::string(argv[2]);
    int iarg = 3, narg = argc - 1;
    if (std::string(argv[iarg]) == "-j") {
        ROOT::EnableImplicitMT(std::stoi(argv[iarg+1])); 
        printf("Enabled Implicit MT with %d threads\n", std::stoi(argv[iarg+1]));
        iarg += 2; narg -= 2;
    }
    std::fstream fin(argv[iarg], std::ios_base::in | std::ios_base::binary);

    int compressionAlgo, compressionLevel = 0;
    if (narg >= 5) {
        std::string compressionName(argv[iarg+2]);
        if (compressionName == "lzma") compressionAlgo = ROOT::RCompressionSetting::EAlgorithm::kLZMA;
        else if (compressionName == "zlib") compressionAlgo = ROOT::RCompressionSetting::EAlgorithm::kZLIB;
        else if (compressionName == "lz4") compressionAlgo = ROOT::RCompressionSetting::EAlgorithm::kLZ4;
        else if (compressionName == "zstd") compressionAlgo = ROOT::RCompressionSetting::EAlgorithm::kZSTD;
        else {
            printf("Unsupported compression algo %s\n", argv[iarg+2]);
            return 1;
        }
        compressionLevel = std::stoi(argv[iarg+3]);
    }
    TFile *fout = TFile::Open(argv[iarg+1], "RECREATE", "", compressionLevel);
    if (compressionLevel) fout->SetCompressionAlgorithm(compressionAlgo);

    TTree *tree = new TTree("Events", "Events");
    if (type == "int" && method == "combined") {
        uint64_t header, data[255];
        uint16_t run, bx;
        uint32_t orbit;
        UInt_t npuppi; // issues with uint8_t that root sees to max at 127
        //puppi candidate info:
        uint16_t pt[255];
        int16_t eta[255], phi[255];
        UInt_t pdgid[255]; //this is really pid but we need to call it pdg
        //charged only:
        int16_t z0[255];
        int8_t dxy[255];
        UInt_t quality[255];
        //neutral only:
        uint16_t wpuppi[255];
        UInt_t id[255];

        tree->Branch("run", &run, "run/s");
        tree->Branch("orbit", &orbit, "orbit/i");
        tree->Branch("bx", &bx, "bx/s");
        tree->Branch("nPuppi", &npuppi, "nPuppi/i");
        tree->Branch("pt", &pt, "pt[nPuppi]/s");
        tree->Branch("eta", &eta, "eta[nPuppi]/S");
        tree->Branch("phi", &phi, "phi[nPuppi]/S");
        tree->Branch("pid", &pdgid, "pid[nPuppi]/b");
        tree->Branch("z0", &z0, "z0[nPuppi]/S");
        tree->Branch("dxy", &dxy, "dxy[nPuppi]/B");
        tree->Branch("quality", &quality, "quality[nPuppi]/b");
        tree->Branch("wpuppi", &wpuppi, "wpuppi[nPuppi]/s");
        tree->Branch("id", &id, "id[nPuppi]/b");

        TStopwatch timer; timer.Start();
        unsigned long entries = 0;
        bool print = true;
        while(fin.good()) {
            readevent(fin, header, run, bx, orbit, npuppi, data, pt, eta, phi, pdgid, z0, dxy, quality, wpuppi, id);
            tree->Fill();
            entries++;
        }
        tree->Write();
        fout->Close();
        timer.Stop();
        double tcpu = timer.CpuTime(), treal = timer.RealTime();
        report(tcpu, treal, std::string(argv[iarg]), std::string(argv[iarg + 1]), entries);
        return 0;
    } else if (type == "int" && method == "separate") {
        uint64_t header, data[255];
        uint16_t run, bx;
        uint32_t orbit;
        UInt_t npuppi, npuppi_c, npuppi_n; // issues with uint8_t that root sees to max at 127
        //puppi candidate info:
        uint16_t pt_c[255], pt_n[255];
        int16_t eta_c[255], eta_n[255], phi_c[255], phi_n[255];
        UInt_t pdgid_c[255], pdgid_n[255]; //this is really pid but we need to give it the same name
        //charged only:
        int16_t z0[255];
        int8_t dxy[255];
        UInt_t quality[255];
        //neutral only:
        uint16_t wpuppi[255];
        UInt_t id[255];

        tree->Branch("run", &run, "run/s");
        tree->Branch("orbit", &orbit, "orbit/i");
        tree->Branch("bx", &bx, "bx/s");
        tree->Branch("nPuppi", &npuppi, "nPuppi/i");
        //charged branches
        tree->Branch("nPuppi_c", &npuppi_c, "nPuppi_c/i");
        tree->Branch("pt_c", &pt_c, "pt_c[nPuppi_c]/s");
        tree->Branch("eta_c", &eta_c, "eta_c[nPuppi_c]/S");
        tree->Branch("phi_c", &phi_c, "phi_c[nPuppi_c]/S");
        tree->Branch("pid_c", &pdgid_c, "pid_c[nPuppi_c]/b");
        tree->Branch("z0", &z0, "z0[nPuppi_c]/S");
        tree->Branch("dxy", &dxy, "dxy[nPuppi_c]/B");
        tree->Branch("quality", &quality, "quality[nPuppi]/b");
        //neutral branches
        tree->Branch("nPuppi_n", &npuppi_n, "nPuppi_n/i");
        tree->Branch("pt_n", &pt_n, "pt_n[nPuppi_n]/s");
        tree->Branch("eta_n", &eta_n, "eta_n[nPuppi_n]/S");
        tree->Branch("phi_n", &phi_n, "phi_n[nPuppi_c]/S");
        tree->Branch("pid_n", &pdgid_n, "pid_n[nPuppi_n]/b");
        tree->Branch("wpuppi", &wpuppi, "wpuppi[nPuppi_n]/s");
        tree->Branch("id", &id, "id[nPuppi]/b");

        TStopwatch timer; timer.Start();
        unsigned long entries = 0;
        bool print = true;
        while(fin.good()) {
            readevent(fin, header, run, bx, orbit, npuppi, npuppi_c, npuppi_n,
            data, pt_c, pt_n, eta_c, eta_n, phi_c, phi_n, pdgid_c, pdgid_n,
            z0, dxy, quality, wpuppi, id);
            tree->Fill();
            entries++;
        }
        tree->Write();
        fout->Close();
        timer.Stop();
        double tcpu = timer.CpuTime(), treal = timer.RealTime();
        report(tcpu, treal, std::string(argv[iarg]), std::string(argv[iarg + 1]), entries);
        return 0;
    }
    else if (type == "float" && method == "combined") {
            uint64_t header, data[255];
            uint16_t run, bx;
            uint32_t orbit;
            UInt_t npuppi; // issues with uint8_t that root sees to max at 127
            //puppi candidate info:
            float pt[255];
            float eta[255], phi[255];
            short int pdgid[255];
            //charged only:
            float z0[255];
            int8_t dxy[255];
            UInt_t quality[255];
            //neutral only:
            float wpuppi[255];
            UInt_t id[255];

            tree->Branch("run", &run, "run/s");
            tree->Branch("orbit", &orbit, "orbit/i");
            tree->Branch("bx", &bx, "bx/s");
            tree->Branch("nPuppi", &npuppi, "nPuppi/i");
            tree->Branch("pt", &pt, "pt[nPuppi]/f");
            tree->Branch("eta", &eta, "eta[nPuppi]/f");
            tree->Branch("phi", &phi, "phi[nPuppi]/f");
            tree->Branch("pdgid", &pdgid, "pdgid[nPuppi]/b");
            tree->Branch("z0", &z0, "z0[nPuppi]/f");
            tree->Branch("dxy", &dxy, "dxy[nPuppi]/f");
            tree->Branch("wpuppi", &wpuppi, "wpuppi[nPuppi]/f");
            tree->Branch("id", &id, "id[nPuppi]/b");

            TStopwatch timer; timer.Start();
            unsigned long entries = 0;
            bool print = true;
            while(fin.good()) {
                readevent(fin, header, run, bx, orbit, npuppi, data, pt, eta, phi, pdgid, z0, dxy, quality, wpuppi, id);
                tree->Fill();
                entries++;
            }
            tree->Write();
            fout->Close();
            timer.Stop();
            double tcpu = timer.CpuTime(), treal = timer.RealTime();
            report(tcpu, treal, std::string(argv[iarg]), std::string(argv[iarg + 1]), entries);
            return 0;
    } else if (type == "float" && method == "separate") {
        uint64_t header, data[255];
        uint16_t run, bx;
        uint32_t orbit;
        UInt_t npuppi, npuppi_c, npuppi_n; // issues with uint8_t that root sees to max at 127
        //puppi candidate info:
        float pt_c[255], pt_n[255], eta_c[255], eta_n[255], phi_c[255], phi_n[255];
        short int pdgid_c[255], pdgid_n[255];
        //charged only:
        float z0[255];
        int8_t dxy[255];
        UInt_t quality[255];
        //neutral only:
        float wpuppi[255];
        UInt_t id[255];

        tree->Branch("run", &run, "run/s");
        tree->Branch("orbit", &orbit, "orbit/i");
        tree->Branch("bx", &bx, "bx/s");
        tree->Branch("nPuppi", &npuppi, "nPuppi/i");
        //charged branches
        tree->Branch("nPuppi_c", &npuppi_c, "nPuppi_c/i");
        tree->Branch("pt_c", &pt_c, "pt_c[nPuppi_c]/f");
        tree->Branch("eta_c", &eta_c, "eta_c[nPuppi_c]/f");
        tree->Branch("phi_c", &phi_c, "phi_c[nPuppi_c]/f");
        tree->Branch("pdgid_c", &pdgid_c, "pdgid_c[nPuppi_c]/b");
        tree->Branch("quality", &quality, "quality[nPuppi_c]/b");
        tree->Branch("z0", &z0, "z0[nPuppi_c]/f");
        tree->Branch("dxy", &dxy, "dxy[nPuppi_c]/f");
        //neutral branches
        tree->Branch("nPuppi_n", &npuppi_n, "nPuppi_n/i");
        tree->Branch("pt_n", &pt_n, "pt_n[nPuppi_n]/f");
        tree->Branch("eta_n", &eta_n, "eta_n[nPuppi_n]/f");
        tree->Branch("phi_n", &phi_n, "phi_n[nPuppi_n]/f");
        tree->Branch("pdgid_n", &pdgid_n, "pdgid_n[nPuppi_n]/b");;
        tree->Branch("wpuppi", &wpuppi, "wpuppi[nPuppi_n]/f");
        tree->Branch("id", &id, "id[nPuppi_n]/b");

        TStopwatch timer; timer.Start();
        unsigned long entries = 0;
        bool print = true;
        while(fin.good()) {
            readevent(fin, header, run, bx, orbit, npuppi, npuppi_c, npuppi_n,
            data, pt_c, pt_n, eta_c, eta_n, phi_c, phi_n, pdgid_c, pdgid_n,
            z0, dxy, quality, wpuppi, id);
            tree->Fill();
            entries++;
        }
        tree->Write();
        fout->Close();
        timer.Stop();
        double tcpu = timer.CpuTime(), treal = timer.RealTime();
        report(tcpu, treal, std::string(argv[iarg]), std::string(argv[iarg + 1]), entries);
        return 0;
    } else printf("ERROR: method or type not recognized. Method should be either combined or separate. Type should be either float or int. \n");
}

/* NOTES
compile with: gcc -lstdc++ unpack.cc -o unpack.exe $(root-config --cflags --ldflags --libs)
for example, run with (on lxplus8s10): 
    ./unpack.exe combined float /scratch/gpetrucc/fewPuppi.dump /eos/home-m/millerca/float_combined.root         (for a test)
    ./unpack.exe separate int /scratch/gpetrucc/fewPuppi.dump /eos/home-m/millerca/float_separate.root
    to test run speed:
        ./unpack.exe combined float -j 100 /scratch/gpetrucc/puppi.bx_0_of_6.dump /run/user/$UID/float_tree.root
 */
