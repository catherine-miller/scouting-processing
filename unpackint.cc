/* ARGUMENTS:
(for multithreading: "-j #threads")
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
#include <filesystem>
//#include <ap_fixed.h>

using std::ifstream;

void readheader(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, UInt_t &npuppi) {
        fin.read(reinterpret_cast<char*>(&header), sizeof(uint64_t));
        npuppi = header          & 0xFFF;
        bx     = (header >> 12)  & 0xFFF;
        orbit  = (header >> 24)  & 0X3FFFFFFF;
        run    = (header >> 54);
}
void readcharged(uint64_t data, int16_t & z0, int8_t &dxy, UInt_t &quality, uint16_t &wpuppi, UInt_t &id) {
    z0 = ((data >> 49) & 1) ? ((data >> 40) | (- 0x200)) : ((data >> 40) & 0x3FF); //10 bits
    dxy = ((data >> 57) & 1) ? ((data >> 50) | (- 0x100)) : ((data >> 50) & 0xFF); //8 bits
    quality = (data >> 58) & 0x7; //3 bits
    //set all the neutral properties to 0
    wpuppi = 0;
    id = 0;
}
void readneutral(uint64_t data, int16_t & z0, int8_t &dxy, UInt_t &quality, uint16_t &wpuppi, UInt_t &id) {
    wpuppi = (data >> 23) & 0x3FF; //10 bits
    id = (data >> 13) & 0x3F; //6 bits
    //set all charged properties to 0
    z0 = 0; dxy = 0; quality = 0;
}
void readevent(std::fstream &fin, uint64_t &header, uint16_t &run, uint16_t &bx, uint32_t &orbit, 
                UInt_t &npuppi, uint64_t (&data)[255],
                uint16_t (&pt)[255], int16_t (&eta)[255], int16_t (&phi)[255], UInt_t (&pid)[255],
                int16_t (&z0)[255], int8_t (&dxy)[255], UInt_t (&quality)[255],
                uint16_t (&wpuppi)[255], UInt_t (&id)[255]
                ) {
    readheader(fin, header, run, bx, orbit, npuppi);
    if (npuppi) {
            fin.read(reinterpret_cast<char*>(&data[0]), npuppi*sizeof(uint64_t));
        }
    for (uint i = 0; i < npuppi; ++i) {
        pt[i] = data[i] & 0x3FFF;
        eta[i] = ((data[i] >> 25) & 1) ? ((data[i] >> 14) | (- 0x800)) : ((data[i] >> 14) & (0xFFF)); 
        phi[i] = ((data[i] >> 36) & 1) ? ((data[i] >> 26) | (- 0x400)) : ((data[i] >> 26) & (0x7FF));
        pid[i] = ((data[i]) >> 37) & 0x7;
        if (pid[i] > 1) readcharged(data[i], z0[i], dxy[i], quality[i], wpuppi[i], id[i]);
        else readneutral(data[i], z0[i], dxy[i], quality[i], wpuppi[i], id[i]);
    }
}

int main(int argc, char**argv) {
    uint64_t header, data[255];
    uint16_t run, bx;
    uint32_t orbit;
    UInt_t npuppi; // issues with uint8_t that root sees to max at 127
    //puppi candidate info:
    uint16_t pt[255];
    int16_t eta[255], phi[255];
    UInt_t pid[255];
    //charged only:
    int16_t z0[255];
    int8_t dxy[255];
    UInt_t quality[255];
    //neutral only:
    uint16_t wpuppi[255];
    UInt_t id[255];

    int iarg = 1, narg = argc - 1;
    if (std::string(argv[iarg]) == "-j") {
        ROOT::EnableImplicitMT(std::stoi(argv[iarg+1])); 
        printf("Enabled Implicit MT with %d threads\n", std::stoi(argv[iarg+1]));
        iarg += 2; narg -= 2;
    }
    std::fstream fin(argv[iarg], std::ios_base::in | std::ios_base::binary);

    int compressionAlgo, compressionLevel = 0;
    if (narg >= 4) {
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
    TTree *tree = new TTree("Events","Events");
    tree->Branch("run", &run, "run/s");
    tree->Branch("orbit", &orbit, "orbit/i");
    tree->Branch("bx", &bx, "bx/s");
    tree->Branch("nPuppi", &npuppi, "nPuppi/i");
    tree->Branch("pt", &pt, "pt[nPuppi]/s");
    tree->Branch("eta", &eta, "eta[nPuppi]/S");
    tree->Branch("phi", &phi, "phi[nPuppi]/S");
    tree->Branch("pid", &pid, "pid[nPuppi]/b");
    tree->Branch("z0", &z0, "z0[nPuppi]/S");
    tree->Branch("dxy", &dxy, "dxy[nPuppi]/B");
    tree->Branch("quality", &quality, "quality[nPuppi]/b");
    tree->Branch("wpuppi", &wpuppi, "wpuppi[nPuppi]/s");
    tree->Branch("id", &id, "id[nPuppi]/b");

    TStopwatch timer; timer.Start();
    unsigned long entries = 0;
    bool print = true;
    while(fin.good()) {
        readevent(fin, header, run, bx, orbit, npuppi, data, pt, eta, phi, pid, z0, dxy, quality, wpuppi, id);
        tree->Fill();
        /*if ((++entries) % 100000 == 0) {
            printf("Processed %8lu entries\n",entries);
        }*/
    }
    tree->Write();
    fout->Close();
    timer.Stop();
    double tcpu = timer.CpuTime(), treal = timer.RealTime();
    printf("Done in %.2fs (cpu), %.2fs (real). Event rate: %.1f kHz\n", tcpu, treal, entries/treal/1000.);
    //float insize = std::filesystem::file_size(argv[iarg]) / 1024. / 1024. ;
    //float outsize = std::filesystem::file_size(argv[iarg+1]) / 1024. / 1024. ;
    //printf("Input file size: %d MB, Output file size: %d MB\n", insize, outsize);
    return 0;
}

/* NOTES
compile with: gcc -lstdc++ unpackint.cc -o unpackint.exe $(root-config --cflags --ldflags --libs)
for example, run with (on lxplus8s10): 
    ./unpackint.exe /scratch/gpetrucc/fewPuppi.dump /eos/home-m/millerca/int_tree.root         (for a test)
    ./unpack.exe -j 100 /scratch/gpetrucc/puppi.bx_0_of_6.dump /run/user/$UID/out_tree.root
 */
