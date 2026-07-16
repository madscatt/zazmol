#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "sasatomBase.h"

/********* declarations   ******************/



namespace sasio
{
    struct DcdReadHandle
    {
        FILE *fp;
        int natoms;
        int nframes;
        int reverseEndian;
        int charmm;
        int istart;
        int nsavc;
        double delta;
        int namnf;

        DcdReadHandle();
    };

    /// @brief (Brief description)<br>
    ///
    /// (Detailed descriptions)<br>
    class Files: public sasmol::SasAtomBase
    {
        public:
            void read_xyz(const std::string &filename) ; ///< (Brief description)

            void write_dcd(const std::string &filename); ///< (Brief description)

            void read_pdb(const std::string &filename); ///< (Brief description)

            void write_pdb(const std::string &filename, int &frame); ///< (Brief description)

            void write_dcd_step(FILE *outfile, int &frame, int &step); ///< (Brief description)

            void read_dcd_step(FILE *dcd_infile, int &frame, int &natoms, int &reverseEndian) ; ///< (Brief description)

            void read_dcd_step(DcdReadHandle &dcd_handle, int &frame) ; ///< (Brief description)

            void read_dcd(std::string &dcd_input_filename); ///< (Brief description)

            void resize_array() ; ///< (Brief description)

            void _set_unique_attributes() ; ///< (Brief description)
    };

    /// @note Where is my home?
    FILE *open_write_dcd_file(const std::string &filename, int &natoms, int &nset); ///< (Brief description)
    int write_dcd_header(FILE *outfile, const std::string &filename, int &natoms, int &nset); ///< (Brief description)
    DcdReadHandle open_dcd_read_handle(std::string &filename); ///< (Brief description)
    int read_dcd_step(DcdReadHandle &dcd_handle, int &frame, float *x, float *y, float *z); ///< (Brief description)
    int close_dcd_read(DcdReadHandle &dcd_handle); ///< (Brief description)
    FILE *open_dcd_read(std::string &filename, int &natoms, int &nset, int &reverseEndian) ; ///< (Brief description)
}
