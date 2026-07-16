#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include <dcdio.h>
#include <sasio.h>

namespace {

int expect_true(bool condition, const std::string &message)
{
    if(!condition)
    {
        std::cerr << message << std::endl;
        return 1;
    }

    return 0;
}

int write_xplor_fixture(const std::string &filename)
{
    char writable_filename[256];
    std::snprintf(writable_filename, sizeof(writable_filename), "%s", filename.c_str());

    FILE *outfile = open_dcd_write(writable_filename);
    if(outfile == 0)
    {
        return 1;
    }

    int natoms = 2;
    int nframes = 1;
    int istart = 0;
    int nsavc = 1;
    double delta = 1.0;

    if(write_dcdheader(outfile, writable_filename, natoms, nframes, istart, nsavc, delta) != 1)
    {
        close_dcd_write(outfile);
        return 1;
    }

    float x[2] = {1.0f, 2.0f};
    float y[2] = {3.0f, 4.0f};
    float z[2] = {5.0f, 6.0f};

    if(write_dcdstep(outfile, natoms, x, y, z, 1) != 1)
    {
        close_dcd_write(outfile);
        return 1;
    }

    return close_dcd_write(outfile);
}

} // namespace

int main()
{
    int failures = 0;

    std::string charmm_fixture = "../../../src/python/test_sasmol/data/dcd_common/1ATM.dcd";
    sasio::DcdReadHandle charmm_handle = sasio::open_dcd_read_handle(charmm_fixture);

    failures += expect_true(charmm_handle.fp != 0, "CHARMM fixture did not open");
    if(charmm_handle.fp == 0)
    {
        return 1;
    }
    failures += expect_true(charmm_handle.natoms == 1, "CHARMM fixture natoms was not preserved");
    failures += expect_true(charmm_handle.nframes == 2, "CHARMM fixture nframes was not preserved");
    failures += expect_true(charmm_handle.reverseEndian == 0, "CHARMM fixture endian flag changed");
    failures += expect_true(charmm_handle.charmm == (DCD_IS_CHARMM | DCD_HAS_EXTRA_BLOCK),
                            "CHARMM fixture flags were not preserved");
    failures += expect_true(charmm_handle.istart == 0, "CHARMM fixture istart changed");
    failures += expect_true(charmm_handle.nsavc == 1, "CHARMM fixture nsavc changed");
    failures += expect_true(std::fabs(charmm_handle.delta - 1.0) < 1.0e-6, "CHARMM fixture delta changed");
    failures += expect_true(charmm_handle.namnf == 0, "CHARMM fixture namnf changed");

    std::vector<float> x(charmm_handle.natoms);
    std::vector<float> y(charmm_handle.natoms);
    std::vector<float> z(charmm_handle.natoms);
    int frame = 0;
    int read_result = sasio::read_dcd_step(charmm_handle, frame, &x[0], &y[0], &z[0]);
    failures += expect_true(read_result == 0, "CHARMM fixture frame read failed through descriptor state");
    failures += expect_true(std::fabs(x[0] - 76.944f) < 1.0e-5f, "CHARMM fixture x coordinate changed");
    failures += expect_true(std::fabs(y[0] - 41.799f) < 1.0e-5f, "CHARMM fixture y coordinate changed");
    failures += expect_true(std::fabs(z[0] - 41.652f) < 1.0e-5f, "CHARMM fixture z coordinate changed");
    failures += expect_true(sasio::close_dcd_read(charmm_handle) == 0, "CHARMM fixture close failed");
    failures += expect_true(charmm_handle.fp == 0, "CHARMM fixture close did not clear file pointer");

    std::string xplor_fixture = "test_dcd_read_handle_xplor.dcd";
    failures += expect_true(write_xplor_fixture(xplor_fixture) == 0, "X-PLOR fixture write failed");

    sasio::DcdReadHandle xplor_handle = sasio::open_dcd_read_handle(xplor_fixture);
    failures += expect_true(xplor_handle.fp != 0, "X-PLOR fixture did not open");
    failures += expect_true(xplor_handle.natoms == 2, "X-PLOR fixture natoms was not preserved");
    failures += expect_true(xplor_handle.nframes == 1, "X-PLOR fixture nframes was not preserved");
    failures += expect_true(xplor_handle.reverseEndian == 0, "X-PLOR fixture endian flag changed");
    failures += expect_true(xplor_handle.charmm == 0, "X-PLOR fixture should not be marked CHARMM");
    failures += expect_true(xplor_handle.istart == 0, "X-PLOR fixture istart changed");
    failures += expect_true(xplor_handle.nsavc == 1, "X-PLOR fixture nsavc changed");
    failures += expect_true(std::fabs(xplor_handle.delta - 1.0) < 1.0e-12, "X-PLOR fixture delta changed");
    failures += expect_true(xplor_handle.namnf == 0, "X-PLOR fixture namnf changed");
    failures += expect_true(sasio::close_dcd_read(xplor_handle) == 0, "X-PLOR fixture close failed");

    std::remove(xplor_fixture.c_str());

    return failures == 0 ? 0 : 1;
}
