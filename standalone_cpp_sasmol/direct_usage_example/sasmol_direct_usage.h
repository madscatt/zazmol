#ifndef SASMOL_DIRECT_USAGE_H
#define SASMOL_DIRECT_USAGE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int sasmol_pdb_to_dcd(const char* pdb_path, const char* dcd_path,
                      char* error_message, size_t error_message_size,
                      int* atom_count, int* frame_count);

#ifdef __cplusplus
}
#endif

#endif
