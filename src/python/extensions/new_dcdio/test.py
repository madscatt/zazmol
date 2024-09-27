#import sys
#sys.path.append('./')
import dcdio_module

print('dcdio_module.__file__:', dcdio_module.__file__)

filename = ('hiv1_gag_200_frames.dcd')

file_ptr = dcdio_module.open_dcd_read(filename)
if file_ptr:
    print(f"Successfully opened {filename}")
else:
    print(f"Failed to open {filename}")
#print(file_ptr.read())
#file_ptr.close()