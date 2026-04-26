import os

###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###
###     BEGIN SYSADMIN EDIT ###

arch = "cluster"
# arch = "mac"

__cuda__ = True
__cuda__ = False

global __logging__
__logging__ = False

# __logging__ controls Python-side debug verbosity only: DCD progress dots,
# logging_utilites debug/info chatter, and the default _debug flag on system
# objects. It does not change numerical/scientific behavior or C debug macros.


###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###
###     END SYSADMIN EDIT ###

global __sasmol_id__
__sasmol_id__ = -1


__logging_level__ = 'ERROR'

if __logging__:
    __logging_level__ = 'DEBUG'


if arch == "cluster":

    installation_bin_path = ['share', 'apps', 'local', 'bin']
    __core_libraries_include__ = [os.path.join(
        os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(
        os.path.sep, 'share', 'apps', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(
            os.path.sep, 'share', 'apps', 'local', 'cuda-6.0')
else:

    installation_bin_path = ['usr', 'local', 'bin']
    __core_libraries_include__ = [os.path.join(
        os.path.sep, 'usr', 'local', 'core_libraries', 'include')]
    __core_libraries_lib__ = [os.path.join(
        os.path.sep, 'usr', 'local', 'core_libraries', 'lib')]

    if __cuda__:
        __cuda_path__ = os.path.join(os.path.sep, 'usr', 'local', 'cuda')

__bin_path__ = os.path.join(os.path.sep, *installation_bin_path)
