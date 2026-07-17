"""Build metadata for the packaged modern native C++ SasMol library."""

import os


_NATIVE_API_VERSION = 1
_CAPABILITIES = frozenset(
    {
        "modern_dcd_reader",
        "dcd_read_header",
        "dcd_read_next_frame_coordinates",
    }
)


class NativeConfigError(RuntimeError):
    """Raised when the installed native SasMol payload is incomplete."""


def _native_root():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "native")


def _require_path(path, description):
    if not os.path.exists(path):
        raise NativeConfigError(
            "installed sasmol does not include %s at %s" % (description, path)
        )
    return path


def include_dirs():
    """Return include directories for downstream native extension builds."""

    include_dir = _require_path(
        os.path.join(_native_root(), "include"), "native headers")
    return [include_dir]


def library_dirs():
    """Return library directories for downstream build systems that need them."""

    library_dir = _require_path(
        os.path.join(_native_root(), "lib"), "native libraries")
    return [library_dir]


def extra_objects():
    """Return the exact static archive downstream extensions should link."""

    archive = _require_path(
        os.path.join(_native_root(), "lib", "libsasmol.a"), "libsasmol.a")
    return [archive]


def native_api_version():
    """Return the version of the packaged native build contract."""

    return _NATIVE_API_VERSION


def capabilities():
    """Return named native capabilities available in this installation."""

    return tuple(sorted(_CAPABILITIES))


def has_capability(name):
    """Return whether a named native capability is available."""

    return name in _CAPABILITIES


def has_dcd_read_handle_api():
    """Return False because the removed legacy sasio handle ABI is not shipped."""

    return False


def has_modern_dcd_reader_api():
    """Return whether the C++20 DcdReader streaming API is shipped."""

    return has_capability("modern_dcd_reader")
