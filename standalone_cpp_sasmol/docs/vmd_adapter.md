# VMD Adapter

Python `sasmol.view.View.send_coordinates_to_vmd(port, flag)` extracts frame 0
coordinates as float32 X/Y/Z arrays and passes them to the compiled VMD IMD
socket helper. The standalone C++ port keeps that same split:

- `prepare_vmd_coordinate_arrays(molecule, frame)` builds SasMol-owned X/Y/Z
  `coord_type` arrays from one molecule frame
- `send_coordinate_arrays_to_vmd(x, y, z, port, clear_socket)` sends existing
  arrays to the active VMD transport
- `send_coordinates_to_vmd(molecule, port, clear_socket, frame)` is the
  Python-style convenience wrapper

The adapter is optional. Portable core builds leave
`SASMOL_ENABLE_VMD_ADAPTER` off and return `ViewCode::unsupported` from the
default sender. Builds that opt in with:

```bash
cmake -S standalone_cpp_sasmol -B standalone_cpp_sasmol/build-vmd -DSASMOL_ENABLE_VMD_ADAPTER=ON
```

compile the promoted legacy IMD/VMD C transport in `src/vmd_legacy`.

Tests do not require VMD to be open. They validate coordinate extraction,
array-shape rejection, transport-status propagation, and opt-out behavior using
a mock sender.
