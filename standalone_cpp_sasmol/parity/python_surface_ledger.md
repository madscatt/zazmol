# Python zazmol Surface Survey

Source: `/Users/curtisj/git_working_copies/zazmol/src/python`

| Module | Kind | Name | Line | C++ Status | Notes |
| --- | --- | --- | ---: | --- | --- |
| `calculate` | class | `Calculate` | 40 | question | |
| `calculate` | method | `Calculate.calculate_mass(self, **kwargs)` | 71 | question | |
| `calculate` | method | `Calculate.calculate_center_of_mass(self, frame, **kwargs)` | 127 | question | |
| `calculate` | method | `Calculate.calculate_radius_of_gyration(self, frame, **kwargs)` | 174 | question | |
| `calculate` | method | `Calculate.calculate_root_mean_square_deviation(self, other, **kwargs)` | 210 | question | |
| `calculate` | method | `Calculate.calculate_principal_moments_of_inertia(self, frame, **kwargs)` | 256 | question | |
| `calculate` | method | `Calculate.calculate_minimum_and_maximum(self, **kwargs)` | 323 | question | |
| `calculate` | method | `Calculate.calculate_minimum_and_maximum_all_steps(self, trajectory_filename=..., **kwargs)` | 400 | question | |
| `calculate` | method | `Calculate.calc_minmax_all_steps(self, dcdfilename, **kwargs)` | 475 | question | |
| `calculate` | method | `Calculate.calculate_residue_charge(self, **kwargs)` | 496 | question | |
| `calculate` | method | `Calculate.calculate_molecular_formula(self, **kwargs)` | 553 | question | |
| `charmm_topology` | class | `CharmmTopology` | 33 | question | |
| `charmm_topology` | method | `CharmmTopology.charmm_names(self, **kwargs)` | 52 | question | |
| `charmm_topology` | method | `CharmmTopology.add(self, dictionary, key, value, **kwargs)` | 104 | question | |
| `charmm_topology` | method | `CharmmTopology.read_charmm_topology(self, topology_file_path=..., topology_file_name=..., **kwargs)` | 155 | question | |
| `charmm_topology` | method | `CharmmTopology.patch_charmm_residue_atoms(self, residue, patch, **kwargs)` | 389 | question | |
| `charmm_topology` | method | `CharmmTopology.setup_cys_patch_atoms_simple(self, **kwargs)` | 468 | question | |
| `charmm_topology` | method | `CharmmTopology.setup_charmm_residue_atoms(self, **kwargs)` | 523 | question | |
| `charmm_topology` | method | `CharmmTopology.compare_list_ignore_order(self, l1, l2, **kwargs)` | 557 | question | |
| `charmm_topology` | method | `CharmmTopology.check_charmm_atomic_order_reorganize(self, **kwargs)` | 595 | question | |
| `dcd_io` | class | `DCD` | 61 | question | |
| `dcd_io` | method | `DCD.open_dcd_read(self, filename)` | 76 | question | |
| `dcd_io` | method | `DCD.open_dcd_write(self, filename)` | 111 | question | |
| `dcd_io` | method | `DCD.write_dcd_header(self, filepointer, nset)` | 141 | question | |
| `dcd_io` | method | `DCD.write_dcd_step(self, filepointer, frame, step)` | 156 | question | |
| `dcd_io` | method | `DCD.write_dcd_frames(self, filename, start, end)` | 172 | question | |
| `dcd_io` | method | `DCD.close_dcd_write(self, filepointer)` | 207 | question | |
| `dcd_io` | method | `DCD.close_dcd_read(self, dcd_file)` | 228 | question | |
| `dcd_io` | method | `DCD.write_dcd(self, filename)` | 242 | question | |
| `dcd_io` | method | `DCD.read_single_dcd_step(self, filename, frame)` | 288 | question | |
| `dcd_io` | method | `DCD.read_dcd_step(self, dcdfile, frame, **kwargs)` | 332 | question | |
| `dcd_io` | method | `DCD.read_dcd(self, filename)` | 388 | question | |
| `file_io` | class | `Files` | 47 | question | |
| `file_io` | method | `Files.open_file(filename)` | 76 | question | |
| `linear_algebra` | function | `cross_product(a, b)` | 35 | question | |
| `linear_algebra` | function | `matrix_multiply(a, b)` | 70 | question | |
| `linear_algebra` | function | `find_u(x, y)` | 129 | question | |
| `linear_algebra` | function | `vec_sub(a, b, c)` | 218 | question | |
| `linear_algebra` | function | `vec_scale(a, b, c)` | 245 | question | |
| `linear_algebra` | function | `cmp(a, b)` | 272 | question | |
| `linear_algebra` | function | `signed_angle(a, b, c)` | 279 | question | |
| `linear_algebra` | function | `dihedral_angle(a1, a2, a3, a4)` | 333 | question | |
| `linear_algebra` | function | `calculate_angle(a, b, c)` | 388 | question | |
| `logging_utilites` | class | `run_utils` | 17 | question | |
| `logging_utilites` | method | `run_utils.write_json_file(self)` | 30 | question | |
| `logging_utilites` | method | `run_utils.general_setup(self, other_self)` | 65 | question | |
| `logging_utilites` | method | `run_utils.preamble(self)` | 107 | question | |
| `logging_utilites` | method | `run_utils.setup_logging(self, other_self)` | 121 | question | |
| `logging_utilites` | method | `run_utils.print_gui(self, message)` | 161 | question | |
| `logging_utilites` | method | `run_utils.clean_up(self, log)` | 171 | question | |
| `logging_utilites` | method | `run_utils.capture_exception(self, message)` | 181 | question | |
| `multiprocessing_sasmol` | class | `Multiprocessing_SasMol` | 41 | question | |
| `multiprocessing_sasmol` | method | `Multiprocessing_SasMol.get_frame_lists(self, number_of_frames, number_of_batches, **kwargs)` | 80 | question | |
| `multiprocessing_sasmol` | method | `Multiprocessing_SasMol.example_worker(self, i, molecule, **kwargs)` | 122 | question | |
| `multiprocessing_sasmol` | method | `Multiprocessing_SasMol.divide_molecule(self, molecule, number_of_batches, **kwargs)` | 149 | question | |
| `multiprocessing_sasmol` | method | `Multiprocessing_SasMol.submit_jobs(self, molecules, target_method, number_of_jobs, **kwargs)` | 215 | question | |
| `operate` | class | `Move` | 45 | question | |
| `operate` | method | `Move.mass_check(self, **kwargs)` | 86 | question | |
| `operate` | method | `Move.translate(self, frame, value, **kwargs)` | 117 | question | |
| `operate` | method | `Move.center(self, frame, **kwargs)` | 188 | question | |
| `operate` | method | `Move.align(self, other, self_basis=..., other_basis=..., mode=..., align_variables=..., **kwargs)` | 240 | question | |
| `operate` | method | `Move.rotate(self, frame, axis, theta, **kwargs)` | 370 | question | |
| `operate` | method | `Move.rotate_general_axis(self, frame, theta, unit_axis, **kwargs)` | 427 | question | |
| `operate` | method | `Move.rotate_euler(self, frame, phi, theta, psi, **kwargs)` | 490 | question | |
| `operate` | method | `Move.align_pmi_on_cardinal_axes(self, frame)` | 557 | question | |
| `operate` | method | `Move.align_pmi_on_axis(self, frame, pmi_eigenvector, alignment_vector_axis)` | 567 | question | |
| `operate` | function | `set_average_vdw(mol)` | 607 | question | |
| `pdb_io` | class | `PDBElementResolutionError` | 55 | question | |
| `pdb_io` | class | `PDB` | 58 | question | |
| `pdb_io` | method | `PDB.print_error(self, name, my_message)` | 78 | question | |
| `pdb_io` | method | `PDB.check_error(self, error)` | 105 | question | |
| `pdb_io` | method | `PDB.element_filter(self)` | 120 | question | |
| `pdb_io` | method | `PDB.get_element(self, name, resname)` | 144 | question | |
| `pdb_io` | method | `PDB.write_pdb(self, filename, frame, flag, **kwargs)` | 200 | question | |
| `pdb_io` | method | `PDB.get_resnames(self)` | 327 | question | |
| `pdb_io` | method | `PDB.initialize_children(self)` | 344 | question | |
| `pdb_io` | method | `PDB.check_for_all_zero_columns(self, coor, frame=...)` | 415 | question | |
| `pdb_io` | method | `PDB.read_pdb(self, filename, **kwargs)` | 437 | question | |
| `pdb_io` | method | `PDB.create_conect_pdb_lines(self)` | 816 | question | |
| `properties` | class | `Atomic` | 36 | question | |
| `properties` | method | `Atomic.amu(self, **kwargs)` | 93 | question | |
| `properties` | method | `Atomic.amino_acid_sld(self, **kwargs)` | 152 | question | |
| `properties` | method | `Atomic.element_scattering_lengths(self, **kwargs)` | 209 | question | |
| `properties` | method | `Atomic.element_sl(self, **kwargs)` | 248 | question | |
| `properties` | method | `Atomic.nucleotide_scattering_lengths(self, **kwargs)` | 251 | question | |
| `properties` | method | `Atomic.nucleotide_sl(self, **kwargs)` | 279 | question | |
| `properties` | method | `Atomic.dna_scattering_lengths(self, **kwargs)` | 282 | question | |
| `properties` | method | `Atomic.dna_sl(self, **kwargs)` | 314 | question | |
| `properties` | method | `Atomic.rna_scattering_lengths(self, **kwargs)` | 317 | question | |
| `properties` | method | `Atomic.rna_sl(self, **kwargs)` | 345 | question | |
| `properties` | method | `Atomic.protein_scattering_lengths(self, **kwargs)` | 348 | question | |
| `properties` | method | `Atomic.protein_sl(self, **kwargs)` | 411 | question | |
| `properties` | method | `Atomic.van_der_waals_radii(self, **kwargs)` | 414 | question | |
| `properties` | method | `Atomic.van_der_Waals_radii(self, **kwargs)` | 465 | question | |
| `subset` | class | `Mask` | 52 | question | |
| `subset` | method | `Mask.get_dihedral_subset_mask(self, flexible_residues, mtype)` | 70 | question | |
| `subset` | method | `Mask.init_child(self, descriptor)` | 98 | question | |
| `subset` | method | `Mask.get_subset_mask(self, basis_filter)` | 194 | question | |
| `subset` | method | `Mask.merge_two_molecules(self, mol1, mol2, **kwargs)` | 322 | question | |
| `subset` | method | `Mask.copy_molecule_using_mask(self, other, mask, frame)` | 515 | question | |
| `subset` | method | `Mask.duplicate_molecule(self, number_of_duplicates, **kwargs)` | 695 | question | |
| `subset` | method | `Mask.get_indices_from_mask(self, mask)` | 744 | question | |
| `subset` | method | `Mask.get_coor_using_mask(self, frame, mask)` | 776 | question | |
| `subset` | method | `Mask.set_coor_using_mask(self, other, frame, mask)` | 828 | question | |
| `subset` | method | `Mask.set_descriptor_using_mask(self, mask, descriptor, value)` | 899 | question | |
| `subset` | method | `Mask.apply_biomt(self, frame, selection, U, M, **kwargs)` | 943 | question | |
| `subset` | method | `Mask.copy_apply_biomt(self, other, frame, selection, U, M, **kwargs)` | 999 | question | |
| `system` | class | `Atom` | 57 | question | |
| `system` | method | `Atom.setFilename(self, newValue)` | 165 | question | |
| `system` | method | `Atom.filename(self)` | 168 | question | |
| `system` | method | `Atom.setId(self, newValue)` | 225 | question | |
| `system` | method | `Atom.id(self)` | 228 | question | |
| `system` | method | `Atom.number_of_frames(self)` | 233 | question | |
| `system` | method | `Atom.setNumber_of_frames(self, newValue)` | 236 | question | |
| `system` | method | `Atom.atom(self)` | 239 | question | |
| `system` | method | `Atom.setAtom(self, newValue)` | 242 | question | |
| `system` | method | `Atom.index(self)` | 245 | question | |
| `system` | method | `Atom.setIndex(self, newValue)` | 248 | question | |
| `system` | method | `Atom.original_index(self)` | 251 | question | |
| `system` | method | `Atom.setOriginal_index(self, newValue)` | 254 | question | |
| `system` | method | `Atom.original_resid(self)` | 257 | question | |
| `system` | method | `Atom.setOriginal_resid(self, newValue)` | 260 | question | |
| `system` | method | `Atom.name(self)` | 263 | question | |
| `system` | method | `Atom.setName(self, newValue)` | 266 | question | |
| `system` | method | `Atom.loc(self)` | 269 | question | |
| `system` | method | `Atom.setLoc(self, newValue)` | 272 | question | |
| `system` | method | `Atom.resname(self)` | 275 | question | |
| `system` | method | `Atom.setResname(self, newValue)` | 278 | question | |
| `system` | method | `Atom.chain(self)` | 281 | question | |
| `system` | method | `Atom.setChain(self, newValue)` | 284 | question | |
| `system` | method | `Atom.resid(self)` | 287 | question | |
| `system` | method | `Atom.setResid(self, newValue)` | 290 | question | |
| `system` | method | `Atom.coor(self)` | 293 | question | |
| `system` | method | `Atom.setCoor(self, newValue)` | 296 | question | |
| `system` | method | `Atom.rescode(self)` | 299 | question | |
| `system` | method | `Atom.setRescode(self, newValue)` | 302 | question | |
| `system` | method | `Atom.occupancy(self)` | 305 | question | |
| `system` | method | `Atom.setOccupancy(self, newValue)` | 308 | question | |
| `system` | method | `Atom.beta(self)` | 311 | question | |
| `system` | method | `Atom.setBeta(self, newValue)` | 314 | question | |
| `system` | method | `Atom.segname(self)` | 317 | question | |
| `system` | method | `Atom.setSegname(self, newValue)` | 320 | question | |
| `system` | method | `Atom.element(self)` | 323 | question | |
| `system` | method | `Atom.setElement(self, newValue)` | 326 | question | |
| `system` | method | `Atom.charge(self)` | 329 | question | |
| `system` | method | `Atom.setCharge(self, newValue)` | 332 | question | |
| `system` | method | `Atom.atom_charge(self)` | 335 | question | |
| `system` | method | `Atom.setAtom_charge(self, newValue)` | 338 | question | |
| `system` | method | `Atom.atom_vdw(self)` | 341 | question | |
| `system` | method | `Atom.setAtom_vdw(self, newValue)` | 344 | question | |
| `system` | method | `Atom.formula(self)` | 347 | question | |
| `system` | method | `Atom.setFormula(self, newValue)` | 350 | question | |
| `system` | method | `Atom.mass(self)` | 353 | question | |
| `system` | method | `Atom.setMass(self, newValue)` | 356 | question | |
| `system` | method | `Atom.total_mass(self)` | 359 | question | |
| `system` | method | `Atom.setTotal_mass(self, newValue)` | 364 | question | |
| `system` | method | `Atom.natoms(self)` | 367 | question | |
| `system` | method | `Atom.setNatoms(self, newValue)` | 370 | question | |
| `system` | method | `Atom.number_of_atoms(self)` | 373 | question | |
| `system` | method | `Atom.setNumber_of_atoms(self, newValue)` | 376 | question | |
| `system` | method | `Atom.moltype(self)` | 379 | question | |
| `system` | method | `Atom.setMoltype(self, newValue)` | 382 | question | |
| `system` | method | `Atom.check_integrity(self, fast_check=...)` | 386 | question | |
| `system` | method | `Atom.validate_integrity(self, fast_check=...)` | 390 | question | |
| `system` | method | `Atom.record(self)` | 438 | question | |
| `system` | method | `Atom.setRecord(self, newValue)` | 441 | question | |
| `system` | method | `Atom.altloc(self)` | 444 | question | |
| `system` | method | `Atom.setAltloc(self, newValue)` | 447 | question | |
| `system` | method | `Atom.icode(self)` | 450 | question | |
| `system` | method | `Atom.setIcode(self, newValue)` | 453 | question | |
| `system` | class | `Molecule` | 457 | question | |
| `system` | method | `Molecule.fasta(self)` | 489 | question | |
| `system` | method | `Molecule.setFasta(self, newValue)` | 492 | question | |
| `system` | method | `Molecule.setCharmm_type(self, newValue)` | 495 | question | |
| `system` | method | `Molecule.charmm_type(self)` | 498 | question | |
| `system` | method | `Molecule.header(self)` | 501 | question | |
| `system` | method | `Molecule.setHeader(self, newValue)` | 504 | question | |
| `system` | method | `Molecule.unitcell(self)` | 507 | question | |
| `system` | method | `Molecule.setUnitcell(self, newValue)` | 510 | question | |
| `system` | method | `Molecule.rg(self)` | 513 | question | |
| `system` | method | `Molecule.setRg(self, newValue)` | 516 | question | |
| `system` | method | `Molecule.pmi(self)` | 519 | question | |
| `system` | method | `Molecule.setPmi(self, newValue)` | 522 | question | |
| `system` | method | `Molecule.minimum(self)` | 525 | question | |
| `system` | method | `Molecule.setMinimum(self, newValue)` | 528 | question | |
| `system` | method | `Molecule.maximum(self)` | 531 | question | |
| `system` | method | `Molecule.setMaximum(self, newValue)` | 534 | question | |
| `system` | method | `Molecule.one_letter_resname(self)` | 537 | question | |
| `system` | method | `Molecule.setOne_letter_resname(self, newValue)` | 540 | question | |
| `system` | method | `Molecule.center_of_mass(self)` | 543 | question | |
| `system` | method | `Molecule.setCenter_of_mass(self, newValue)` | 546 | question | |
| `system` | method | `Molecule.names(self)` | 549 | question | |
| `system` | method | `Molecule.setNames(self, newValue)` | 552 | question | |
| `system` | method | `Molecule.resnames(self)` | 555 | question | |
| `system` | method | `Molecule.setResnames(self, newValue)` | 558 | question | |
| `system` | method | `Molecule.resids(self)` | 561 | question | |
| `system` | method | `Molecule.setResids(self, newValue)` | 564 | question | |
| `system` | method | `Molecule.chains(self)` | 567 | question | |
| `system` | method | `Molecule.setChains(self, newValue)` | 570 | question | |
| `system` | method | `Molecule.segnames(self)` | 573 | question | |
| `system` | method | `Molecule.setSegnames(self, newValue)` | 576 | question | |
| `system` | method | `Molecule.occupancies(self)` | 579 | question | |
| `system` | method | `Molecule.setOccupancies(self, newValue)` | 582 | question | |
| `system` | method | `Molecule.betas(self)` | 585 | question | |
| `system` | method | `Molecule.setBetas(self, newValue)` | 588 | question | |
| `system` | method | `Molecule.moltypes(self)` | 591 | question | |
| `system` | method | `Molecule.setMoltypes(self, newValue)` | 594 | question | |
| `system` | method | `Molecule.elements(self)` | 597 | question | |
| `system` | method | `Molecule.setElements(self, newValue)` | 600 | question | |
| `system` | method | `Molecule.names_mask(self)` | 603 | question | |
| `system` | method | `Molecule.setNames_mask(self, newValue)` | 606 | question | |
| `system` | method | `Molecule.resnames_mask(self)` | 609 | question | |
| `system` | method | `Molecule.setResnames_mask(self, newValue)` | 612 | question | |
| `system` | method | `Molecule.resids_mask(self)` | 615 | question | |
| `system` | method | `Molecule.setResids_mask(self, newValue)` | 618 | question | |
| `system` | method | `Molecule.chains_mask(self)` | 621 | question | |
| `system` | method | `Molecule.setChains_mask(self, newValue)` | 624 | question | |
| `system` | method | `Molecule.occupancies_mask(self)` | 627 | question | |
| `system` | method | `Molecule.setOccupancies_mask(self, newValue)` | 630 | question | |
| `system` | method | `Molecule.betas_mask(self)` | 633 | question | |
| `system` | method | `Molecule.setBetas_mask(self, newValue)` | 636 | question | |
| `system` | method | `Molecule.elements_mask(self)` | 639 | question | |
| `system` | method | `Molecule.setElements_mask(self, newValue)` | 642 | question | |
| `system` | method | `Molecule.segnames_mask(self)` | 645 | question | |
| `system` | method | `Molecule.setSegnames_mask(self, newValue)` | 648 | question | |
| `system` | method | `Molecule.number_of_names(self)` | 651 | question | |
| `system` | method | `Molecule.setNumber_of_names(self, newValue)` | 654 | question | |
| `system` | method | `Molecule.number_of_resnames(self)` | 657 | question | |
| `system` | method | `Molecule.setNumber_of_resnames(self, newValue)` | 660 | question | |
| `system` | method | `Molecule.number_of_resids(self)` | 663 | question | |
| `system` | method | `Molecule.setNumber_of_resids(self, newValue)` | 666 | question | |
| `system` | method | `Molecule.number_of_chains(self)` | 669 | question | |
| `system` | method | `Molecule.setNumber_of_chains(self, newValue)` | 672 | question | |
| `system` | method | `Molecule.number_of_segnames(self)` | 675 | question | |
| `system` | method | `Molecule.setNumber_of_segnames(self, newValue)` | 678 | question | |
| `system` | method | `Molecule.number_of_occupancies(self)` | 681 | question | |
| `system` | method | `Molecule.setNumber_of_occupancies(self, newValue)` | 684 | question | |
| `system` | method | `Molecule.number_of_betas(self)` | 687 | question | |
| `system` | method | `Molecule.setNumber_of_betas(self, newValue)` | 690 | question | |
| `system` | method | `Molecule.number_of_moltypes(self)` | 693 | question | |
| `system` | method | `Molecule.setNumber_of_moltypes(self, newValue)` | 696 | question | |
| `system` | method | `Molecule.number_of_elements(self)` | 699 | question | |
| `system` | method | `Molecule.setNumber_of_elements(self, newValue)` | 702 | question | |
| `system` | method | `Molecule.residue_charge(self)` | 705 | question | |
| `system` | method | `Molecule.setResidue_charge(self, newValue)` | 708 | question | |
| `system` | method | `Molecule.conect(self)` | 711 | question | |
| `system` | method | `Molecule.setConect(self, newValue)` | 714 | question | |
| `system` | method | `Molecule.residue_flag(self)` | 717 | question | |
| `system` | method | `Molecule.setResidue_flag(self, newValue)` | 724 | question | |
| `system` | method | `Molecule.set_average_vdw(self)` | 727 | question | |
| `system` | class | `System` | 733 | question | |
| `system` | class | `Molecule_Maker` | 777 | question | |
| `system` | method | `Molecule_Maker.residue_flag(self)` | 914 | question | |
| `topology` | class | `Topology` | 33 | question | |
| `topology` | method | `Topology.renumber(self, **kwargs)` | 43 | question | |
| `topology` | method | `Topology.make_constraint_pdb(self, filename, basis_type, **kwargs)` | 150 | question | |
| `topology` | method | `Topology.make_backbone_pdb_from_fasta(self, filename, moltype, **kwargs)` | 285 | question | |
| `topology` | method | `Topology.create_fasta(self, **kwargs)` | 369 | question | |
| `utilities` | function | `find_unique(this_list)` | 40 | question | |
| `utilities` | class | `Copy_Using_Mask` | 57 | question | |
| `utilities` | method | `Copy_Using_Mask.from_sasmol(cls, class_instance, mask, **kwargs)` | 65 | question | |
| `utilities` | function | `duplicate_molecule(molecule, number_of_duplicates)` | 180 | question | |
| `utilities` | class | `Element` | 315 | question | |
| `utilities` | method | `Element.getweight(self)` | 326 | question | |
| `utilities` | method | `Element.addsyms(self, weight, result)` | 330 | question | |
| `utilities` | function | `build_dict(s)` | 345 | question | |
| `utilities` | class | `ElementSequence` | 367 | question | |
| `utilities` | method | `ElementSequence.append(self, thing)` | 376 | question | |
| `utilities` | method | `ElementSequence.getweight(self)` | 380 | question | |
| `utilities` | method | `ElementSequence.set_count(self, n)` | 387 | question | |
| `utilities` | method | `ElementSequence.addsyms(self, weight, result)` | 394 | question | |
| `utilities` | method | `ElementSequence.displaysyms(self, sym2elt)` | 400 | question | |
| `utilities` | class | `Tokenizer` | 410 | question | |
| `utilities` | method | `Tokenizer.gettoken(self)` | 419 | question | |
| `utilities` | method | `Tokenizer.error(self, msg)` | 442 | question | |
| `utilities` | function | `parse(s, sym2elt)` | 450 | question | |
| `utilities` | function | `parse_sequence(sym2elt)` | 476 | question | |
| `utilities` | function | `get_chemical_formula(formula_string)` | 517 | question | |
| `utilities` | function | `parse_fasta(fasta_sequence, **kwargs)` | 586 | question | |
| `utilities` | function | `check_integrity(obj, fast_check=..., warn=...)` | 654 | question | |
| `view` | class | `View` | 51 | question | |
| `view` | method | `View.send_coordinates_to_vmd(self, port, flag, **kwargs)` | 87 | question | |
