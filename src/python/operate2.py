



class Move2():

    def initialize_subsets(self, other, self_basis, other_basis, frame=0):
        """
        Initialize the subsets for self and other based on the given basis.

        Parameters
        ----------
        other : Molecule
            The reference molecule.
        self_basis : str
            The basis for the self molecule.
        other_basis : str
            The basis for the other molecule.
        frame : int, optional
            The trajectory frame number to use (default is 0).

        Returns
        -------
        tuple
            A tuple containing the subsets for self and other.
        """
        # Initialize subset for other (molecule_1)
        error, other_mask = other.get_subset_mask(other_basis)
        subset_other = sasmol.system.Molecule()
        error = other.copy_molecule_using_mask(subset_other, other_mask, frame)
        subset_other.center(frame)

        # Initialize subset for self (molecule_2)
        error, self_mask = self.get_subset_mask(self_basis)
        subset_self = sasmol.system.Molecule()
        error = self.copy_molecule_using_mask(subset_self, self_mask, frame)
        subset_self.center(frame)

        return subset_self, subset_other


    def align(self, other, self_basis, other_basis, subset_self=None, subset_other=None, frame=0, **kwargs):
        """
        Alignment of one object on top of another.

        Parameters
        ----------
        other : Molecule
            The reference molecule.
        self_basis : str
            The basis for the self molecule.
        other_basis : str
            The basis for the other molecule.
        subset_self : Molecule, optional
            Precomputed subset for self (default is None).
        subset_other : Molecule, optional
            Precomputed subset for other (default is None).
        frame : int, optional
            The trajectory frame number to use (default is 0).

        Returns
        -------
        None
            Updated self._coor.
        """
        if subset_self is None or subset_other is None:
            subset_self, subset_other = self.initialize_subsets(other, self_basis, other_basis, frame)

        com_subset_other = subset_other.calculate_center_of_mass(frame)
        coor_subset_other = subset_other.coor()[frame]

        com_subset_self = subset_self.calculate_center_of_mass(frame)


        u = linear_algebra.find_u(coor_subset_self, coor_subset_other)

        tao = numpy.transpose(self.coor()[frame] - com_subset_other)

        error, nat2 = linear_algebra.matrix_multiply(u, tao)

        ncoor = numpy.transpose(nat2) + com_subset_other

        self._coor[frame, :] = ncoor

        return

