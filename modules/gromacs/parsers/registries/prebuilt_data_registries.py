acpype_generated_data = {
    "molecules": ["Compound", "nmols"],  # itp, top
    "system": ["Compound"],  # top
    "moleculetype": ["name", "nrexcl"],  # itp
    "atomtypes": [
        "name",
        "bond_type",
        "mass",
        "charge",
        "ptype",
        "sigma",
        "epsilon",
    ],  # itp
    "atoms": ["nr", "type", "resi", "res", "atom", "cgnr", "charge", "mass"],  # itp
    "bonds": ["ai", "aj", "funct", "r", "k"],  # itp
    "pairs": ["ai", "aj", "funct"],
    "angles": ["ai", "aj", "ak", "funct", "theta", "cth"],
    "dihedrals": ["i", "j", "k", "l", "funct", "phase", "kd", "pn"],
}
