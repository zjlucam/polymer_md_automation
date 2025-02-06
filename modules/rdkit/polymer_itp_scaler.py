import pandas as pd
import pandas as pd
<<<<<<< HEAD
import string
=======
>>>>>>> 91758eb (cleaned up)
import re


class PolymerITPScaler:
    num_cols = {
        "bonds": ["ai", "aj"],
        "pairs": ["ai", "aj"],
        "angles": ["ai", "aj", "ak"],
        "dihedrals": ["ai", "aj", "ak", "al"],
        "impropers": ["ai", "aj", "ak", "al"],
        "atoms": ["nr", "cgnr"],
    }
    SECTION_HEADERS = {
        "atoms": ["nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass"],
        "bonds": ["ai", "aj", "funct", "r", "k"],
        "pairs": ["ai", "aj", "funct"],
        "angles": ["ai", "aj", "ak", "funct", "theta", "k_theta"],
        "dihedrals": ["ai", "aj", "ak", "al", "funct", "phi", "k_phi", "pn"],
        "impropers": ["ai", "aj", "ak", "al", "funct", "phi", "k_phi", "pn"],
    }
    comment_sections = {
        k: v for k, v in num_cols.items() if k != "atoms"
    }  # Exclude "atoms"

    def __init__(self, itp_path, short_cg_map, n_repeat, atom_start_index=None):
        """
        Args:
            itp_path (str): Path to .itp file.
            cg_map (list): Mapping of residue indices.
            n_repeat (int): Number of times to repeat the middle section.
            atom_start_index (int or None): Determines atom name counting style.
                - None → "C", "H", etc.
                - 0 → "C0", "H0", etc.
                - 1 → "C1", "H1", etc.
        """
        self.itp_path = itp_path
        self.cg_map = short_cg_map
        self.n_repeat = n_repeat
        self.atom_start_index = atom_start_index  # Determines how atom indices start
        self.sections = {}

        self._middle_length = self._calculate_middle_length()
        self._parse_itp()
        self._process_sections()

    def _parse_itp(self) -> None:
        """Parses .itp file into DataFrames while preserving interaction types."""
        itp_file = self.itp_path
        sections = {}
        current_section = None
        section_data = []

        with open(itp_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                # **Remove in-line comments** before processing
                line = re.split(r"\s*;\s*", line)[0]

                # Detect section headers
                match = re.match(r"^\[\s*(\w+)\s*\](.*)$", line)
                if match:
                    if current_section and section_data:
                        sections[current_section] = self._create_dataframe(
                            current_section, section_data
                        )
                        section_data = []

                    section_name, extra_info = match.groups()
                    section_name = section_name.lower()

                    if "improper" in extra_info.lower():
                        current_section = "impropers"
                    else:
                        current_section = section_name

                    continue

                if current_section:
                    section_data.append(line.split())

        # Convert last section into DataFrame
        if current_section and section_data:
            sections[current_section] = self._create_dataframe(
                current_section, section_data
            )

        self.sections = sections

    def _create_dataframe(self, section_name, data):
        """Creates a DataFrame with predefined headers, handling extra/missing columns."""
        headers = self.SECTION_HEADERS.get(section_name, None)
        if headers:

            df = pd.DataFrame(data, columns=headers)

            # Fill missing headers with NaN
            for col in self.SECTION_HEADERS.get(section_name, []):
                if col not in df.columns:
                    df[col] = None
        else:
            df = pd.DataFrame(data)  # Generic DataFrame for unknown sections

        return df

    def _calculate_middle_length(self):
        """Determines the number of atoms in the repeating middle section."""
        middle_atoms = sum(len(entry["atom_indices"]) for entry in self.cg_map[1:-1])
        return middle_atoms

    def _shift_atom_indices(self):
        """Shifts atom names based on the defined counting style."""
        if "atoms" not in self.sections:
            return

        df = self.sections["atoms"]
        atom_counts = {}  # Track occurrences of each atom type

        for i, row in df.iterrows():
            atom_name = row["atom"]
            match = re.match(
                r"([A-Za-z]+)(\d*)", atom_name
            )  # Extract base name & number
            if match:
                base_name, number = match.groups()
                if base_name not in atom_counts:
<<<<<<< HEAD
                    atom_counts[base_name] =  self.atom_start_index or 0 
                    

                atom_index, atom_counts[base_name] = self._num_to_alphabet_gromacs_name(atom_counts[base_name])
                df.at[i, "atom"] = f"{base_name}{atom_index}"
=======
                    atom_counts[base_name] = self.atom_start_index or 0  # Start count

                df.at[i, "atom"] = f"{base_name}{atom_counts[base_name]}"
>>>>>>> 91758eb (cleaned up)
                atom_counts[base_name] += 1  # Increment count

        self.sections["atoms"] = df

<<<<<<< HEAD
    def _num_to_alphabet_gromacs_name(self, index):
        if index <100:
            return index, index
        else:
            base_index = (index-100)//10
            sub_index = (index-100)%10
            letter = string.ascii_uppercase[base_index]
            return f"{letter}{sub_index}", index
=======
>>>>>>> 91758eb (cleaned up)
    def _process_sections(self):
        """Processes all relevant sections, shifts indices, and duplicates the middle section."""
        num_cols = {
            "bonds": ["ai", "aj"],
            "pairs": ["ai", "aj"],
            "angles": ["ai", "aj", "ak"],
            "dihedrals": ["ai", "aj", "ak", "al"],
            "impropers": ["ai", "aj", "ak", "al"],
            "atoms": ["nr", "cgnr"],
        }

        for section, cols in num_cols.items():
            if section in self.sections:
                df = self.sections[section]

                # Ensure numeric conversion
                for col in cols:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

                # Identify start, middle, and end
                start_indices = set(self.cg_map[0]["atom_indices"])
                middle_indices = set(
                    idx for entry in self.cg_map[1:-1] for idx in entry["atom_indices"]
                )

                start_df = df[df[cols[0]].isin(start_indices)].copy()
                middle_df = df[df[cols[0]].isin(middle_indices)].copy()

                assigned_indices = start_indices.union(middle_indices)
                all_indices = set(df[cols[0]].unique())

                # End should contain everything that isn't in start or middle
                end_indices = all_indices - assigned_indices
                end_df = df[df[cols[0]].isin(end_indices)].copy()

                # Duplicate middle n times and shift indices accordingly
                all_middle = []
                for i in range(self.n_repeat + 1):
                    temp_df = middle_df.copy()
                    for col in cols:
                        temp_df[col] += (
                            i * self._middle_length
                        )  # Shift middle section indices
                    all_middle.append(temp_df)

                # 🔄 Apply index shift to `end_df`
                for col in cols:
                    end_df[col] += (
                        self.n_repeat * self._middle_length
                    )  # Shift end section indices

                # Combine everything back together
                df = pd.concat([start_df] + all_middle + [end_df], ignore_index=True)
                self.sections[section] = df

        self._shift_atom_indices()
        self._add_comments()

    def _add_comments(self):
        """Adds comments to relevant sections based on atom names."""
        if "atoms" not in self.sections:
            return

        atom_map = {
            row["nr"]: row["atom"] for _, row in self.sections["atoms"].iterrows()
        }  # nr -> atom_name

        for section, cols in self.comment_sections.items():
            if section in self.sections:
                df = self.sections[section]
                df["comment"] = df.apply(
                    lambda row: ";"
                    + " - ".join(atom_map.get(row[col], "?") for col in cols),
                    axis=1,
                )

    def write_itp(self, output_path):
        """Writes the modified ITP file to the specified output path."""
        with open(output_path, "w") as f:
            for section, df in self.sections.items():
                f.write(f"\n[ {section} ]\n")
                for _, row in df.iterrows():
                    line = " ".join(map(str, row.dropna().tolist()))
                    f.write(f"{line}\n")
