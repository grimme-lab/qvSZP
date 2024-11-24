"""
This script reads a CSV file containing the paths to directories of compounds
and extracts the atomic charges from the
files in the directories. The atomic charges are then written to a new CSV file.
"""

from abc import ABC, abstractmethod
from pathlib import Path
import json
import warnings
import argparse

import pandas as pd
from tqdm import tqdm

PSE_SYMBOLS = {
    0: "X",
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Nh",
    114: "Fl",
    115: "Mc",
    116: "Lv",
    117: "Ts",
    118: "Og",
}
PSE_NUMBERS: dict[str, int] = {k.lower(): v for v, k in PSE_SYMBOLS.items()}


# Define the abstract base class for charge parsers
class ChargeParser(ABC):
    """
    Abstract base class for charge parsers.
    """

    def __init__(self, dirname: str | None):
        self.dirname = dirname

    @property
    @abstractmethod
    def filename(self) -> str:
        """
        Return the name of the charge file.
        """
        return ""

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        return f"{self.__class__.__name__}(filename={self.filename})"

    # define a an abstract property that contains the name of the method
    @property
    @abstractmethod
    def method(self) -> str:
        """
        Return the method name.
        """
        return ""

    def get_file_path(self, cid: Path) -> Path:
        """
        Get the full path to the charge file.
        """
        if self.dirname:
            return cid / self.dirname / self.filename
        else:
            return cid / self.filename

    @abstractmethod
    def parse_charges(self, cid: Path) -> dict[int, float]:
        """
        Parse the charges from the files in the directory.
        """
        return {}


# Implement specific parsers for each charge directory
class TbliteParser(ChargeParser):
    """
    Parser for tblite charges.
    """

    @property
    def filename(self) -> str:
        return "tblite.json"

    def parse_charges(self, cid: Path) -> dict[int, float]:
        """
        The charge file looks as follows:
        {
          "version":  "0.3.0",
          "charges": [
            -5.3002100975529953E-01,
             1.1785065999162758E-01,
             ...
             1.0694730077050983E-01
          ]
        }
        """
        charges: dict[int, float] = {}
        chargefile = self.get_file_path(cid)
        with chargefile.open("r", encoding="utf-8") as f:
            charge_data = json.load(f)
            charges = {
                i + 1: float(charge) for i, charge in enumerate(charge_data["charges"])
            }

        return charges


class XtbParser(ChargeParser):
    """
    Parser for xTB charges.
    """

    @property
    def filename(self) -> str:
        return "xtbout.json"

    def parse_charges(self, cid: Path) -> dict[int, float]:
        """
        The charge file looks as follows:
        {
           "total energy":         -50.44446907,
           "HOMO-LUMO gap / eV":           4.58016104,
           "electronic energy":         -51.12528068,
           "dipole / a.u.": [    -0.85296406,     0.22431893,    -0.28731140],
           "partial charges": [
               -0.56254128,
                0.24017395,
               -0.02909668,
               -0.07455592,
               -0.07071172,
               -0.04560688,
               -0.04913459,
                0.00329375,
               -0.04736118,
                0.00388507,
               -0.02524975,
               -0.09087441,
               ...

        """
        charges: dict[int, float] = {}
        chargefile = self.get_file_path(cid)
        with chargefile.open("r", encoding="utf-8") as f:
            charge_data = json.load(f)
            charges = {
                i + 1: charge for i, charge in enumerate(charge_data["partial charges"])
            }

        return charges


class OrcaParser(ChargeParser):
    """
    Parser for ORCA Hirshfeld charges.
    """

    @property
    def filename(self) -> str:
        return "orca.out"

    def parse_charges(self, cid: Path, numatoms: int | None = None) -> dict[int, float]:
        """
        The charge file looks as follows (this is the cutout of a larger file):
        ...
        ------------------
        HIRSHFELD ANALYSIS
        ------------------

        Total integrated alpha density =     90.999981137
        Total integrated beta density  =     90.999981137

          ATOM     CHARGE      SPIN
           0 F   -0.105360    0.000000
           1 F   -0.110575    0.000000
           ...
          39 H    0.048434    0.000000

          TOTAL   0.000038    0.000000

        -------
        TIMINGS
        -------

        Total SCF time: 0 days 1 hours 34 min 4 sec

        """
        charges: dict[int, float] = {}
        chargefile = self.get_file_path(cid)
        with chargefile.open("r", encoding="utf-8") as f:
            lines = f.readlines()
            found = False
            for k, line in enumerate(lines):
                if "HIRSHFELD ANALYSIS" in line:
                    found = True
                    break
            if not found:
                raise ValueError("No Hirshfeld analysis found in file")
            if numatoms is not None:
                for line in lines[k + 7 : k + 7 + numatoms]:  # pylint: disable=W0631
                    index, _, charge, _ = line.strip().split()
                    charges[int(index) + 1] = float(charge)
            else:
                raise ValueError("Number of atoms not provided")
                # TODO: Implement a way to read the number of atoms from the file

        return charges


class EspalomaParser(ChargeParser):
    """
    Parser for the Espaloma charges.

    The charge file looks as follows:
    -0.413498 -0.145383 -0.613541 -0.325130 -0.362620  0.106125 -0.005218 -0.005218
    -0.051348 -0.051348 -0.071277  0.210935  0.743412  0.057223 -0.083591  0.776921
    0.199649  0.312315 -0.116748 -0.161659
    """

    @property
    def filename(self) -> str:
        return "espaloma.charge"

    @property
    def method(self) -> str:
        return "EspalomaCharge"

    def parse_charges(
        self,
        cid: Path,
        atomic_numbers: list[int] | None = None,
        total_charge: float = 0.0,
    ) -> dict[int, float]:
        # variable to store the charges
        charges: dict[int, float] = {}

        # catch improper input
        if atomic_numbers is None:
            raise ValueError("No atomic numbers provided!")
        if total_charge != 0.0:
            warnings.warn(
                f"Total charge of {total_charge} is not zero! The charges might be incorrect. Using dummy charges..."
            )
            charges = {i + 1: -10.0 for i in range(len(atomic_numbers))}
        chargefile = self.get_file_path(cid)

        try:
            with chargefile.open("r", encoding="utf-8") as f:
                lines = f.readlines()
                chargelist: list[float] = []
                for i, line in enumerate(lines):
                    chargelist.extend([float(charge) for charge in line.split()])

                # Check if the number of charges is equal to the number of heavy (non-hydrogen) atoms
                number_heavy_atoms = len(
                    [atnum for atnum in atomic_numbers if atnum != 1]
                )
                if len(chargelist) != number_heavy_atoms:
                    warnings.warn(
                        f"Number of Espaloma charges ({len(chargelist)}) does not match the number of heavy atoms ({number_heavy_atoms})"
                        + f" in {cid.name}"
                    )
                    charges = {i + 1: -10.0 for i in range(len(atomic_numbers))}
                else:
                    # Go though all atoms of the molecule
                    # If it is not a hydrogen: assign the charge to the atom from the next entry in the file
                    # Else: append a dummy value to the charges
                    for i, atnum in enumerate(atomic_numbers):
                        if atnum != 1:
                            charges[i + 1] = float(chargelist.pop(0))
                        else:
                            charges[i + 1] = -10.0
        except FileNotFoundError:
            # raise a warning if the file does not exist
            warnings.warn(f"File {chargefile} not found! Using dummy charges...")
            # use dummy values for the charges
            charges = {i + 1: -10.0 for i in range(len(atomic_numbers))}

        return charges


class HirsheldAPCParser(ChargeParser):
    """
    Parser for the Hirshfeld charges in apc format (like in acqm).

    The charge file looks as follows
    ╰─(chempy) ⠠⠵ cat hirshfeld.apc
    -1.289991
    -0.777432
    -0.932609
    """

    @property
    def method(self) -> str:
        return "HirshfeldAPC"

    @property
    def filename(self) -> str:
        return "hirshfeld.apc"

    def parse_charges(self, cid: Path) -> dict[int, float]:
        """
        The charge file looks as follows:
        -1.289991
        -0.777432
        -0.932609
        """
        charges: dict[int, float] = {}
        chargefile = self.get_file_path(cid)
        with chargefile.open("r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                charges[i + 1] = float(line.strip())

        return charges


class EeqParser(TbliteParser):
    """
    Parser for the EEQ charges.
    """

    @property
    def method(self) -> str:
        return "EEQ"


class EeqMultichargeParser(EeqParser):
    """
    Parser for the EEQ charges.
    """

    @property
    def filename(self) -> str:
        return "multicharge.json"


class CehV2Parser(TbliteParser):
    """
    Parser for the CEH charges.
    """

    @property
    def method(self) -> str:
        return "CEH-v2"


class CehV1Parser(TbliteParser):
    """
    Parser for the CEH charges.
    """

    @property
    def method(self) -> str:
        return "CEH-v1"


class Gfn1XtbParser(XtbParser):
    """
    Parser for the GFN1-xTB charges.
    """

    @property
    def method(self) -> str:
        return "GFN1-xTB"


class Gfn2XtbParser(XtbParser):
    """
    Parser for the GFN2-xTB charges.
    """

    @property
    def method(self) -> str:
        return "GFN2-xTB"


class WB97MVParser(OrcaParser):
    """
    Parser for the wB97M-V charges.
    """

    @property
    def method(self) -> str:
        return "wB97M-V"


# Dictionary to map directory names to their respective parser classes
PARSER_CLASSES: dict[type[ChargeParser], dict[str, str | None]] = {
    EeqParser: {"dirname": "eeq"},
    CehV2Parser: {"dirname": "ceh-v2"},
    CehV1Parser: {"dirname": "ceh-v1"},
    Gfn1XtbParser: {"dirname": "gfn1xtb"},
    Gfn2XtbParser: {"dirname": "gfn2xtb"},
    WB97MVParser: {"dirname": "wB97M-V_def2-TZVPPD_Hirshfeld"},
    EspalomaParser: {"dirname": "espaloma"},
    HirsheldAPCParser: {"dirname": None},
    EeqMultichargeParser: {"dirname": "eeq"},
}


def read_compounds_csv(file_path: Path):
    """
    Read the CSV file containing the paths to the directories of compounds.
    It looks as follows:
    CID	Name	Structure Path	Number of Atoms	Total Charge	HOMO-LUMO gap
    6695	[9-(2-carboxyphenyl)-6-(diethylamino)xanthen-3-ylidene]-diethylazanium	/Users/marcelmueller/source/randommolecules/example/6695/6695_opt.xyz	64	1	1.561
    43467	diphenyl propan-2-yl phosphate	/Users/marcelmueller/source/randommolecules/example/43467/43467_opt.xyz	37	0	4.520
    44961	[4-(dimethylamino)-3-methylphenyl] N,N-dimethylcarbamate	/Users/marcelmueller/source/randommolecules/example/44961/44961_opt.xyz	34	0	3.511
    60641	4-[2-(4-methylphenyl)sulfonylhydrazinyl]-4-oxobut-2-enoic acid	/Users/marcelmueller/source/randommolecules/example/60641/60641_opt.xyz	31	0	2.193
    76470	bis(2-methylpropyl)carbamothioylsulfanyl N,N-bis(2-methylpropyl)carbamodithioate	/Users/marcelmueller/source/randommolecules/example/76470/76470_opt.xyz	60	0	2.258
    94837	(4-hydroxy-2-methyl-5-propan-2-ylphenyl) thiocyanate	/Users/marcelmueller/source/randommolecules/example/94837/94837_opt.xyz	273.760
    ...
    """
    # read the CSV into a pandas dataframa
    return pd.read_csv(file_path, sep="\t")


def atomic_numbers_from_xyz(file_path: Path) -> list[int]:
    """
    Read the atomic numbers from an XYZ file.
    """
    with file_path.open("r", encoding="utf-8") as f:
        lines = f.readlines()
        num_atoms = int(lines[0])
        atomic_numbers: list[int] = [
            PSE_NUMBERS[line.split()[0].lower()] for line in lines[2 : num_atoms + 2]
        ]
    return atomic_numbers


def get_compounds_from_directory(search_key: str) -> pd.DataFrame:
    """
    Get the compounds from the directory structure.
    """
    # open the file and read all lines
    compounds = pd.DataFrame()
    # search current working directory for directories containing the search_key
    for calc_dir in Path(".").glob(f"**/{search_key}"):
        calc_dir = calc_dir.resolve()
        # if one of the directories in the path starts with a dot, skip it (so we avoid using hidden directories)
        if any(part.startswith(".") for part in calc_dir.parts):
            continue
        # the compound directory is the parent directory of the calculation directory
        directory = calc_dir.parent
        struc_files = list(directory.glob("*.xyz"))
        if len(struc_files) > 1:
            raise ValueError(f"More than one structure file found in {directory}")
        else:
            struc_file = struc_files[0]
        atnums = atomic_numbers_from_xyz(struc_file)
        # get the charge from the .CHRG file
        mol_charge = 0.0
        try:
            with open(directory / ".CHRG", "r", encoding="utf-8") as f:
                mol_charge = float(f.read())
        except FileNotFoundError:
            warnings.warn(f"File {directory / '.CHRG'} not found! Using zero charge...")
        # append the compound to the dataframe
        compounds = pd.concat(
            [
                compounds,
                pd.DataFrame(
                    {
                        "CID": [directory.name],
                        "Name": [directory.name],
                        "Structure Path": [struc_file],
                        "Number of Atoms": [len(atnums)],
                        "Total Charge": [mol_charge],
                    }
                ),
            ]
        )
    # set the CID as the index
    compounds.set_index("CID", inplace=True)

    return compounds


def get_args() -> argparse.Namespace:
    """
    Parse the command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract atomic charges from directories."
    )
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="Path to the CSV file containing the Compounds dataframe.",
        required=False,
        default=None,
    )
    parser.add_argument(
        "--min-num-atoms",
        type=int,
        help="Minimum number of atoms in the structure file.",
        required=False,
    )
    return parser.parse_args()


def main():
    """
    Main function to extract the atomic charges from the directories and write them to a new CSV file.
    """

    # Parse the command line arguments
    args = get_args()

    # Define the paths and files
    if args.file:
        # Check if the files exist
        if not Path(args.file).exists():
            raise FileNotFoundError(f"File {args.file} not found!")
        dir_list_file = Path(args.file).resolve()
        if dir_list_file.suffix == ".csv":
            compounds = read_compounds_csv(dir_list_file)
            # for the first compound (assuming the same behavior for all), check if the structure path exists
            if not Path(compounds["Structure Path"].iloc[0]).exists():
                # raise a warning if the file does not exist
                warnings.warn(
                    f"File {compounds['Structure Path'].iloc[0]} not found! Replacing with local path structure..."
                )
                # check if a directory with name of CID exists
                tmp_cid = str(compounds["CID"].iloc[0])
                if Path(tmp_cid + "/" + tmp_cid + "_opt.xyz").exists():
                    # replace the structure path with the local path structure
                    compounds["Structure Path"] = compounds["CID"].apply(
                        lambda x: Path(str(x) + "/" + str(x) + "_opt.xyz").resolve()
                    )
                else:
                    raise FileNotFoundError(
                        f"Directory {compounds['CID'].iloc[0]} not found!"
                    )
        else:
            raise ValueError(f"File {dir_list_file} is not a CSV file!")
    else:
        # dirname of first parser in PARSER_CLASSES is used as search key
        search_key = PARSER_CLASSES[list(PARSER_CLASSES.keys())[0]]["dirname"]
        warnings.warn(
            f"Using all directories containing '{search_key}' as compound directories..."
        )
        compounds = get_compounds_from_directory(search_key)
    # in both cases, write the compounds to a new CSV file
    compounds.to_csv("compounds_rev.csv", index=False, sep="\t")

    # Only consider compounds above the minimum number of atoms
    if args.min_num_atoms:
        compounds = compounds[compounds["Number of Atoms"] >= args.min_num_atoms]
        compounds.to_csv("compounds_rev_reduced.csv", index=False, sep="\t")

    # Prepare a new pandas dataframe consisting of the columns "CID", "Method", "Charge", "Atom number", "Atom type"
    # the new data frame shall have the following structure (with Atom type being the ordinal number of the atom in the periodic table):
    # CID | Method | Charge | Atom number | Atom type
    # 6695 | "EEQ" | -0.53 | 1 | 6
    # 6695 | "EEQ" | 0.11 | 2 | 1
    # ...
    # 6695 | "GFN2" | 0.05 | 64 | 7
    # 43467 | "GFN2" | -0.43 | 1 | 8
    # 43467 | "GFN2" | 0.16 | 2 | 7
    # ...

    def process_directory(parser: ChargeParser, compound: pd.Series) -> pd.DataFrame:
        """
        Process the directory and extract the charges.
        """
        # get directory names from column "CID" in the dataframe row
        directory = Path(compound["Structure Path"]).parent
        atnums = atomic_numbers_from_xyz(Path(compound["Structure Path"]))
        if len(atnums) != compound["Number of Atoms"]:
            raise ValueError(
                f"Number of atoms in the structure file ({len(atnums)}) does not match the number of atoms in the dataframe ({compound['Number of Atoms']})"
                + f" in {directory.name}"
            )
        args: list = [directory]
        if isinstance(parser, OrcaParser):
            args.append(len(atnums))
        if isinstance(parser, EspalomaParser):
            args.append(atnums)
            try:
                args.append(float(compound["Total Charge"]))
            except IndexError as e:
                print(f"Total charge not found for {directory.name}")
                # print dataframe cutout for specific for debugging
                print(compounds.loc[compounds["CID"] == directory.name])
                raise e
        try:
            atcharges = parser.parse_charges(*args)
        except FileNotFoundError:
            warnings.warn(
                f"File '{parser.get_file_path(directory)}' for parser '{parser}' not found! Inserting 'NaN' charges..."
            )
            atcharges = {i + 1: float("NaN") for i in range(len(atnums))}
        # Check if the number of charges is equal to the number of atoms
        if len(atcharges) != compound["Number of Atoms"]:
            raise ValueError(
                f"Number of charges ({len(atcharges)}) does not match the number of atoms ({compound['Number of Atoms']})"
                + f" in {directory.name}"
            )
        data = {
            "CID": [directory.name] * compound["Number of Atoms"],
            "Method": [parser.method] * compound["Number of Atoms"],
            "Charge": list(atcharges.values()),
            "Atom number": list(range(1, compound["Number of Atoms"] + 1)),
            "Atom type": atnums,
        }
        return pd.DataFrame(data)

    charges = pd.DataFrame()
    total_iterations = len(PARSER_CLASSES) * len(compounds)
    with tqdm(total=total_iterations, desc="Parsed charges...") as pbar:
        for parser_class, config in PARSER_CLASSES.items():
            parser = parser_class(**config)
            for _, compound in compounds.iterrows():
                df = process_directory(parser, compound)
                charges = pd.concat([charges, df])
                pbar.update(1)

    # write the charges to a new CSV file
    charges.to_csv("charges.csv", index=False)


if __name__ == "__main__":
    main()
