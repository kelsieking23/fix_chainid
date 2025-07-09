import os
import argparse
from dataclasses import dataclass
from typing import Union
import sys

__version__ = 1.0

RESIDUES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HSD', 'HSE']

class StructureFile:

    def __init__(self, inp, filetype=None, atom_data=[]):
        self.inp = inp
        self.atoms = atom_data
        if self.inp is not None:
            if isinstance(self.inp, str):
                _, ext = os.path.splitext(self.inp)
                self.ext = ext[1:]
                if self.ext != 'gro':
                    self.atoms = [atom for atom in self.read()]
            else:
                self.ext = filetype
                self.atoms = [atom for atom in self.read()]
        else:
            self.ext = None
        # self.atoms = [atom for atom in self.read()]

    def _structureReader(self):
        if isinstance(self.inp, str):
            with open(self.inp, 'r') as f:
                for line in f:
                    yield line
        elif isinstance(self.inp, list):
            for line in self.inp:
                yield line
        else:
            raise ValueError('Cannot read from input')

    def _atomDataIterator(self):
        for atom in self.atoms:
            yield atom
        return self.atoms

    def read(self):
        '''
        Iterates through structure file and yeilds lines.
        '''
        if self.ext is None:
            return self._atomDataIterator()
        if self.ext == 'pdb':
            return self.pdb()


    def write(self, out):
        if self.ext == None:
            ext = os.path.splitext(out)[1][1:]
        else:
            ext = self.ext 
        if ext == 'pdb':
            self.write_pdb(self.atoms, out)
    def get_atom_data_pdb(self, line):
        pass
    
    def crd(self, writer=False):
        atom_index = 0
        residue_index = -1
        last_residue = None
        last_chain = None
        chain_index = -1
        chain = 'A'
        charge = ''
        box = (0,0,0)
        model = 1
        for line in self._structureReader():
            if (line.startswith('*')) or ('EXT' in line) or (line.strip() == ''):
                continue
            else:
                line_parts = line.strip().split()
                atom_number = line_parts[0]
                atom_name = line_parts[3]
                x = float(line_parts[4])*10
                y = float(line_parts[5])*10
                z = float(line_parts[6])*10
                segid = line_parts[7]
                residue_name = line_parts[2]
                residue_number = int(line_parts[1])
                residue_id = residue_name + residue_name + str(residue_number)
                if residue_id != last_residue:
                        residue_index += 1
                elem = atom_name[0]
                charge = ''
                if segid != last_chain:
                    chain_index += 1
                    if chain != 'A':
                        chain = chr(ord(chain) + 1)
                temp = 0.00
                occ = 0.00
                if not writer:
                    atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, chain, chain_index, 
                                    residue_number, residue_index, x, y, z, occ, temp, segid, elem, charge, model, box, 
                                    line, False, '', 'crd')
                    yield atom
                else:
                    atom = ['ATOM', str(atom_number), atom_name, residue_name, chain, str(residue_number), x, y, z, '1.00', '0.00', segid, elem]
                    yield atom

    def pdb(self):
        atoms = []
        model = 1
        last_chain = None
        last_residue = None
        atom_index = 0
        chain_index = -1
        residue_index = -1
        box = (0,0,0)
        for line in self._structureReader():
            line_parts = line.split()
            if len(line_parts) == 0:
                continue
            if line.startswith('ENDMDL'):
                model = model + 1
                atom_index = 0
                chain_index = 0
                residue_index = 0
                last_chain = None
                last_residue = None
            elif (line.startswith('ATOM')) or (line.startswith('HETATM')):
                atom_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:21].strip()
                chain = line[21]
                residue_number = int(line[22:26].strip())
                residue_id = residue_name + str(residue_number)
                if (residue_name + str(residue_number) != last_residue):
                    residue_index += 1
                if (last_chain != chain):
                    chain_index += 1
                if (chain == '') or (chain == ' '):
                    _chain_index = -1
                else:
                    _chain_index = chain_index
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                occ = float(line[54:60].strip())
                temp = float(line[60:66].strip())
                segid = line[72:76].strip()
                elem = line[76:78].strip()
                _charge = line.strip()[-1]
                if (_charge == '+') or (_charge == '-'):
                    charge = _charge
                else:
                    charge = ''
                atom = AtomData(atom_number, atom_index, atom_name, residue_name, residue_id, 
                                chain, _chain_index, residue_number, residue_index,
                                x, y, z, occ, temp, segid, elem, charge, model, box, line, True, line[0:6].strip(), 'pdb', self)
                atom_index += 1
                last_chain = chain
                last_residue = residue_name + str(residue_number)
                yield atom
            else:
                continue
        return atoms

    def gro(self):
        atoms = []
        atom_index = 0
        residue_index = 0
        i = 0
        for line in self._structureReader():
            if i < 2:
                i += 1
                continue
            try:
                residue_number = int(line[0:5].strip())
                residue_name = line[5:10].strip()
                atom_name = line[10:15].strip()
                atom_number = int(line[15:20].strip())
                residue_id = '{}{}'.format(residue_name, residue_number)
                x = float(line[20:29].strip())
                y = float(line[29:37].strip())
                z = float(line[37:46].strip())
                if (len(line) > 45):
                    vx = float(line[46:53].strip())
                    vy = float(line[53:61].strip())
                    vz = float(line[61:].strip())
                else:
                    vx = 0.0
                    vy = 0.0
                    vz = 0.0
                atom = GroAtomData(atom_number, atom_index, atom_name, residue_name, residue_id,
                                   residue_number, residue_index, x, y, z, vx, vy, vz, line)
                atom_index += 1
                residue_index += 1
                yield atom
            except:
                pass
            i += 1
        return atoms

    def sdf(self):
        return []
    
    def mol(self):
        return []

    @staticmethod
    def write_gro(atom_data, out, title='Title', box=[0.0, 0.0, 0.0]):
        with open(out, 'w') as f:
            f.write(f'{title}\n')
            f.write(f' {len(atom_data)}\n')
            for atom in atom_data:
                ld = [atom.residue_number, atom.residue_name, atom.atom_name, atom.atom_number, atom.x, atom.y, atom.z,
                    atom.vx, atom.vy, atom.vz]
                line = '{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(*ld)
                f.write(line)
            f.write('   {:.7f}   {:.7f}   {:.7f}\n'.format(*box))
    
    @staticmethod
    def write_pdb(atom_data, out):
        with open(out, 'w') as f:
            for atom in atom_data:
                f.write(atom.line)
            
    @classmethod
    def fromAtomData(cls, atom_data):
        return cls(inp=None, atom_data=atom_data)
    
    @classmethod
    def fromList(cls, lst, filetype='pdb'):
        return cls(inp=lst, filetype=filetype)

@dataclass
class GroAtomData:
    atom_number: int
    atom_index: int
    atom_name: str
    residue_name: str
    residue_id: str
    residue_number: int
    residue_index: int
    x: float
    y: float
    z: float
    vx: float
    vy: float
    vz: float
    line: str

    def update_line(self):
        ld = [self.residue_number, self.residue_name, self.atom_name, self.atom_number, self.x, self.y, self.z,
              self.vx, self.vy, self.vz]
        line = '{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.4f}{:>8.4f}{:>8.4f}\n'.format(*ld)
        self.line = line
        return self


class AtomData:

    def __init__(self, atom_number: int, atom_index: int, atom_name: str, residue_name: str, residue_id: str,
                 chain: str, chain_index: int, residue_number: int, residue_index: int, x: float, y: float, z: float,
                 occ:float, temp: float, segid: str, elem: str, charge: str, model: float, box: tuple, line: str, 
                 is_pdb: bool, pdb_label: str, ext: str, parent: Union[StructureFile, None]):
        self._atom_number = atom_number
        self.atom_index = atom_index
        self.atom_name = atom_name
        self.residue_name = residue_name
        self.residue_id = residue_id
        self.chain = chain
        self.chain_index = chain_index
        self._residue_number = residue_number
        self.residue_index = residue_index
        self.x = x
        self.y = y
        self.z = z
        self.occ = occ
        self.temp = temp
        self.segid = segid
        self.elem = elem
        self.charge = charge
        self.model = model
        self.box = box
        self._line = line
        self.is_pdb = is_pdb
        self.pdb_label = pdb_label
        self.ext = ext
        self.parent = parent
        self.update_dict: dict={
                'pdb':self._update_pdb,
                'crd':self._update_crd}
        self._cannonical = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'GLY', 'PRO', 
                            'SER', 'THR', 'ASN', 'GLN', 'CYS', 'HIS', 'HSD', 'ASP', 'GLU', 'ARG', 'LYS']
        
    def __str__(self):
        return f'<pymd.structure.structure_file.AtomData Object>: {self.atom_number} {self.atom_name} {self.residue_number} {self.residue_name} {self.chain}'

    def _update_parent(self):
        self._line = self.update_dict[self.ext]
        self.parent.atoms[self.atom_index] = self

    @property
    def atom_number(self):
        return self._atom_number
    
    @atom_number.setter
    def atom_number(self, number):
        try:
            self._atom_number = int(number)
        except Exception:
            raise ValueError(f'Error setting attribute "atom_number"')
        self._update_parent()
    
    @property
    def residue_number(self):
        return self._residue_number
    
    @residue_number.setter
    def residue_number(self, number):
        try:
            self._residue_number = int(number)
        except Exception:
            raise ValueError(f'Error setting attribute "residue_number"')
        self._update_parent()

    @property
    def line(self):
        if self.ext == 'pdb':
            self._line = self._update_pdb()
            if isinstance(self.parent, StructureFile):
                self.parent.atoms[self.atom_index] = self
        return self._line
    
    def _update_pdb(self):
        line = [self.pdb_label, self.atom_number, self.atom_name, self.residue_name, self.chain, self.residue_number, self.x, self.y, self.z,
                self.occ, self.temp, self.segid, self.elem]
        string = "{:6s}{:5d} {:^4s} {:^4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}{:>10s}  {:<3s}\n".format(*line)
               # "****  NNNNN xxxxI "
        return string

    def _update_crd(self):
        return ''
    
    def _in_rtp(self, rtp):
        with open(rtp, 'r') as f:
            for line in f:
                if line.startswith('[ '):
                    entry = line.strip()[2:-2]
                    if entry.upper() == self.residue_name.upper():
                        return True
        return False
    
    def mol_type(self, ff='guess'):
        '''
        at this moment super rudimentary, but just want to be able to ID if it is protein or not protein. 
        ff (str) : can be 'guess', which means just go off of _cannonical() which is hard coded
                   ff can also be: 
                    1) a path to a force field folder that contains .rtp files
                    2) a path to a specific .rtp file
        '''
        if ff == 'guess':
            if self.residue_name in self._cannonical:
                return 'protein'
            else:
                return 'other'
        else:
            if os.path.isdir(ff):
                if os.path.isfile(os.path.join(ff, 'aminoacids.rtp')):
                    if self._in_rtp(os.path.join(ff, 'aminoacids.rtp')):
                        return 'protein'
                for file in os.listdir(ff):
                    if file.endswith('rtp'):
                        if self._in_rtp(os.path.join(ff, file)):
                            return os.path.splitext(file)[0]
                return 'other'
            elif os.path.isfile(ff):
                if self._in_rtp(ff):
                    return True
                else:
                    return False
            else:
                raise ValueError('dont know!!! RIP')
            return 'other'


def extract_model(pdb):
    '''
    Extract models from cluster.pdb file
    '''
    contents = {}
    f = open(pdb, 'r')
    outputs = []
    mod = 1
    for line in f:
        if mod not in contents.keys():
            contents[mod] = []
        
        if line.startswith('ENDMDL'):
            contents[mod].append(line)
            mod += 1
        else:
            contents[mod].append(line)
    f.close()
    return contents

def fix_chainID(content_list, output=None, start=1, atom_start=1, renumber_chain=True):
    struct = StructureFile.fromList(content_list)
    new_contents = ['REMARK    CHAIN ID EDITED BY BY fix_chainid.py\n']
    res_counter = start
    atom_counter = atom_start
    last_chain_index = None
    last_res_id = None
    last_res_num = None
    chain_id = 'A'
    _real_chain_id = 'A'
    for i, ld in enumerate(struct.pdb()):
        # write non-atom strings
        if isinstance(ld, str):
            new_contents.append(ld)
            continue
        
        # get res id
        res_id = ld.residue_name + str(ld.residue_number)
        if ld.residue_name not in RESIDUES:
            new_contents.append(ld.line)
            continue

        #check for residue change
        if last_res_id is not None:
            if last_res_id != res_id:
                res_counter += 1

        if last_res_num is not None:
            if (ld.residue_number != last_res_num + 1) and (ld.residue_number != last_res_num):
                if chain_id != '':
                    chain_id = chr(ord(chain_id) + 1)
                    _real_chain_id = chr(ord(_real_chain_id) + 1)
                    if renumber_chain:
                        res_counter = start
                else:
                    _real_chain_id = chr(ord(_real_chain_id) + 1)
            elif (ld.residue_id != last_res_id) and (ld.residue_name not in RESIDUES):
                if chain_id != '':
                    chain_id = chr(ord(chain_id) + 1)
                    _real_chain_id = chr(ord(_real_chain_id) + 1)
                    if renumber_chain:
                        res_counter = start
                else:
                    _real_chain_id = chr(ord(_real_chain_id) + 1)

        # check if residue counter > 9999
        if res_counter > 9999:
            res_counter = 0

        # check if atom_counter > 99999
        if atom_counter > 99999:
            atom_counter = 0

        # check chain index
        if ld.chain_index != last_chain_index:
            res_counter = start

        if ld.residue_name not in RESIDUES:
            chain_id = ''
        # format and save line for writing
        line_parts = ['ATOM', str(atom_counter), ld.atom_name, ld.residue_name, chain_id, str(res_counter), 
              ld.x, ld.y, ld.z, ld.occ, ld.temp, ld.segid, ld.elem, ld.charge]
        if len(ld.residue_name) <= 3:
            if len(ld.atom_name) > 3:
                line = '{:<4s}{:>7s} {:<4s} {:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>1.2f}  {:>1.2f}{:>10s} {}{}\n'.format(*line_parts)
            else:
                line = '{:<4s}{:>7s}  {:<4s}{:>3s} {:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>1.2f}  {:>1.2f}{:>10s} {}{}\n'.format(*line_parts)
        else:
            if len(ld.atom_name) > 3:
                line = '{:<4s}{:>7s} {:<4s} {:>3s}{:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>1.2f}  {:>1.2f}{:>10s} {}{}\n'.format(*line_parts)
            else:
                line = '{:<4s}{:>7s}  {:<4s}{:>3s}{:1s}{:>4s}    {:>8.3f}{:>8.3f}{:>8.3f}  {:>1.2f}  {:>1.2f}{:>10s} {}{}\n'.format(*line_parts)

        new_contents.append(line)

        # update last data
        last_chain_index = ld.chain_index
        last_res_id = res_id
        last_res_num = ld.residue_number
        atom_counter += 1

    # write file
    if output is None:
        output = filename
    f = open(output, 'w')
    for line in new_contents:
        f.write(line)
    f.close() 


if __name__ == '__main__':
    description = 'Script to fix GROMACS-written .pdb files that do not include chain IDs.\n--------------------------------------------------------------------------'
    epilog = '--------------------------------------------------------------------------\n'
    epilog = epilog + 'Usage examples:\n1. Overwrite input .pdb with correct chain information:\n>> python fix_chainid.py incorrect_chain.pdb\n\n'
    epilog = epilog + '2. Write a new .pdb file with correct chain information:\n>> python fix_chainid.py incorrect_chain.pdb -o corrected_file.pdb\n\n'
    epilog = epilog + '3. Extract and correct the first 3 models from a multi-model .pdb:\n >>python fix_chainid.py incorrect_multimodel.pdb -extract 3'
    parser = argparse.ArgumentParser(description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inp', type=str, help='(positional argument, required) input .pdb file')
    parser.add_argument('-extract', type=int, default=1, help='(optional, default=1) number of models to extract, starting from model 1.')
    parser.add_argument('-o', type=str, default='', help='(optional) output filename. if not specified, will overwrite the input file. if extracting multiple files, will add numbers to the file name. if -o not specified, will write to the directory the input file is located. ')
    parser.add_argument('-resstart', type=int, default=1, help='(optional, default=1) renumber residues, with residue numbering starting with specified number')
    parser.add_argument('-atomstart', type=int, default=1, help='(optional, default=1) renumber atoms, with atom numbering starting with specified number')
    parser.add_argument('-v', '--version', action='store_true', help='print version and quit')
    args = parser.parse_args()

    if args.version:
        print(__version__)
        sys.exit(0)
    if args.o == '':
        output = args.inp
    else:
        output = args.o
    models = extract_model(args.inp)
    if (isinstance(args.extract, int))  or (isinstance(args.extract, float)):
        try:
            n_models = int(args.extract)
        except:
            raise ValueError('{} not recognized for n_models. -n_models must be an integer or all')
        
        model_keys = range(1, n_models+1)
    else:
        raise ValueError('{} not recognized for n_models. -n_models must be an integer or all')
    if len(model_keys) > 1:
        for key in model_keys:
            contents = models[key]
            path = os.path.dirname(output)
            base = os.path.splitext(os.path.basename(output))[0]
            out = os.path.join(path, '{}.{}.pdb'.format(base, str(key).zfill(3)))
            fix_chainID(contents, out, start=args.resstart, atom_start=args.atomstart)
    else:
        fix_chainID(models[1], output, start=args.resstart, atom_start=args.atomstart)
