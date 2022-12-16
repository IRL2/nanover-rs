#!/usr/bin/env python
"""
Read components.cif and write a rust module describing bond templates.
"""

from typing import Iterable, NamedTuple, Optional
from enum import Enum, IntEnum
from dataclasses import dataclass, field
import sys

ORDER_CORRESPONDANCE = {'SING': 1.0, 'DOUB': 2.0, 'TRIP': 3.0}


class Context(Enum):
    IDLE = 1
    LOOP_KEY = 2
    LOOP = 3


class ResidueType(IntEnum):
    NON_POLYMER = 0
    PEPTIDE = 1
    DNA = 2
    RNA = 3
    UNSUPPORTED = 4


class BondTemplate(NamedTuple):
    from_atom: str
    to_atom: str
    order: float

@dataclass
class Template:
    name: Optional[str] = None
    residue_type: ResidueType = ResidueType.UNSUPPORTED
    raw_residue_type: str = ""
    bonds: list[BondTemplate] = field(default_factory=lambda: list())


def combine_line(line: str, keys: Iterable[str]) -> dict[str, str]:
    return dict(zip(keys, line.split()))


def do_bond_line(line: str, loop_keys: Iterable[str], template: Template):
    combined = combine_line(line, loop_keys)
    from_atom = combined['atom_id_1']
    to_atom = combined['atom_id_2']
    order = ORDER_CORRESPONDANCE[combined['value_order']]
    bond = BondTemplate(from_atom, to_atom, order)
    template.bonds.append(bond)


def parse_components(lines: Iterable[str]) -> list[BondTemplate]:
    data_blocks = []
    context = Context.IDLE
    loop_keys: list[str] = []
    loop_name: Optional[str] = None
    for line in lines:
        tokens = line.strip().split(maxsplit=1)
        if not tokens:
            continue
        first = tokens[0]
        if context == Context.IDLE:
            if first.startswith('data_'):
                data_blocks.append(Template())
            elif first == 'loop_':
                loop_keys = []
                loop_name = None
                context = Context.LOOP_KEY
            elif first == '_chem_comp.id':
                if not data_blocks:
                    raise RuntimeError('Data item outside of a data block.')
                if len(tokens) != 2:
                    raise RuntimeError('Unexpected number of tokens.')
                template = data_blocks[-1]
                template.name = tokens[1]
            elif first == '_chem_comp.type':
                value = tokens[1].upper()
                if 'PEPTIDE' in value:
                    residue_type = ResidueType.PEPTIDE
                elif 'DNA' in value:
                    residue_type = ResidueType.DNA
                elif 'RNA' in value:
                    residue_type = ResidueType.RNA
                elif value == 'NON-POLYMER':
                    residue_type = ResidueType.NON_POLYMER
                else:
                    residue_type = ResidueType.UNSUPPORTED
                template.residue_type = residue_type
                template.raw_residue_type = value
        elif context == context.LOOP_KEY:
            if first.startswith('_'):
                name, key = first.split('.', maxsplit=1)
                if loop_name is not None and loop_name != name:
                    raise RuntimeError(f'Inconsistent keys in a loop. Found {name}, expected {loop_name}.')
                loop_keys.append(key)
                loop_name = name
            elif first.startswith('#'):
                context = Context.IDLE
            else:
                if loop_name == '_chem_comp_bond':
                    if not data_blocks:
                        raise RuntimeError('Data loop outside of a data block.')
                    template = data_blocks[-1]
                    do_bond_line(line, loop_keys, template)
                context = Context.LOOP
        elif context == context.LOOP:
            if first.startswith('#'):
                context = Context.IDLE
            elif loop_name == '_chem_comp_bond':
                if not data_blocks:
                    raise RuntimeError('Data loop outside of a data block.')
                template = data_blocks[-1]
                do_bond_line(line, loop_keys, template)
    return data_blocks


def integrate_exceptions(data_blocks: Iterable[Template]):
    exceptions = {
        ResidueType.PEPTIDE: {
            "H": "H1",
            "OXT": "O2",
        },
    }
    for block in data_blocks:
        for bond in block.bonds:
            from_atom = exceptions.get(block.residue_type, {}).get(bond.from_atom, bond.from_atom)
            to_atom = exceptions.get(block.residue_type, {}).get(bond.to_atom, bond.to_atom)
            new_bond = BondTemplate(from_atom, to_atom, bond.order)
            if bond != new_bond:
                block.bonds.append(new_bond)


def add_extra_bonds(data_blocks: Iterable[Template]):
    for block in data_blocks:
        if block.residue_type == ResidueType.PEPTIDE:
            new_bond = BondTemplate("N", "H3", 1.0)
            block.bonds.append(new_bond)


def prune_components(data_blocks: Iterable[Template]) ->  list[Template]:
    keep = (
        'A', 'ACE', 'AIB', 'ALA', 'ARG', 'ASN', 'ASP', 'C', 'CYS',
        'DA', 'DC', 'DG', 'DT', 'FOR', 'G', 'GLN', 'GLU', 'GLY',
        'HIS', 'HOH', 'ILE', 'LEU', 'LYS', 'MET', 'ORN', 'PCA', 'PHE',
        'PRO', 'SER', 'THR', 'TRP', 'TYR', 'U', 'VAL',
    )
    return [block for block in data_blocks if block.name in keep]


def write_components_as_bytes(outfile, data_blocks):
    for block in data_blocks:
        if len(block.name) != 3:
            continue
        froms = (len(bond.from_atom) <= 3 for bond in block.bonds)
        tos = (len(bond.to_atom) <= 3 for bond in block.bonds)
        if not all(froms) or not all(tos):
            continue
        if not block.bonds:
            continue

        chunk = b""
        chunk += block.name.encode()
        chunk += int(block.residue_type).to_bytes(1, 'little')
        chunk += len(block.bonds).to_bytes(2, 'little')
        outfile.write(chunk)
        for bond in block.bonds:
            chunk = b""
            chunk += "{:3s}".format(bond.from_atom).encode('ascii')
            chunk += "{:3s}".format(bond.to_atom).encode('ascii')
            chunk += int(bond.order).to_bytes(1, 'little')
            if len(chunk) != 7:
                raise RuntimeError(f"Wrong size for this block {block}. {len(chunk)} bytes instead of 7 bytes.")
            outfile.write(chunk)


def main():
    path = sys.argv[1]
    outpath = sys.argv[2]
    with open(path) as infile:
        templates = parse_components(infile)
    integrate_exceptions(templates)
    add_extra_bonds(templates)
    templates = prune_components(templates)
    with open(outpath, 'wb') as outfile:
        write_components_as_bytes(outfile, templates)
    
    print([block.name for block in templates])


if __name__ == '__main__':
    main()