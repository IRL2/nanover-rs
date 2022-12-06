#!/usr/bin/env python
"""
Read components.cif and write a rust module describing bond templates.
"""

from typing import Iterable, NamedTuple, Optional
from enum import Enum
from dataclasses import dataclass, field
import sys

ORDER_CORRESPONDANCE = {'SING': 1.0, 'DOUB': 2.0, 'TRIP': 3.0}


class Context(Enum):
    IDLE = 1
    LOOP_KEY = 2
    LOOP = 3


class BondTemplate(NamedTuple):
    from_atom: str
    to_atom: str
    order: float

@dataclass
class Template:
    name: Optional[str] = None
    bonds: list = field(default_factory=lambda: list())


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


def escape_name(name: str) -> str:
    return name.replace('"', '\\"')

def write_rust_module(data_block: list[Template]):
    header = """use phf::phf_map;

pub struct BondTemplate<'a> {
    pub from: &'a str,
    pub to: &'a str,
    pub order: f32,
}
pub static BOND_TEMPLATES: phf::Map<&str, &'static [BondTemplate]> = phf_map! {"""
    footer = "};"
    item = '    "{name}" => &[\n{templates}\n   ],'
    bond = '        BondTemplate{{from: "{}", to: "{}", order: {}}}'

    print(header)
    for block in data_block:
        if block.name is None or len(block.name) != 3:
            continue
        bonds = ",\n".join([
            bond.format(escape_name(b.from_atom), escape_name(b.to_atom), b.order)
            for b in block.bonds
        ])
        print(item.format(name=block.name, templates=bonds))
    print(footer)


def main():
    path = sys.argv[1]
    with open(path) as infile:
        templates = parse_components(infile)
    write_rust_module(templates)


if __name__ == '__main__':
    main()