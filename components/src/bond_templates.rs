use std::collections::HashMap;
use std::str;

pub enum ResidueType {
    NonPolymer,
    Peptide,
    DNA,
    RNA,
    Unsupported,
}

pub struct ResidueTemplate<'a> {
    pub residue_type: ResidueType,
    pub bonds: Vec<BondTemplate<'a>>,
}

#[derive(Debug)]
pub struct BondTemplate<'a> {
    pub from: &'a str,
    pub to: &'a str,
    pub order: f32,
}


pub fn get_bond_templates<'a>() -> HashMap<&'a str, ResidueTemplate<'a>> {
    let bytes_per_bond = 7;
    let bytes_in_header = 6;
    let bytes = include_bytes!("components.dat");
    let mut templates = HashMap::new();
    let mut index = 0;
    while index < bytes.len() {
        let name = str::from_utf8(&bytes[index..(index + 3)]).unwrap();
        let residue_type = match bytes[index + 4] {
            0 => ResidueType::NonPolymer,
            1 => ResidueType::Peptide,
            2 => ResidueType::DNA,
            3 => ResidueType::RNA,
            _ => ResidueType::Unsupported
        };
        let number_of_bonds: usize = u16::from_le_bytes([bytes[index + 4], bytes[index + 5]]).into();
        let bonds: Vec<BondTemplate> = bytes[(index + bytes_in_header)..(index + bytes_in_header + number_of_bonds * bytes_per_bond)]
            .chunks(7)
            .map(|c| {BondTemplate{
                from: str::from_utf8(&c[0..3]).unwrap(),
                to: str::from_utf8(&c[3..6]).unwrap(),
                order: c[6].into(),
            }})
            .collect();
        templates.insert(name, ResidueTemplate{residue_type, bonds});
        index += bytes_in_header + number_of_bonds * bytes_per_bond;
    }
    templates
}