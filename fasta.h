#pragma once
#include <string>
#include <vector>
#include <unordered_map>

// should not use namespace here because it forces that choice on any file including this one
 
struct FastaRecord {
    std::string header;
    std::string sequence;
};

// Parse the file and return the record vectors
std::vector<FastaRecord> parseFasta(const std::string& filename);

// Build and return codon table
const std::unordered_map<std::string, std::string> CODONTABLE = {
    
        {"AUU", "ILE (I)"},
		{"AUC", "ILE (I)"},
		{"AUA", "ILE (I)"},
		{"AUG", "MET (M)"},
		
		{"ACU", "THR (T)"},
		{"ACC", "THR (T)"},
		{"ACA", "THR (T)"},
		{"ACG", "THR (T)"},
		
		{"AAU", "ASN (N)"},
		{"AAC", "ASN (N)"},
		{"AAA", "LYS (K)"},
		{"AAG", "LYS (K)"},
		
		{"AGU", "SER (S)"},
		{"AGC", "SER (S)"},
		{"AGA", "ARG (R)"},
		{"AGG", "ARG (R)"},

		{"GUU", "VAL (V)"},
		{"GUC", "VAL (V)"},
		{"GUA", "VAL (V)"},
		{"GUG", "VAL (V)"},
		
		{"GCU", "ALA (A)"},
		{"GCC", "ALA (A)"},
		{"GCA", "ALA (A)"},
		{"GCG", "ALA (A)"},
		
		{"GAU", "ASP (D)"},
		{"GAC", "ASP (D)"},
		{"GAA", "GLU (E)"},
		{"GAG", "GLU (E)"},
		
		{"GGU", "GLY (G)"},
		{"GGC", "GLY (G)"},
		{"GGA", "GLY (G)"},
		{"GGG", "GLY (G)"},

		{"UUU", "PHE (F)"},
		{"UUC", "PHE (F)"},
		{"UUA", "LEU (L)"},
		{"UUG", "LEU (L)"},
	
		{"UCU", "SER (S)"},
		{"UCC", "SER (S)"},
		{"UCA", "SER (S)"},
		{"UCG", "SER (S)"},

		{"UAU", "TYR (Y)"},
		{"UAC", "TYR (Y)"},
		{"UAA", "STOP"},
		{"UAG", "STOP"},
		
		{"UGU", "CYS (C)"},
		{"UGC", "CYS (C)"},
		{"UGA", "STOP"},
		{"UGG", "TRP (W)"},
		
		{"CUU", "LEU (L)"},
		{"CUC", "LEU (L)"},
		{"CUA", "LEU (L)"},
		{"CUG", "LEU (L)"},

		{"CCU", "PRO (P)"},
		{"CCC", "PRO (P)"},
		{"CCA", "PRO (P)"},
		{"CCG", "PRO (P)"},

		{"CAU", "HIS (H)"},
		{"CAC", "HIS (H)"},
		{"CAA", "GLN (Q)"},
		{"CAG", "GLN (Q)"},

		{"CGU", "ARG (R)"},
		{"CGC", "ARG (R)"},
		{"CGA", "ARG (R)"},
		{"CGG", "ARG (R)"},

    };

// Convert DNA to RNA
std::string transcribeDNA(const std::string& sequence);