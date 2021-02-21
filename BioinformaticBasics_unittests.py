import unittest
from BioinformaticsBasics import Bioinfo
from BioinformaticsBasics import Bioinfo_DNA
from BioinformaticsBasics import Bioinfo_Protein


#Variables for testing
test_tuple = ("AAT","CAG","TTA","CGC")
test_list = ["AAT","CAG","TTA","CGC"]
test_dict = {1: "ich", 2:"bin", 3:"ein", 4:"dict"}
test_number = 1
test_string = "AACTTAGACGACATGAGGATCAC"
test_fasta_dna = r".\Unittests_testfiles\testfile_dna.fasta"
test_fasta_empty = r".\Unittests_testfiles\testfile_empty.fasta"
test_fasta_protein = r".\Unittests_testfiles\testfile_protein.fasta"
test_csv = r".\Unittests_testfiles\testfile_dna.csv"
test_txt_dna = r".\Unittests_testfiles\testfile_dna.txt"
test_txt_lowercase = r".\Unittests_testfiles\testfile_lowercase.txt"
test_txt_protein = r".\Unittests_testfiles\testfile_protein.txt"
test_txt_nonsense = r".\Unittests_testfiles\testfile_nonsense.txt"
test_dict = {"T":"U"}


class UnitTests(unittest.TestCase):

# Unittests for Bioinfo.__init__() methode:
# checks if 1) argument is string and 2) string is path for .fasta or .txt file 
# 
    def test_Bioinfo_Input_csv(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_csv)
        self.assertEqual(cm.exception.code, 1)
            
    def test_Bioinfo_Input_number(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_number)
        self.assertEqual(cm.exception.code, 1)

    def test_Bioinfo_Input_string(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_string)
        self.assertEqual(cm.exception.code, 1)

    def test_Bioinfo_Input_list(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_list)
        self.assertEqual(cm.exception.code, 1)

    def test_Bioinfo_Input_dict(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_dict)
        self.assertEqual(cm.exception.code, 1) 

    def test_Bioinfo_Input_fasta(self):
        self.assertTrue(Bioinfo(test_fasta_dna))

    def test_Bioinfo_Input_text(self):
        self.assertTrue(Bioinfo(test_txt_dna))


# Unittest for Bioinfo.sequence_identifer():
# Checks if string is 1) empty, 2) includes lowercases 3) is DNA or 4) is Protein 

    def test_sequence_identifer_empty(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_empty)
            bio.sequence_identifer()
        self.assertEqual(cm.exception.code, 1)

    def test_sequence_identifer_lower_case(self):
        with self.assertRaises(SystemExit):
            bio = Bioinfo(test_txt_lowercase)
            self.assertFalse(bio.sequence_identifer())

    def test_sequence_identifer_dna(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.sequence_identifer())

    def test_sequence_identifer_protein(self):
        bio = Bioinfo(test_fasta_protein)
        self.assertTrue(bio.sequence_identifer())

    def test_sequence_identifer_nonsense(self):
        bio = Bioinfo(test_txt_nonsense)
        self.assertFalse(bio.sequence_identifer())


# Unittests for Bioinfo.count_characters():
# Count amount of specific characters in string

    def test_count_characters_search_for_Strings(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.count_characters("A", "T", "G", "C", "AAA"))

    def test_count_characters_search_for_not_Strings(self):
        with self.assertRaises(TypeError) as cm:
            bio = Bioinfo(test_fasta_dna)
            bio.count_characters(1,2,3,4)
        self.assertEqual(type(cm.exception), TypeError)


# Unittests for Bioinfo.change_characters()
# Changes characters in a string

    def test_change_characters_dict_as_translate_dict(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.change_characters(test_string, test_dict))

    def test_change_characters_str_as_translate_dict(self):
            with self.assertRaises(SystemExit) as cm:
                bio = Bioinfo(test_fasta_dna)
                bio.change_characters(test_string, test_string)
            self.assertEqual(cm.exception.code, 1)

    def test_change_characters_list_as_translate_dict(self):
            with self.assertRaises(SystemExit) as cm:
                bio = Bioinfo(test_fasta_dna)
                bio.change_characters(test_string, test_list)
            self.assertEqual(cm.exception.code, 1)

    def test_change_characters_dict_as_template(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_dna)
            bio.change_characters(test_dict, test_dict)
        self.assertEqual(cm.exception.code, 1)

    def test_change_characters_string_as_template(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.change_characters(test_string, test_dict))

    def test_change_characters_list_as_template(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.change_characters(test_list, test_dict))

    def test_change_characters_tuple_as_template(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.change_characters(test_tuple, test_dict))


# Unittest for Bioinfo.splice_string_by_character()

    def test_splice_string_by_charater_invalid_charater(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_empty)
            bio.splice_string_by_charaters("A", "")
        self.assertEqual(cm.exception.code, 1)

    def test_splice_string_by_charater_invalid_charater2(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_empty)
            bio.splice_string_by_charaters("", "A")
        self.assertEqual(cm.exception.code, 1)

    def test_splice_string_by_charater(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertIsInstance(bio.splice_string_by_charaters("A", "T"), list)


# Unittests for Bioinfo.splice_string_by_position()

    def test_splice_string_by_position_slice_to_big(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_empty)
            bio.splice_string_by_position(1)
        self.assertEqual(cm.exception.code, 1)

    def test_splice_string_by_position_step_to_big(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo(test_fasta_empty)
            bio.splice_string_by_position(0,2)
        self.assertEqual(cm.exception.code, 1)

    def test_splice_string_by_position(self):
        bio = Bioinfo(test_fasta_dna)
        self.assertTrue(bio.splice_string_by_position(3), list)


# Unittests for Bionfo_DNA

    def test_Bioinfo_DNA_Input_DNA(self):
        self.assertTrue(Bioinfo_DNA(test_fasta_dna))
        
    def test_Bioinfo_DNA_Input_Protein(self):
        with self.assertRaises(SystemExit) as cm:
            bio = Bioinfo_DNA(test_fasta_protein)
        self.assertEqual(cm.exception.code, 1)

    def test_Bioinfo_DNA_return_DNA(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.return_DNA(), str)

    def test_Bioinfo_DNA_return_RNA(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.return_RNA(), str)

    def test_Bioinfo_DNA_reverse_sequence(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.reverse_sequence(), str)

    def test_Bioinfo_DNA_complementary_sequence(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.complementary_sequence(), str)

    def test_Bioinfo_DNA_nucleotide_content(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.nucleotide_content(), str)

    def test_Bioinfo_DNA_GC_content(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.GC_content(), str)

    def test_Bioinfo_DNA_count_motif(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.count_motif("CAG"), dict)

    def test_Bioinfo_DNA_codons(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.codons(), list)

    def test_Bioinfo_DNA_translation(self):
        bio = Bioinfo_DNA(test_fasta_dna)
        self.assertIsInstance(bio.translation(bio.codons()), str)


# Unittests for Bionfo_Protein

    def test_Bioinfo_Protein_Input_DNA(self):
        self.assertTrue(Bioinfo_Protein(test_fasta_dna))

    def test_Bioinfo_Protein_Input_Protein(self):
        self.assertTrue(Bioinfo_Protein(test_fasta_protein))

    def test_Bioinfo_Protein_return_sequence(self):
        bio = Bioinfo_Protein(test_fasta_protein)
        self.assertIsInstance(bio.return_AA(), str)

    def test_Bioinfo_Protein_find_gene(self):
        bio = Bioinfo_Protein(test_fasta_dna)
        self.assertIsInstance(bio.find_gene(), list)


if __name__=='__main__': 
    unittest.main()
