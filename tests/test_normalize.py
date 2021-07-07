from qebil.normalize import (
    qiimp_parser,
    apply_validation,
    normalize_lat_lon,
    add_emp_info,
    split_lat_lon
)
import unittest
from pandas.testing import assert_frame_equal
from os import path
import pandas as pd

from qebil.tools.util import setup_output_dir

_THIS_DIR, _THIS_FILENAME = path.split(__file__)

_TEST_SUPPORT_DIR = path.join(_THIS_DIR, "support_files")

_TEST_OUTPUT_DIR = path.join(_THIS_DIR, "test_output/")

setup_output_dir(_TEST_OUTPUT_DIR)


class NormalizeTest(unittest.TestCase):
    
    def setUp(self):
        self.maxDiff = None
        
    def test_qiimp_parser(self):
        """Tests parsing of Qiimp-format yaml or xlsx files to obtain a dictionary of
        rules for metadata normalization
        """

        expected_yml = {
            "sample_name": {
                "empty": False,
                "is_phi": False,
                "regex": "^[a-zA-Z0-9\\.]+$",
                "required": True,
                "type": "string",
                "unique": True,
            }
        }
        failed_yml = {}
        test_yml = _TEST_SUPPORT_DIR + "/test.yml"
        test_xlsx = _TEST_SUPPORT_DIR + "/test.xlsx"

        self.assertEqual(expected_yml, qiimp_parser(test_yml))
        self.assertEqual(expected_yml, qiimp_parser(test_xlsx))
        self.assertEqual(failed_yml, qiimp_parser("not a file"))

    def test_apply_validation(self):

        test_yml = {
            "collection_device": {
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "collection_method": {
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "collection_timestamp": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {
                        "regex": "^([0-9]{4})(?:-([0-1][0-2])(?:-([0-3][0-9])(?: ([0-2][0-9])(?::([0-5][0-9])(?::([0-5][0-9]))?)?)?)?)?$",
                        "type": "datetime",
                    },
                ],
                "empty": False,
                "field_desc": "The day and time of sampling as a single point in time expressed in 24-hour time format, e.g. 2016-11-22.",
                "is_phi": False,
                "required": True,
            },
            "description": {
                "empty": False,
                "field_desc": "A description of the sample that can include site, subject, and sample material.",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "dna_extracted": {
                "anyof": [
                    {
                        "allowed": ["not collected", "not provided"],
                        "type": "string",
                    },
                    {"allowed": ["TRUE", "FALSE"], "type": "string"},
                ],
                "empty": False,
                "field_desc": "Whether the DNA been extracted from the sample.",
                "is_phi": False,
                "required": True,
            },
            "elevation": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {"min": "-413.0", "type": "number"},
                ],
                "empty": False,
                "field_desc": "Height of land above sea level in meters at the sampling site",
                "is_phi": False,
                "required": True,
                "units": "elevation_units",
            },
            "elevation_units": {
                "allowed": ["meters"],
                "default": "meters",
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "empo_1": {
                "allowed": ["Host-associated", "Free-living"],
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "empo_2": {
                "allowed": [
                    "Animal",
                    "Fungus",
                    "Plant",
                    "Saline",
                    "Non-saline",
                ],
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "empo_3": {
                "allowed": [
                    "Plant secretion",
                    "Plant rhizosphere",
                    "Plant surface",
                    "Plant corpus",
                    "Fungus surface",
                    "Fungus corpus",
                    "Animal surface",
                    "Animal proximal gut",
                    "Animal corpus",
                    "Animal distal gut",
                    "Animal secretion",
                    "Sub-surface (non-saline)",
                    "Surface (non-saline)",
                    "Water (non-saline)",
                    "Soil (non-saline)",
                    "Aerosol (non-saline)",
                    "Sediment (non-saline)",
                    "Surface (saline)",
                    "Water (saline)",
                    "Hypersaline (saline)",
                    "Aerosol (saline)",
                    "Sediment (saline)",
                ],
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "env_biome": {
                "empty": False,
                "field_desc": "Classification of the location where the sample was obtained, from the Environmental Ontology (ENVO).",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "env_feature": {
                "empty": False,
                "field_desc": "Classification of a specific feature in the biome, from the Environmental Ontology (ENVO).",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "env_material": {
                "empty": False,
                "field_desc": "Classification of the material being sampled, from the Environmental Ontology (ENVO).",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "env_package": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {
                        "allowed": [
                            "air",
                            "built environment",
                            "host-associated",
                            "human-associated",
                            "human-skin",
                            "human-oral",
                            "human-gut",
                            "human-vaginal",
                            "microbial mat/biofilm",
                            "misc environment",
                            "plant-associated",
                            "sediment",
                            "soil",
                            "wastewater/sludge",
                            "water",
                        ],
                        "type": "string",
                    },
                ],
                "empty": False,
                "field_desc": "Environment where the sample was obtained.",
                "is_phi": False,
                "required": True,
            },
            "geo_loc_name": {
                "empty": False,
                "field_desc": "The geographical origin of the sample as defined by the country or location as chosen from the GAZ ontology, e.g. USA:CA:San Diego .",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "host_subject_id": {
                "empty": False,
                "field_desc": "A label that applies to all samples belonging to one subject or one sample",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "latitude": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {"max": "90.00", "min": "-90.00", "type": "number"},
                ],
                "empty": False,
                "field_desc": "Location of the sample site by latitude in decimal degrees, e.g. 32.72 for San Diego, USA.",
                "is_phi": False,
                "required": True,
                "units": "latitude_units",
            },
            "latitude_units": {
                "allowed": ["Decimal degrees"],
                "default": "Decimal degrees",
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "longitude": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {"max": "180.00", "min": "-180.00", "type": "number"},
                ],
                "empty": False,
                "field_desc": "Location of the sample site by longitude in decimal degrees, e.g. -117.16 for San Diego, USA.",
                "is_phi": False,
                "required": True,
                "units": "longitude_units",
            },
            "longitude_units": {
                "allowed": ["Decimal degrees"],
                "default": "Decimal degrees",
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "physical_specimen_location": {
                "empty": False,
                "field_desc": "Where the sample is stored (if there is any left).",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "physical_specimen_remaining": {
                "anyof": [
                    {
                        "allowed": [
                            "not collected",
                            "not provided",
                            "restricted access",
                        ],
                        "type": "string",
                    },
                    {"allowed": ["TRUE", "FALSE"], "type": "string"},
                ],
                "empty": False,
                "field_desc": "Whether there is any raw specimen left to sample.",
                "is_phi": False,
                "required": True,
            },
            "sample_name": {
                "empty": False,
                "is_phi": False,
                "regex": "^[a-zA-Z0-9\\.]+$",
                "required": True,
                "type": "string",
                "unique": True,
            },
            "sample_type": {
                "allowed": ["other"],
                "default": "other",
                "empty": False,
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "scientific_name": {
                "allowed": ["metagenome"],
                "default": "metagenome",
                "empty": False,
                "field_desc": "The scientific name of the metagenome associated with this sample, as defined in the NCBI Taxonomy. Be aware that this is NOT the scientific name of the host! Note that this value must correspond to the taxon id entered in the taxon_id field; for example, for a sample from the human gut, the scientific_name should be human gut microbiome, corresponding to taxon id 408170.",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "taxon_id": {
                "allowed": ["256318"],
                "default": "256318",
                "empty": False,
                "field_desc": "The NCBI Taxonomy id of the metagenome associated with this sample. Be aware that this is NOT the taxon id of the host! Note that this value must correspond to the scientific name entered in the scientific_name field; for example, for a sample from the human gut, the taxon_id should be 408170, corresponding to the scientific name human gut metagenome.",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "title": {
                "empty": False,
                "field_desc": "Title of study.",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
            "tube_id": {
                "empty": False,
                "field_desc": "The identification of the tube containing the sample.",
                "is_phi": False,
                "required": True,
                "type": "string",
            },
        }
        expected_df = pd.DataFrame.from_dict(
            {
                "qebil_collection_device": ["not provided","not provided"],
                "qebil_collection_method": ["not provided","not provided"],
                "qebil_collection_timestamp": ["not provided","not provided"],
                "qebil_description": ["not provided","not provided"],
                "qebil_dna_extracted": ["not provided","not provided"],
                "qebil_elevation": ["not provided","not provided"],
                "qebil_elevation_units": ["meters","meters"],
                "qebil_empo_1": ["not provided","not provided"],
                "qebil_empo_2": ["not provided","not provided"],
                "qebil_empo_3": ["not provided","not provided"],
                "qebil_env_biome": ["not provided","not provided"],
                "qebil_env_feature": ["not provided","not provided"],
                "qebil_env_material": ["not provided","not provided"],
                "qebil_env_package": ["not provided","not provided"],
                "qebil_geo_loc_name": ["not provided","not provided"],
                "qebil_host_subject_id": ["not provided","not provided"],
                "qebil_latitude": ["not provided","not provided"],
                "qebil_latitude_units": ["Decimal degrees","Decimal degrees"],
                "qebil_longitude": ["not provided","not provided"],
                "qebil_longitude_units": ["Decimal degrees","Decimal degrees"],
                "qebil_physical_specimen_location": ["not provided","not provided"],
                "qebil_physical_specimen_remaining": ["not provided","not provided"],
                "prep_file": ["AMBIGUOUS_0", "Metagenomic_0"],
                "sample_name": ["sample1", "sample2"],
                "qebil_sample_type": ["other","other"],
                "qebil_scientific_name": ["metagenome","metagenome"],
                "qebil_taxon_id": ["256318","256318"],
                "qebil_title": ["not provided","not provided"],
                "qebil_tube_id": ["not provided","not provided"],
            }
        )
        test_df = pd.DataFrame.from_dict(
            {
                "sample_name": ["sample1", "sample2"],
                "prep_file": ["AMBIGUOUS_0", "Metagenomic_0"],
            }
        )
        expected_msg = (
            "collection_device not found in metadata. Attempting to generate qebil_collection_device from available sources.\n"
            + "collection_device has no default in yaml template. Encoding as 'not provided'\n"
            + "collection_method not found in metadata. Attempting to generate qebil_collection_method from available sources.\n"
            + "collection_method has no default in yaml template. Encoding as 'not provided'\n"
            + "collection_timestamp not found in metadata. Attempting to generate qebil_collection_timestamp from available sources.\n"
            + "collection_timestamp has no default in yaml template. Encoding as 'not provided'\n"
            + "description not found in metadata. Attempting to generate qebil_description from available sources.\n"
            + "description has no default in yaml template. Encoding as 'not provided'\n"
            + "dna_extracted not found in metadata. Attempting to generate qebil_dna_extracted from available sources.\n"
            + "dna_extracted has no default in yaml template. Encoding as 'not provided'\n"
            + "elevation not found in metadata. Attempting to generate qebil_elevation from available sources.\n"
            + "elevation has no default in yaml template. Encoding as 'not provided'\n"
            + "elevation_units not found in metadata. Attempting to generate qebil_elevation_units from available sources.\n"
            + "Setting qebil_elevation_units to meters\n"
            + "empo_1 not found in metadata. Attempting to generate qebil_empo_1 from available sources.\n"
            + "empo_1 has no default in yaml template. Encoding as 'not provided'\n"
            + "empo_2 not found in metadata. Attempting to generate qebil_empo_2 from available sources.\n"
            + "empo_2 has no default in yaml template. Encoding as 'not provided'\n"
            + "empo_3 not found in metadata. Attempting to generate qebil_empo_3 from available sources.\n"
            + "empo_3 has no default in yaml template. Encoding as 'not provided'\n"
            + "env_biome not found in metadata. Attempting to generate qebil_env_biome from available sources.\n"
            + "env_biome has no default in yaml template. Encoding as 'not provided'\n"
            + "env_feature not found in metadata. Attempting to generate qebil_env_feature from available sources.\n"
            + "env_feature has no default in yaml template. Encoding as 'not provided'\n"
            + "env_material not found in metadata. Attempting to generate qebil_env_material from available sources.\n"
            + "env_material has no default in yaml template. Encoding as 'not provided'\n"
            + "env_package not found in metadata. Attempting to generate qebil_env_package from available sources.\n"
            + "env_package has no default in yaml template. Encoding as 'not provided'\n"
            + "geo_loc_name not found in metadata. Attempting to generate qebil_geo_loc_name from available sources.\n"
            + "geo_loc_name has no default in yaml template. Encoding as 'not provided'\n"
            + "host_subject_id not found in metadata. Attempting to generate qebil_host_subject_id from available sources.\n"
            + "host_subject_id has no default in yaml template. Encoding as 'not provided'\n"
            + "latitude not found in metadata. Attempting to generate qebil_latitude from available sources.\n"
            + "latitude has no default in yaml template. Encoding as 'not provided'\n"
            + "latitude_units not found in metadata. Attempting to generate qebil_latitude_units from available sources.\n"
            + "Setting qebil_latitude_units to Decimal degrees\n"
            + "longitude not found in metadata. Attempting to generate qebil_longitude from available sources.\n"
            + "longitude has no default in yaml template. Encoding as 'not provided'\n"
            + "longitude_units not found in metadata. Attempting to generate qebil_longitude_units from available sources.\n"
            + "Setting qebil_longitude_units to Decimal degrees\n"
            + "physical_specimen_location not found in metadata. Attempting to generate qebil_physical_specimen_location from available sources.\n"
            + "physical_specimen_location has no default in yaml template. Encoding as 'not provided'\n"
            + "physical_specimen_remaining not found in metadata. Attempting to generate qebil_physical_specimen_remaining from available sources.\n"
            + "physical_specimen_remaining has no default in yaml template. Encoding as 'not provided'\n"
            + "sample_type not found in metadata. Attempting to generate qebil_sample_type from available sources.\n"
            + "Setting qebil_sample_type to other\n"
            + "scientific_name not found in metadata. Attempting to generate qebil_scientific_name from available sources.\n"
            + "Setting qebil_scientific_name to metagenome\n"
            + "taxon_id not found in metadata. Attempting to generate qebil_taxon_id from available sources.\n"
            + "Setting qebil_taxon_id to 256318\n"
            + "title not found in metadata. Attempting to generate qebil_title from available sources.\n"
            + "title has no default in yaml template. Encoding as 'not provided'\n"
            + "tube_id not found in metadata. Attempting to generate qebil_tube_id from available sources.\n"
            + "tube_id has no default in yaml template. Encoding as 'not provided'\n"
                       )
        result_df, result_msg = apply_validation(test_df, test_yml)

        self.assertEqual(result_msg, expected_msg)
        assert_frame_equal(expected_df.sort_index(axis=1), result_df.sort_index(axis=1))

        # Addtional potential tests:
        # need a column that will test if the key from the validator is not in the test_df
        # look for the same thing but with 'qebil_ at the start
        # the final normalized column should be prepended with qebil_
        # going to omit testing the output message for now due to complexity
        # but do need to have at least one column with a source and have it map.
        # and one column with a map and a mapping dictionary
        # need to test that anyof rule is applied for a string
        # need to test that any of rule is applied safely with a string and numeric mix
        # need to test that numeric rules apply for min, max, min_exclusive, and max_exclusive
        # need to test when not anyof, but rather type is applied that strings are checked
        # and that numbers are checked against min, max, min_exclusive, max_exclusive
        #

    def test_normalize_lat_lon(self):
        test_df = pd.DataFrame.from_dict(
            {"sample_name": ["sample1"], "lat_lon": ["30.01 N 90.01 W"]}
        )
        expected_df = pd.DataFrame.from_dict(
            {
                "sample_name": ["sample1"],
                "lat_lon": ["30.01 N 90.01 W"],
                "qebil_latitude": ["30.01"],
                "qebil_longitude": ["-90.01"],
            }
        )
        test_df_no_lat_lon = pd.DataFrame.from_dict(
            {
                "sample_name": ["sample1"],
                "latitude": ["30.01"],
                "longitude": ["-90.01"],
            }
        )
        assert_frame_equal(expected_df, normalize_lat_lon(test_df))
        assert_frame_equal(
            test_df_no_lat_lon, normalize_lat_lon(test_df_no_lat_lon)
        )

    def test_add_emp_info(self):
        """
        Tests method to automatically adds/replace information with the
        standard EMP protocol information for 16S rRNA gene sequencing.
        """
        test_df = pd.DataFrame.from_dict(
            {
                "library_strategy": ["AMPLICON", "WGS"],
                "sample_name": ["sample1", "sample2"],
                "qebil_prep_file": ["AMBIGUOUS_0", "Metagenomic_0"],
            }
        )
        expected_df = pd.DataFrame.from_dict(
            {
                "library_strategy": ["AMPLICON", "WGS"],
                "sample_name": ["sample1", "sample2"],
                "qebil_prep_file": ["16S_0", "Metagenomic_0"],
                "target_gene": ["16S rRNA", "not applicable"],
                "target_subfragment": ["V4", "not applicable"],
                "library_construction_protocol": [
                    "Illumina EMP protocol 515fbc, 806r amplification of 16SrRNA V4",
                    "not applicable",
                ],
                "pcr_primers": [
                    "FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT",
                    "not applicable",
                ],
                "primer": ["GTGTGCCAGCMGCCGCGGTAA", "not applicable"],
                "sequencing_meth": [
                    "Sequencing by synthesis",
                    "not applicable",
                ],
            }
        )
        assert_frame_equal(expected_df, add_emp_info(test_df))

    def test_split_lat_lon(self):
        """ """
        test_dict = {
            "30.01 N 90.01 W": {"lat": "30.01", "long": "-90.01"},
            "30.01 S 90.01 W": {"lat": "-30.01", "long": "-90.01"},
            "30.01 N 90.01 E": {"lat": "30.01", "long": "90.01"},
            "30.01 S 90.01 E": {"lat": "-30.01", "long": "90.01"},
        }
        test_non_standard = "30.01 -90.01"
        results_dict = {}
        for k in test_dict.keys():
            self.assertEqual(k, split_lat_lon(k))
            self.assertEqual(test_dict[k]["lat"], split_lat_lon(k, "lat"))
            self.assertEqual(
                test_dict[k]["long"], split_lat_lon(k, "long")
            )

        self.assertEqual(
            test_non_standard, split_lat_lon(test_non_standard, "lat")
        )
        self.assertEqual(
            test_non_standard, split_lat_lon(test_non_standard, "long")
        )


if __name__ == "__main__":
    # begin the unittest.main()
    unittest.main()
