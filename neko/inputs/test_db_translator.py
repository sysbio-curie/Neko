import unittest
import pandas as pd
import os
import tempfile
from unittest.mock import patch, MagicMock
from unipressed.id_mapping.types import From, To
from unipressed import IdMappingClient


class TestIDTranslator(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.input_file = os.path.join(self.temp_dir, 'test_input.csv')
        self.output_file = os.path.join(self.temp_dir, 'test_output.csv')

        # Create a sample input file
        df = pd.DataFrame({
            'source': ['ENSG00000139618', 'ENSG00000141510', 'ENSG00000157764'],
            'target': ['ENSG00000155657', 'ENSG00000121879', 'ENSG00000171862']
        })
        df.to_csv(self.input_file, index=False)

        self.translator = IDTranslator(
            input_file=self.input_file,
            output_file=self.output_file,
            source_type=From.ENSEMBL,
            dest_type=To.UNIPROTKB,
            pickle_dir=self.temp_dir
        )

    def tearDown(self):
        # Clean up temporary files
        for file in os.listdir(self.temp_dir):
            os.remove(os.path.join(self.temp_dir, file))
        os.rmdir(self.temp_dir)

    def test_load_data(self):
        self.translator.load_data()
        self.assertIsNotNone(self.translator.df)
        self.assertEqual(len(self.translator.df), 3)
        self.assertEqual(len(self.translator.unique_ids), 6)

    @patch('id_translator.IdMappingClient.submit')
    def test_get_id(self, mock_submit):
        mock_submit.return_value = MagicMock()
        mock_submit.return_value.each_result.return_value = [{'to': 'P12345'}]

        result = self.translator.get_id(('ENSG00000139618', From.ENSEMBL, To.UNIPROTKB, 3))
        self.assertEqual(result, ('ENSG00000139618', 'P12345'))

    @patch('id_translator.IdMappingClient.submit')
    def test_translate(self, mock_submit):
        mock_submit.return_value = MagicMock()
        mock_submit.return_value.each_result.return_value = [{'to': 'P12345'}]

        self.translator.translate()
        self.assertTrue(os.path.exists(self.output_file))

        df = pd.read_csv(self.output_file)
        self.assertEqual(len(df), 3)
        self.assertTrue('source_UniProtKB' in df.columns)
        self.assertTrue('target_UniProtKB' in df.columns)

    def test_analyze_results(self):
        # Create a sample translated DataFrame
        df = pd.DataFrame({
            'source': ['ENSG1', 'ENSG2', 'ENSG3'],
            'target': ['ENSG4', 'ENSG5', 'ENSG6'],
            'source_UniProtKB': ['P1', 'P2', None],
            'target_UniProtKB': ['P4', None, 'P6']
        })
        self.translator.df = df

        with self.assertLogs(level='INFO') as cm:
            self.translator.analyze_results()

        self.assertIn("Translation success rate: 66.67%", cm.output[0])
        self.assertIn("Untranslated entries: 2 out of 6", cm.output[1])

    def test_remove_untranslated_entries(self):
        # Create a sample translated DataFrame
        df = pd.DataFrame({
            'source': ['ENSG1', 'ENSG2', 'ENSG3'],
            'target': ['ENSG4', 'ENSG5', 'ENSG6'],
            'source_UniProtKB': ['P1', 'P2', None],
            'target_UniProtKB': ['P4', None, 'P6']
        })
        self.translator.df = df

        cleaned_output = os.path.join(self.temp_dir, 'cleaned_output.csv')
        self.translator.remove_untranslated_entries(cleaned_output)

        cleaned_df = pd.read_csv(cleaned_output)
        self.assertEqual(len(cleaned_df), 1)  # Only one row should remain
        self.assertEqual(cleaned_df.iloc[0]['source'], 'ENSG1')
        self.assertEqual(cleaned_df.iloc[0]['target'], 'ENSG4')


if __name__ == '__main__':
    unittest.main()
