import os
import pandas as pd
import pickle
from multiprocessing import Pool, cpu_count
import time
from tqdm.auto import tqdm
import logging
from typing import List, Dict, Union
from unipressed.id_mapping.types import From, To
from unipressed import IdMappingClient
import itertools
import sys
from IPython.display import display, HTML

class IDTranslator:
    def __init__(self, input_file: str, output_file: str, source_type: From, dest_type: To,
                 pickle_dir: str = './pickles', batch_size: int = 100, processes: int = None,
                 columns: List[str] = None, max_retries: int = 3, clear_progress: bool = False,
                 input_columns: Union[List[str], Dict[str, str]] = None, has_header: bool = True,
                 multiple_mapping_strategy: str = 'expand'):
        self.input_file = input_file
        self.output_file = output_file
        self.source_type = source_type
        self.dest_type = dest_type
        self.pickle_dir = pickle_dir
        self.batch_size = batch_size
        self.processes = processes or max(1, cpu_count() - 1)
        self.columns = columns or ['source', 'target']
        self.max_retries = max_retries
        self.input_columns = input_columns
        self.has_header = has_header
        self.df = None
        self.id_mapping = {}
        self.unique_ids = set()
        self.multiple_mapping_strategy = multiple_mapping_strategy

        # Set up logging
        self.logger = logging.getLogger(f"{self.__class__.__name__}_{id(self)}")
        self.logger.setLevel(logging.INFO)

        # Create formatter
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # Create file handler
        fh = logging.FileHandler('id_translator.log')
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)

        # Create stream handler for notebook display
        class NotebookHandler(logging.Handler):
            def emit(self, record):
                msg = self.format(record)
                display(HTML(f"<pre>{msg}</pre>"))

        nh = NotebookHandler()
        nh.setLevel(logging.INFO)
        nh.setFormatter(formatter)
        self.logger.addHandler(nh)

        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

        os.makedirs(self.pickle_dir, exist_ok=True)

        if clear_progress:
            self.clear_all_progress()

    def load_data(self):
        file_extension = os.path.splitext(self.input_file)[1].lower()

        if file_extension == '.csv':
            self.df = pd.read_csv(self.input_file, header=0 if self.has_header else None)
        elif file_extension in ['.xls', '.xlsx']:
            self.df = pd.read_excel(self.input_file, header=0 if self.has_header else None)
        elif file_extension == '.tsv':
            self.df = pd.read_csv(self.input_file, sep='\t', header=0 if self.has_header else None)
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

        if not self.has_header:
            self.df.columns = [f'col_{i}' for i in range(len(self.df.columns))]

        self.process_columns()

        self.unique_ids = set()
        for column in self.columns:
            self.unique_ids.update(self.df[column].unique())
        self.logger.info(f"Loaded {len(self.df)} rows with {len(self.unique_ids)} unique IDs")

    def process_columns(self):
        if isinstance(self.input_columns, dict):
            self.df.rename(columns=self.input_columns, inplace=True)
            self.columns = list(self.input_columns.values())
        elif isinstance(self.input_columns, list):
            if len(self.input_columns) > len(self.df.columns):
                raise ValueError("More column names provided than exist in the file")
            self.df.columns = self.input_columns + list(self.df.columns[len(self.input_columns):])
            self.columns = self.input_columns[:2]  # Assume the first two columns are of interest
        elif self.input_columns is None:
            if self.columns is None:
                self.columns = list(self.df.columns[:2])  # Assume the first two columns are of interest
        else:
            raise ValueError("input_columns must be a list, dict, or None")

    def load_progress(self):
        pickle_files = [f for f in os.listdir(self.pickle_dir) if f.endswith('.pkl')]
        if pickle_files:
            latest_pickle = max(pickle_files, key=lambda x: os.path.getmtime(os.path.join(self.pickle_dir, x)))
            with open(os.path.join(self.pickle_dir, latest_pickle), 'rb') as f:
                self.id_mapping = pickle.load(f)
            self.logger.info(f"Loaded progress from {latest_pickle}")
        else:
            self.logger.info("No existing progress found. Starting from scratch.")

    def clear_all_progress(self):
        for file in os.listdir(self.pickle_dir):
            if file.endswith('.pkl'):
                os.remove(os.path.join(self.pickle_dir, file))
        self.logger.info("Cleared all progress files.")

    @staticmethod
    def get_ids(args):
        batch, source_type, dest_type, max_retries = args
        results = {}
        for id_ in batch:
            for attempt in range(max_retries):
                try:
                    request = IdMappingClient.submit(source=source_type, dest=dest_type, ids={id_})
                    time.sleep(1)  # Sleep for 1 second
                    result = list(request.each_result())
                    results[id_] = [m['to'] for m in result] if result else []
                    break
                except Exception as e:
                    if attempt == max_retries - 1:
                        logging.error(f"Error processing {id_} after {max_retries} attempts: {str(e)}")
                        results[id_] = []
                    time.sleep(min(2 ** attempt, 10))  # Cap the backoff at 10 seconds
        return results

    def translate(self):
        self.load_data()
        self.load_progress()

        remaining_ids = self.unique_ids - set(self.id_mapping.keys())
        id_batches = [list(remaining_ids)[i:i + self.batch_size] for i in range(0, len(remaining_ids), self.batch_size)]

        with Pool(self.processes) as pool:
            progress_bar = tqdm(total=len(id_batches), desc="Processing batches")
            for i, batch_results in enumerate(pool.imap_unordered(self.get_ids,
                                                                  [(batch, self.source_type, self.dest_type,
                                                                    self.max_retries)
                                                                   for batch in id_batches])):
                self.id_mapping.update(batch_results)
                progress_bar.update(1)
                self.logger.info(f"Processed batch {i + 1}/{len(id_batches)}")

                if (i + 1) % 2 == 0 or i == len(id_batches) - 1:
                    self.save_progress(f'checkpoint_batch_{i + 1}.pkl')

            progress_bar.close()

        self.apply_translation()

    def save_progress(self, filename: str):
        with open(os.path.join(self.pickle_dir, filename), 'wb') as f:
            pickle.dump(self.id_mapping, f)
        self.logger.info(f"Progress saved to {filename}")

    def apply_translation(self):
        self.logger.info("Applying translation to dataframe...")

        if self.multiple_mapping_strategy == 'expand':
            expanded_rows = []
            for _, row in self.df.iterrows():
                source_mappings = self.id_mapping.get(row[self.columns[0]])
                target_mappings = self.id_mapping.get(row[self.columns[1]])

                # If either mapping is None, skip this row
                if source_mappings is None or target_mappings is None:
                    continue

                # If mappings are empty lists, use the original identifiers
                source_mappings = source_mappings if source_mappings else [row[self.columns[0]]]
                target_mappings = target_mappings if target_mappings else [row[self.columns[1]]]

                for source, target in itertools.product(source_mappings, target_mappings):
                    new_row = row.copy()
                    new_row[f"{self.columns[0]}_{self.dest_type}"] = source
                    new_row[f"{self.columns[1]}_{self.dest_type}"] = target
                    expanded_rows.append(new_row)

            if expanded_rows:
                self.df = pd.DataFrame(expanded_rows)
            else:
                logging.warning("No valid translations found. The DataFrame is empty.")
        else:
            for column in self.columns:
                new_column = f"{column}_{self.dest_type}"
                self.df[new_column] = self.df[column].map(lambda x: self.id_mapping.get(x, [x]))
                self.df[new_column] = self.df[new_column].apply(lambda x: ';'.join(x) if isinstance(x, list) else x)

        # Save the translated dataframe
        self._save_dataframe(self.output_file)

    def run(self):
        self.logger.info(f"Starting ID translation process from {self.source_type} to {self.dest_type}")
        start_time = time.time()
        self.translate()
        end_time = time.time()
        self.logger.info(f"ID translation process completed in {end_time - start_time:.2f} seconds")
        self.analyze_results()

    def analyze_results(self):
        original_count = len(self.df)
        translated_count = self.df[[f"{col}_{self.dest_type}" for col in self.columns]].notna().all(axis=1).sum()
        expansion_factor = translated_count / original_count if original_count > 0 else 0

        self.logger.info(f"Original entry count: {original_count}")
        self.logger.info(f"Translated entry count: {translated_count}")
        self.logger.info(f"Expansion factor: {expansion_factor:.2f}")

        success_rate = (translated_count / original_count) * 100 if original_count > 0 else 0
        self.logger.info(f"Translation success rate: {success_rate:.2f}%")

    def remove_untranslated_entries(self, output_file: str = None):
        if self.df is None:
            raise ValueError("No data loaded. Please run the translation process first.")

        original_count = len(self.df)

        # Identify columns that contain translated IDs
        translated_columns = [f"{col}_{self.dest_type}" for col in self.columns]

        # Remove rows where any of the translated columns is null
        self.df = self.df.dropna(subset=translated_columns)

        removed_count = original_count - len(self.df)

        if output_file is None:
            file_name, file_extension = os.path.splitext(self.output_file)
            output_file = f"{file_name}_cleaned{file_extension}"

        # Save the cleaned dataframe
        self._save_dataframe(output_file)

        self.logger.info(f"Removed {removed_count} untranslated entries.")
        self.logger.info(f"Cleaned database saved to {output_file}")
        self.logger.info(f"Original entry count: {original_count}")
        self.logger.info(f"Cleaned entry count: {len(self.df)}")

        # Update the success rate
        self.analyze_results()

    def _save_dataframe(self, output_file):
        if self.df.empty:
            logging.warning("The DataFrame is empty. No file will be saved.")
            return

        file_extension = os.path.splitext(output_file)[1].lower()
        if file_extension == '.csv':
            self.df.to_csv(output_file, index=False)
        elif file_extension in ['.xls', '.xlsx']:
            self.df.to_excel(output_file, index=False)
        elif file_extension == '.tsv':
            self.df.to_csv(output_file, sep='\t', index=False)
        else:
            raise ValueError(f"Unsupported output file format: {file_extension}")

        self.logger.info(f"Results saved to {output_file}")
