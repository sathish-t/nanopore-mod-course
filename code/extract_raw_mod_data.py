import pandas as pd
import sys
import os

def extract_raw_mod_data(input_file, mod_code, position):
    """
    Extracts required columns from the input TSV file based on specified conditions.

    Args:
    - input_file (str): Path to the input TSV file.
    - mod_code (str): Modification code used to filter the data (e.g., m or h).
    - position (str): Indicates the position information to use ('not_use_ref' or 'use_ref').

    Returns:
    - tsv file: Extracted data containing read_id, start, end, mod_qual, and label columns.
    """
    data = pd.read_csv(input_file, sep='\t')

    if mod_code not in data['mod_code'].unique():
        print(f"Modification '{mod_code}' is not present in the data.")
        sys.exit(1)

    filtered_data = data[data['mod_code'] == mod_code].copy()

    if position == 'not_use_ref':
        start_col = 'forward_read_position'
    elif position == 'use_ref':
        start_col = 'ref_position'
    else:
        raise ValueError("Invalid position provided. Position must be: not_use_ref or use_ref")

    output_data = filtered_data[['read_id', start_col, 'mod_qual']].copy()
    output_data['end'] = filtered_data[start_col] + 1
    output_data['label'] = 'rawVal'
    output_data = output_data[['read_id', start_col, 'end', 'mod_qual', 'label']]
    return output_data

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_raw_mod_data.py <mod_code> <position> <input_file> <output_file>")
        sys.exit(1)

    mod_code = sys.argv[1]
    position = sys.argv[2]
    input_file = sys.argv[3]
    output_file = sys.argv[4]

    output_data = extract_raw_mod_data(input_file, mod_code, position)

    output_data.to_csv(output_file, index=False)
