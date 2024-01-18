import sys
import pandas as pd

# Written by Isabel Díez-Santos, Sathish Thiyagarajan, and ChatGPT.

def window_mod_data(threshold, window_size, input_file, output_file):
    def convert_threshold(value):
        return 1 if value >= threshold else 0
    """
    Applies windowing on raw modification data after thresholding the data.

    Args:
    - input_file: TSV file containing comma-separated read_id, start, end, mod_qual, and detect columns.
    - threshold (float): Threshold for calling modifications.
    - window_size (int): Size in bases for data windowing.

    Returns:
    - None.
    """

    # Read the data from the input file into a pandas DataFrame
    data = pd.read_csv(input_file, sep="\t")

    # Apply conversion based on the threshold
    data['mod_qual'] = data['mod_qual'].apply(convert_threshold)

    # Calculate averages within the given window size along with start and end values
    windows = []
    start_values = []
    end_values = []
    read_ids = []
    detect_values = []

    for i in range(0, len(data), window_size):
        window = data.iloc[i:i + window_size]
        start = window['ref_position'].iloc[0] if 'ref_position' in data.columns else window['forward_read_position'].iloc[0]
        end = window['end'].iloc[-1]
        start_values.append(start)
        end_values.append(end)

        read_id = window['read_id'].iloc[0]
        read_ids.append(read_id)

        detect_values.append('winVal')

        average = window['mod_qual'].mean()
        windows.append(average)

    # Create a new DataFrame with the windowed averages, start, end, read_id, and label values
    windowed_data = pd.DataFrame({
        'read_id': read_ids,
        'start': start_values,
        'end': end_values,
        'mod_qual': windows,
        'label': detect_values
    })

    # if there are multiple read ids, throw error
    if windowed_data["read_id"].value_counts() != 1:
        raise ValueError("Input must contain one unique read id")

    # Write the windowed data to the output file
    windowed_data.to_csv(output_file, index=False, sep="\t")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python window_mod_data.py <threshold> <window_size> <input_file> <output_file>")
        sys.exit(1)

   # Check if threshold is provided and within the valid range (0 < threshold < 1)
    try:
        threshold = float(sys.argv[1])
        if not (0 < threshold < 1):
            print("Threshold value must be > 0 and < 1")
            sys.exit(1)
    except ValueError:
        print("Threshold value not provided")
        sys.exit(1)

    # Check if window_size is provided
    try:
        window_size = int(sys.argv[2])
    except ValueError:
        print("Window_size value not provided")
        sys.exit(1)


    input_file = sys.argv[3]
    output_file = sys.argv[4]

    window_mod_data(threshold, window_size, input_file, output_file)
