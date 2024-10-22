# >>>>>>>>> CSV/JSON Functions <<<<<<<<<<<<
# Function to read file as json and return the data as a dictionary
# If file is not present, return None
def read_json_file(file_path):
    import json
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        return None
    
# Function that extracts path from a file path
def extract_path_filename(file_path):
    import os
    path_file_split = os.path.split(file_path)
    return path_file_split[0], path_file_split[1]

# Function that checks if path exists and creates it if it doesn't
def check_create_path(path):
    import os
    if not os.path.exists(path):
        print(f"(check_create_path) Creating directory: {path}")
        os.makedirs(path)

# Function to write data as json to a file
def write_json_file(file_path, data):
    import json
    path, filename = extract_path_filename(file_path)
    check_create_path(path)
    with open(file_path, 'w') as f:
        json.dump(data, f)

# Function to write data as csv to a file
def read_csv_file(file_path):
    import csv
    encoding = detect_text_file_encoding(file_path)
    if encoding is None:
        raise ValueError("Failed to detect the encoding of the csv file.")
    try:
        with open(file_path, 'r', encoding=encoding) as file:
            reader = csv.reader(file)
            data = [row for row in reader]
        return data
    except FileNotFoundError:
        return None

def detect_text_file_encoding(file_path):
    #import csv
    # List of encodings to try, in order of likelihood
    encodings_to_try = ['utf-8', 'cp1252', 'ISO-8859-1', 'latin1']
    for encoding in encodings_to_try:
        try:
            with open(file_path, mode='r', encoding=encoding) as file:
                # Read first line of text file
                first_line = file.readline()
                return encoding
        except UnicodeDecodeError:
            continue  # Try the next encoding if decoding failed
    else:
        print(f"Failed to open the file with any of the tried encodings ({encodings_to_try}).")
        return None