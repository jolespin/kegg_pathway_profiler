#!/usr/bin/env python
import sys, os, time, gzip, bz2, subprocess, pickle, json, logging, hashlib
from datetime import datetime
from collections import defaultdict
import pathlib
from pandas.errors import EmptyDataError

# Read/Write
# ==========
# Get file object
def open_file_reader(filepath: str, compression="auto", binary=False):
    """
    Opens a file for reading with optional compression.

    Args:
        filepath (str): Path to the file.
        compression (str, optional): Type of compression {None, 'gzip', 'bz2'}. Defaults to "auto".
        binary (bool, optional): Whether to open the file in binary mode. Defaults to False.

    Returns:
        file object: A file-like object.
    """
    if filepath == sys.stdin:
        return sys.stdin
    else:
        # Determine compression type based on the file extension if 'auto' is specified
        if compression == "auto":
            ext = filepath.split(".")[-1].lower()
            if ext == "gz":
                compression = "gzip"
            elif ext == "bz2":
                compression = "bz2"
            else:
                compression = None

        # Determine the mode based on the 'binary' flag
        mode = "rb" if binary else "rt"

        # Open the file with or without compression
        if not compression:
            return open(filepath, mode)
        elif compression == "gzip":
            return gzip.open(filepath, mode)
        elif compression == "bz2":
            return bz2.open(filepath, mode)
        else:
            raise ValueError(f"Unsupported compression type: {compression}")
                
# Get file object
def open_file_writer(filepath: str, compression="auto", binary=False):
    """
    Args:
        filepath (str): path/to/file
        compression (str, optional): {None, gzip, bz2}. Defaults to "auto".
        binary (bool, optional): Whether to open the file in binary mode. Defaults to False.
    
    Returns:
        file object
    """
    if compression == "auto":
        ext = filepath.split(".")[-1].lower()
        if ext == "gz":
            compression = "gzip"
        elif ext == "bz2":
            compression = "bz2"
        else:
            compression = None

    if binary:
        mode = "wb"
    else:
        mode = "wt"

    if not compression:
        return open(filepath, mode)
    elif compression == "gzip":
        return gzip.open(filepath, mode)
    elif compression == "bz2":
        return bz2.open(filepath, mode)
    else:
        raise ValueError(f"Unsupported compression type: {compression}")

# Pickle I/O
def read_pickle(filepath, compression="auto"):
    with open_file_reader(filepath, compression=compression, binary=True) as f:
        return pickle.load(f)
    
def write_pickle(obj, filepath, compression="auto"):
    with open_file_writer(filepath, compression=compression, binary=True) as f:
        pickle.dump(obj, f)
        
# Json I/O
def read_json(filepath):
    with open_file_reader(filepath, compression=None, binary=False) as f:
        return json.load(f)
    
def write_json(obj, filepath, indent=4):
    with open_file_writer(filepath, compression=None, binary=False) as f:
        return json.dump(obj, f)
    
from collections import defaultdict

def read_kos(filepath: str, name: str):
    """
    Reads KEGG Orthology (KO) identifiers from a file and maps them to corresponding genomes.

    This function processes a file containing KO identifiers, either for a single genome or multiple genomes.
    The file can either have one KO per line or be a tab-separated table with two columns: 
    [id_genome] <tab> [id_ko]. It returns a dictionary mapping genome identifiers to sets of KO identifiers.

    Parameters
    ----------
    filepath : str
        The path to the input file containing KO identifiers.
    
    name : str
        The identifier for the genome if the file contains KOs for a single genome.

    Returns
    -------
    genome_to_kos : dict
        A dictionary where keys are genome identifiers and values are sets of KO identifiers.
    
    Raises
    ------
    IndexError
        If the input file has a number of columns that is not 1 or 2.
    
    ValueError
        If any KO identifier in the file is not 6 characters long or does not start with 'K'.
    
    Notes
    -----
    - The function assumes that the input file has no header.
    - The function checks the format of KO identifiers and raises an error if they do not conform to the expected format.
    """
    
    # Initialize a dictionary to store genome-to-KO mappings
    genome_to_kos = defaultdict(set)
    
    with open_file_reader(filepath) as f:
        # Read the first line and determine the number of columns
        first_line = next(f).strip()
        fields = first_line.split("\t")
        number_of_fields = len(fields)
        
        # Check for the correct number of columns
        if number_of_fields not in {1, 2}:
            raise IndexError(
                f"Input file must have either 1 KO per line or a tab-separated table with 2 columns: "
                f"[id_genome] <tab> [id_ko]. No header expected. Current table contains {number_of_fields} columns."
            )
        
        # Case 1: Single genome (1 column per line)
        if number_of_fields == 1:
            if name is None:
                raise ValueError("If --kos are a single column of KOs then --name must be provided")   
            genome_to_kos[name].add(first_line)
            for line in f:
                line = line.strip()
                if line:
                    genome_to_kos[name].add(line)
                    
        # Case 2: Multiple genomes (2 columns per line)
        elif number_of_fields == 2:
            id_genome, id_ko = fields
            genome_to_kos[id_genome].add(id_ko)
            for line in f:
                line = line.strip()
                if line:
                    id_genome, id_ko = line.split("\t")
                    genome_to_kos[id_genome].add(id_ko)
                    
    # Validate KO identifiers: ensure they are 6 characters long and start with 'K'
    for id_ko in set.union(*genome_to_kos.values()):
        if len(id_ko) != 6 or not id_ko.startswith("K"):
            raise ValueError("Each KO should be 6 characters long and start with 'K'.")
    
    return genome_to_kos

    
    
# Formatting
# ==========
# Get duration
def format_duration(duration):
    """
    Format the elapsed time since `t0` in hours, minutes, and seconds.
    
    Adapted from @john-fouhy:
    https://stackoverflow.com/questions/538666/python-format-timedelta-to-string
    """
    hours, remainder = divmod(int(duration), 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{hours:02}:{minutes:02}:{seconds:02}"

# Format header for printing
def format_header(text, line_character="=", n=None):
    if n is None:
        n = len(text)
    line = n*line_character
    return "{}\n{}\n{}".format(line, text, line)

# Format memory
def format_bytes(B, unit="auto", return_units=True):
    """
    Return the given bytes as a human-readable string in KB, MB, GB, or TB.
    1 KB = 1024 Bytes

    Adapted from the following source (@whereisalext):
    https://stackoverflow.com/questions/12523586/python-format-size-application-converting-b-to-kb-mb-gb-tb/52379087
    """
    KB = 1024
    MB = KB ** 2  # 1,048,576
    GB = KB ** 3  # 1,073,741,824
    TB = KB ** 4  # 1,099,511,627,776

    def format_with_unit(size, unit_name):
        return f"{size:.2f} {unit_name}" if return_units else size

    unit = unit.lower()
    if unit != "auto":
        unit = unit.lower()
        if unit == "b":
            return format_with_unit(B, "B")
        elif unit == "kb":
            return format_with_unit(B / KB, "KB")
        elif unit == "mb":
            return format_with_unit(B / MB, "MB")
        elif unit == "gb":
            return format_with_unit(B / GB, "GB")
        elif unit == "tb":
            return format_with_unit(B / TB, "TB")
        else:
            raise ValueError(f"Unknown unit: {unit}")
    else:
        if B < KB:
            return format_with_unit(B, "B")
        elif KB <= B < MB:
            return format_with_unit(B / KB, "KB")
        elif MB <= B < GB:
            return format_with_unit(B / MB, "MB")
        elif GB <= B < TB:
            return format_with_unit(B / GB, "GB")
        else:
            return format_with_unit(B / TB, "TB")
        
# Logging
# =======
def build_logger(logger_name=__name__, stream=sys.stdout):
    # Create a logger object
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)  # Set the logging level
    
    # Create a stream handler to output logs to stdout
    stream_handler = logging.StreamHandler(stream)
    stream_handler.setLevel(logging.DEBUG)  # Set the level for the handler
    
    # Create a formatter and set it to the handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    
    # Add the handler to the logger
    logger.addHandler(stream_handler)

    return logger
    
def reset_logger(logger):
    # Remove all existing handlers
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)
        handler.close()
    
    # Set a new handler (for example, to output to stdout)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    # Optionally set a new level
    logger.setLevel(logging.DEBUG)
    
# Timestamp
def get_timestamp(format_string:str="%Y-%m-%d %H:%M:%S"):
    # Get the current date and time
    now =  datetime.now()
    # Create a timestamp string
    return now.strftime(format_string)

# Directory
# =========
def get_file_size(filepath:str, format=False):
    size_in_bytes = os.stat(filepath).st_size
    if format:
        return format_bytes(size_in_bytes)
    else:
        return size_in_bytes
    
def check_file(filepath:str, empty_ok=False, minimum_filesize=1): # Doesn't handle empty gzipped files
    if not os.path.exists(filepath):
        raise FileNotFoundError(filepath)
    if not empty_ok:
        if get_file_size(filepath) < minimum_filesize:
            raise EmptyDataError(filepath)

# md5 hash from file
def get_md5hash_from_file(filepath:str, block_size=65536):
    """
    Calculate the MD5 hash of a file.

    Parameters:
    - file_path: The path to the file.
    - block_size: The size of each block read from the file (default is 64KB).

    Returns:
    - A string containing the MD5 hash.
    """
    md5 = hashlib.md5()
    with open(filepath, 'rb') as f:
        for block in iter(lambda: f.read(block_size), b''):
            md5.update(block)
    return md5.hexdigest()

# md5 hash from directory
def get_md5hash_from_directory(directory:str):
    """
    Calculate the MD5 hash of all files in a directory.

    Parameters:
    - directory_path: The path to the directory.

    Returns:
    - A dictionary where the keys are file paths and the values are their MD5 hashes.
    """
    md5_hashes = {}
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.isfile(file_path):
                file_md5 = get_md5hash_from_file(file_path)
                md5_hashes[file_path] = file_md5
    return md5_hashes

# Get directory tree structure
def get_directory_tree(root, ascii=False):
    if not ascii:
        return DisplayablePath.view(root)
    else:
        return DisplayablePath.get_ascii(root)

# Directory size
def get_directory_size(directory:str='.'):
    """
    Adapted from @Chris:
    https://stackoverflow.com/questions/1392413/calculating-a-directorys-size-using-python
    """

    total_size = 0
    seen = {}
    for dirpath, dirnames, filenames in os.walk(directory):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            try:
                stat = os.stat(fp)
            except OSError:
                continue

            try:
                seen[stat.st_ino]
            except KeyError:
                seen[stat.st_ino] = True
            else:
                continue

            total_size += stat.st_size

    return total_size


# Classes
# =======
class RunShellCommand(object):
    """
    Args: 
        command:str command to be executed
        name:str name associated with command [Default: None]
        shell_executable:str path to executable [Default: /bin/bash]
        
    Usage: 
        cmd = RunShellCommand("time (sleep 5 & echo 'Hello World')", name="Demo")
        cmd.run()
        cmd
        # ================================================
        # RunShellCommand(name:Demo)
        # ================================================
        # (/bin/bash)$ time (sleep 5 & echo 'Hello World')
        # ________________________________________________
        # Properties:
        #     - stdout: 61.00 B
        #     - stderr: 91.00 B
        #     - returncode: 0
        #     - peak memory: 37.22 B
        #     - duration: 00:00:05

    """

    def __init__(
        self, 
        command:str, 
        name:str=None, 
        shell_executable:str="/bin/bash",
        validate_input_filepaths:list=None,
        validate_output_filepaths:list=None,
        ):

        if isinstance(command, str):
            command = [command]
        command = " ".join(list(filter(bool, map(str, command))))
        self.command = command
        self.name = name
        self.shell_executable = shell_executable
        self.validate_input_filepaths = validate_input_filepaths if validate_input_filepaths else list()
        self.validate_output_filepaths = validate_input_filepaths if validate_input_filepaths else list()
        self.executed = False
        
    def run(self, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8", **popen_kws):
        def execute_command(encoding, stdout, stderr):
            # Execute the process
            self.process_ = subprocess.Popen(
                self.command,
                shell=True,
                stdout=stdout,
                stderr=stderr,
                executable=self.shell_executable,
                **popen_kws,
            )
            # Wait until process is complete and return stdout/stderr
            self.stdout_, self.stderr_ = self.process_.communicate()
            self.returncode_ = self.process_.returncode
            
            # Encode
            if encoding:
                if self.stdout_:
                    self.stdout_ = self.stdout_.decode(encoding)
                if self.stderr_:
                    self.stderr_ = self.stderr_.decode(encoding)

        # I/O
        self.redirect_stdout = None
        if isinstance(stdout, str):
            self.redirect_stdout = stdout
            stdout = open(stdout, "wb")

        self.redirect_stderr = None
        if isinstance(stderr, str):
            self.redirect_stderr = stderr
            stderr = open(stderr, "wb")

        # Measure memory usage
        t0 = time.time()
        if self.validate_input_filepaths:
            for filepath in self.validate_input_filepaths:
                check_file(filepath, empty_ok=False)
        self.memory_usage_ = memory_usage((execute_command, (encoding, stdout, stderr,)), max_iterations=1)
        self.duration_ = time.time() - t0

        # # Flush
        # if hasattr(stdout, "flush"):
        #     stdout.flush()
        # if hasattr(stderr, "flush"):
        #     stderr.flush()
            
        # Close
        if hasattr(stdout, "close"):
            stdout.close()
        if hasattr(stderr, "close"):
            stderr.close()

        self.peak_memory_ = max(self.memory_usage_)
        self.executed = True

        return self

    def __repr__(self):
        name_text = "{}(name:{})".format(self.__class__.__name__, self.name)
        command_text = "({})$ {}".format(self.shell_executable, self.command)
        n = max(len(name_text), len(command_text))
        pad = 4
        fields = [
            format_header(name_text,line_character="=", n=n),
            *format_header(command_text, line_character="_", n=n).split("\n")[1:],
            ]
        if self.executed:
            fields += [
            "Properties:",
            ]
            # stdout
            if self.redirect_stdout:
                fields += [
                pad*" " + "- stdout({}): {}".format(
                    self.redirect_stdout,
                    get_file_size(self.redirect_stdout, format=True),
                )
                ]
            else:
                fields += [
                pad*" " + "- stdout: {}".format(format_bytes(sys.getsizeof(self.stdout_))),
                ]
            # stderr
            if self.redirect_stderr:
                fields += [
                pad*" " + "- stderr({}): {}".format(
                    self.redirect_stderr,
                    get_file_size(self.redirect_stderr, format=True),
                )
                ]
            else:
                fields += [
                pad*" " + "- stderr: {}".format(format_bytes(sys.getsizeof(self.stderr_))),
                ]

            fields += [
            pad*" " + "- returncode: {}".format(self.returncode_),
            pad*" " + "- peak memory: {}".format(format_bytes(self.peak_memory_)),
            pad*" " + "- duration: {}".format(format_duration(self.duration_)),
            ]
        return "\n".join(fields)
    
    # Dump stdout, stderr, and returncode
    def dump(self, output_directory:str):    
        # stdout
        with open_file_writer(os.path.join(output_directory, f"{self.name}.o")) as f:
            print(self.stdout_, file=f)
        # stderr
        with open_file_writer(os.path.join(output_directory, f"{self.name}.e")) as f:
            print(self.stderr_, file=f)
        # returncode
        with open_file_writer(os.path.join(output_directory, f"{self.name}.returncode")) as f:
            print(self.returncode_, file=f)
            
    # Check status
    def check_status(self):
        if self.returncode_ != 0:
            raise subprocess.CalledProcessError(
                returncode=self.returncode_,
                cmd="\n".join([
                f"Command Failed: {self.command}",
                f"return code: {self.returncode_}",
                f"stderr:\n{self.stderr_}",
                ]),
            )
        else:
            if self.validate_output_filepaths:
                for filepath in self.validate_output_filepaths:
                    check_file(filepath, empty_ok=False)
            print(f"Command Successful: {self.command}", file=sys.stderr)

# Check argument choices
def check_argument_choice(query, choices:set):
    """_summary_

    Args:
        query (_type_): Query option
        choices (set): Acceptable options

    Raises:
        ValueError: _description_
    """
    choices = set(choices)
    if query not in choices:
        raise ValueError(f"Invalid option '{query}'. Allowed choices are: {choices}")