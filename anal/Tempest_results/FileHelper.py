# ===============================================================================
# Helpful class and functions for dealing with files and directories
# ===============================================================================

class File(file):
    """ An helper class for file reading  """

    def __init__(self, *args, **kwargs):
        super(File, self).__init__(*args, **kwargs)
        self.BLOCKSIZE = 4096

    def head(self, lines_2find=1):
        """Return first few lines of a file"""
        self.seek(0) # Rewind file
        return [super(File, self).next() for x in xrange(lines_2find)]

    def tail(self, lines_2find=1):
        """Return last few lines of a file"""
        self.seek(0, 2) # Go to end of file
        bytes_in_file = self.tell()
        lines_found, total_bytes_scanned = 0, 0
        while (lines_2find + 1 > lines_found and
               bytes_in_file > total_bytes_scanned): 
            byte_block = min(
                self.BLOCKSIZE,
                bytes_in_file - total_bytes_scanned)
            self.seek( -(byte_block + total_bytes_scanned), 2)
            total_bytes_scanned += byte_block
            lines_found += self.read(self.BLOCKSIZE).count('\n')
        self.seek(-total_bytes_scanned, 2)
        line_list = list(self.readlines())
        return line_list[-lines_2find:]

    def backward(self):
        """Read file in reverse"""
        self.seek(0, 2) # Go to end of file
        blocksize = self.BLOCKSIZE
        last_row = ''
        while self.tell() != 0:
            try:
                self.seek(-blocksize, 1)
            except IOError:
                blocksize = self.tell()
                self.seek(-blocksize, 1)
            block = self.read(blocksize)
            self.seek(-blocksize, 1)
            rows = block.split('\n')
            rows[-1] = rows[-1] + last_row
            while rows:
                last_row = rows.pop(-1)
                if rows and last_row:
                    yield last_row
        yield last_row
# ===============================================================================

def list_files(dir, ext):
    """get list files in dir that match ext"""
    import os, fnmatch

    r = []
    for root, dirs, files in os.walk(dir):
        for name in fnmatch.filter(files, ext):
            r.append(os.path.join(root, name))
    return r
# ===============================================================================

def numerical_sort(value):
    """custom sorting function for files"""
    import re

    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts
# ===============================================================================

def get_immediate_subdirectories(a_dir):
    """get only the immediate subdirectories of a directory"""
    import os

    subdirs = []
    for l in os.listdir(a_dir):   
        name = os.path.join(a_dir, l)
        if os.path.isdir(name):
            subdirs.append(name)

    return subdirs
# ===============================================================================

def get_files_with_extension(a_dir, a_ext):
    """get all the files in a directory with a given extension"""
    import os

    files = []
    for l in os.listdir(a_dir):
        name = os.path.join(a_dir, l)
        if name.endswith(a_ext):
            files.append(name)

    return files
# ===============================================================================
