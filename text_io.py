import sys   # system environment

# ----------------------------------------------------------------------
# function to grep lines between two flag lines
#    start1 and start2 are the starting position of flag1 and flag2 in the line, respectively
#    
def grab_lines_between(input, output, flag1, flag2, start1 = None, start2 = None):
    try:
        if start1:
            offset = len(flag1)
            end1 = start1 + offset
            while 1:
                line = input.readline()
                if flag1 == line[start1:end1]: 
                    output.write(line)
                    break
        else:
            while 1:
                line = input.readline()
                if flag1 in line:
                    output.write(line)
                    break
        if start2:
            offset = len(flag2)
            end2 = start2 + offset
            while 1:
                line = input.readline()
                output.write(line)
                if flag2 == line[start2:end2]: break
        else:
            while 1:
                line = input.readline()
                output.write(line)
                if flag2 in line: break
    except EOFError:
        print('Error: met unexpected end-of-file')
        sys.exit(1)

# ----------------------------------------------------------------------
# function to skip specified number of lines
def skip_lines(input,num_lines):
    try:
        for i in range(num_lines):
            input.readline()
    except EOFError:
        print('Error: met unexpected end-of-file')
        sys.exit(1)

# ----------------------------------------------------------------------
# function to skip specified number of lines, blank line does not count
def skip_real_lines(file, num_lines, comment_tag=None):
    try:
        i = 0
        while 1:
            line = read_next_line(file)
            if comment_tag: line = line.split(comment_tag)[0]
            if len(line.strip()) > 0: i += 1
            if i == num_lines: break 
    except EOFError:
        print('Error: met unexpected end-of-file')
        sys.exit(1)
        
# ----------------------------------------------------------------------
# function to skip to specified line
def skip_to(file,flag,start_position = None):
    if start_position is None:
        while 1:
            line = file.readline()
            if line == '':
                print('Error: met unexpected end-of-file')
                return -1
            else:
                if flag in line: break
    else:
        offset = len(flag)
        end_position = start_position + offset
        while 1:
            line = file.readline()
            if line == '':
                print('Error: met unexpected end-of-file')
                return -1
            else:
                if flag == line[start_position:end_position]: break
    return line


# ----------------------------------------------------------------------
# copy nlines data (starting from the current position) from a file to another 
def copy_lines(from_file, to_file, nlines):
    for i in range(nlines):
        line = from_file.readline()
        to_file.write(line)

# ----------------------------------------------------------------------
def read_next_line(f, comment_tag=None):
    ''' read next line of data from text file (f=file handler), ignore blank or comment line
        content after comment_tag is considered comment
        return False when unexpected End-of-file is met
    '''
    while 1:
        try:
            line = f.readline().split(comment_tag)[0].strip() if comment_tag else f.readline().strip()
            if len(line) == 0: 
                continue
            else:
                return line
        except EOFError:
            return False

# ----------------------------------------------------------------------
def convert_str2num(data):
    ''' convert numeric (integer or floating number) string to values
        data can be scalar or a list of numeric string
    '''
    if isinstance(data, list):
        val = []
        for item in data:
            try:
                val.extend([int(item)])
            except ValueError:
                val.extend([float(item)])
    else:
        try:
            val = int(data)
        except ValueError:
            val = float(data)
    return val

        
