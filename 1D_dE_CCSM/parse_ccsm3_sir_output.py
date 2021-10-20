"""Script to parse output files"""

def parse_assignment(line):
    """Parses a line of form x = y, and returns y"""
    return float(line.split('=')[1].strip())


def main(filename):
    """Does the job"""
    with open(filename, 'r') as f:
        lines = f.readlines()

#   Day of year line 6
    print(parse_assignment(lines[5]))
#   Latitude line 7
#    for line in lines[0:10]:
#        print(line)

    return


if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
