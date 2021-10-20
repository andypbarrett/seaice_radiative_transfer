"""Script to parse output files"""

def main(filename):
    """Does the job"""
    with open(filename, 'r') as f:
        lines = f.readlines()

    for line in lines[0:10]:
        print(lines)

    return


if __name__ == '__main__':
    filename = 'ccsm3_sir_de_output.dat'
    main(filename)
