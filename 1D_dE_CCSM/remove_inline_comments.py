"""Removes inline comments from ccsm3 fortran code"""

def strip_inline_comment(s):
    """Removes lines with inline comments and returns
    new line right-stripped"""
    if ("!" in s) & (not s.startswith(tuple("!Cc"))):
        idx = s.find("!")
        return s[:idx].rstrip()+"\n"
    else:
        return s


def remove_inline_comments(filepath):

    with open(filepath) as f:
        lines = f.readlines()

    with open(filepath.replace('.for', '.no_inline_comment.for'), "w") as of:
        for idx, line in enumerate(lines):
            of.write(strip_inline_comment(line))


if __name__ == "__main__":
    filepath = 'ccsm3_sir_de.for'
    remove_inline_comments(filepath)
