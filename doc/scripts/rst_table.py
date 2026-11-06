"""Format documentation tables using Astropy's RST writer."""

from io import StringIO

from astropy.table import Table


def format_rst_table(data_table, has_header=True):
    """Format a two-dimensional array as an RST simple table.

    For a headerless table, AstroPy sizes the columns using the first data
    row as its header.  Removing only the header/body rule then leaves a
    valid headerless simple table with all data rows intact.
    """
    table = Table(rows=data_table[1:], names=data_table[0])
    stream = StringIO()
    table.write(stream, format='ascii.rst')
    lines = stream.getvalue().splitlines()
    if not has_header:
        del lines[2]
    # The writer pads empty final cells with spaces; omit trailing whitespace
    # from the generated documentation without changing the table layout.
    return '\n'.join(line.rstrip() for line in lines) + '\n'
