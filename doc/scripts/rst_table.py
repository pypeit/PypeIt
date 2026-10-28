"""Format documentation tables using Astropy's RST writer."""

from io import StringIO

from astropy.table import Table


def format_rst_table(data_table):
    """Format a two-dimensional array whose first row contains column names."""
    table = Table(rows=data_table[1:], names=data_table[0])
    stream = StringIO()
    table.write(stream, format='ascii.rst')
    # The writer pads empty final cells with spaces; omit trailing whitespace
    # from the generated documentation without changing the table layout.
    return '\n'.join(line.rstrip() for line in stream.getvalue().splitlines()) + '\n'
