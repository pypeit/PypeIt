#!/usr/bin/env python3
"""
Reflow every multi-line ``descr=`` parameter description passed to
:func:`~pypeit.par.parset.set_parameter_definition` so each line fills as
much of the target width as possible (99 columns by default).

This mirrors ``reflow_docstrings.py`` (a similar tool for triple-quoted
docstrings, written for a different repository) but is adapted for
``descr=(...)`` implicit string concatenation instead: a tuple of adjacent
string literals, rather than one triple-quoted block.

Two cases are handled differently:

  - Plain (non f-string) multi-line ``descr`` values: reflowed as ordinary
    string literals, preserving ``r'...'`` raw-string style wherever the
    text contains a backslash (e.g. a ``:math:`...``` span), since
    ``repr()`` would otherwise silently re-escape it and drop the ``r``
    prefix.
  - f-string multi-line ``descr`` values (``ast.JoinedStr``): each
    ``{expr}`` interpolation is treated as one atomic, never-split token
    (using its exact original source text), while the literal text around
    it is split into ordinary fillable words.

By default this only reports what would change; pass ``--apply`` to write
the changes.
"""
import argparse
import ast
import re
import sys
from pathlib import Path


def smart_repr(line):
    """
    Represent `line` as a Python string literal, preferring an
    ``r'...'`` / ``r"..."`` raw-string form whenever `line` contains a
    backslash, so re-wrapping doesn't silently turn a raw string back
    into an escaped one (``repr()`` always escapes backslashes and never
    re-emits the ``r`` prefix).

    Parameters
    ----------
    line : str
        The text to represent as a string literal.

    Returns
    -------
    str
        A valid Python string-literal source for `line`.
    """
    if '\\' in line and not line.endswith('\\'):
        if "'" not in line:
            return "r'" + line + "'"
        if '"' not in line:
            return 'r"' + line + '"'
    return repr(line)


def escape_for_single_quote(s):
    """Escape `s` for embedding inside a plain ``'...'`` string literal."""
    return s.replace('\\', '\\\\').replace("'", "\\'")


def tokenize_plain(text):
    """
    Split plain text into ``('text', word_with_trailing_spaces)`` tokens.

    A leading run of spaces (e.g. text immediately following a ``{expr}``
    interpolation) is kept as its own token rather than silently dropped
    -- ``\\S+`` alone would require a match to start on a non-space
    character.

    Parameters
    ----------
    text : str
        The literal text to tokenize.

    Returns
    -------
    list
        ``('text', substring)`` tokens that concatenate back to `text`.
    """
    return [('text', tok) for tok in re.findall(r'^ +|\S+ *', text)]


def tokenize_joined(node, src):
    """
    Turn an ``ast.JoinedStr`` (f-string) into a flat list of
    ``('text', ...)`` / ``('expr', ...)`` tokens: literal segments are
    word-split like :func:`tokenize_plain`; each ``{expression}`` becomes
    one atomic ``('expr', '{...}')`` token built from its exact original
    source text (so nested quoting, calls, etc. survive untouched).

    Parameters
    ----------
    node : ast.JoinedStr
        The f-string node to tokenize.
    src : str
        The full source text `node` was parsed from (needed to recover
        the exact text of each ``{expr}``).

    Returns
    -------
    list
        A mix of ``('text', ...)`` and ``('expr', ...)`` tokens.
    """
    tokens = []
    for part in node.values:
        if isinstance(part, ast.Constant):
            tokens.extend(tokenize_plain(part.value))
        elif isinstance(part, ast.FormattedValue):
            expr_src = ast.get_source_segment(src, part.value)
            conv = f'!{chr(part.conversion)}' if part.conversion != -1 else ''
            spec = ''
            if part.format_spec is not None:
                spec_src = ast.get_source_segment(src, part.format_spec)
                spec = f':{spec_src}'
            tokens.append(('expr', '{' + expr_src + conv + spec + '}'))
        else:
            raise ValueError(f'unexpected JoinedStr part: {part!r}')
    return tokens


def pack_lines(tokens, indent, width=99):
    """
    Greedily pack a token stream (from :func:`tokenize_plain` or
    :func:`tokenize_joined`) into physical lines, each rendered as a
    Python string literal (plain, raw, or f-string as needed), with
    total column width (indentation plus quotes/prefix) at most `width`.

    An ``'expr'`` token is never split and always forces its line to be
    an f-string. A line containing a backslash needs an ``r'...'`` /
    ``r"..."`` prefix (see :func:`smart_repr`) -- that extra character is
    accounted for here so the final rendered line can never overflow
    `width`.

    Parameters
    ----------
    tokens : list
        Tokens as produced by :func:`tokenize_plain` or
        :func:`tokenize_joined`.
    indent : int
        The column the string literal(s) start at in the source.
    width : int, optional
        Target total column width (indentation included).

    Returns
    -------
    list
        One already-indented Python string-literal source line per
        entry.
    """
    avail_base = width - indent - 2  # 2 for the quote chars

    def prefix_len(has_expr, has_backslash):
        # 'f' or 'r' -- both single characters; expr (f-string) takes
        # precedence since that's the more common/deliberate case here.
        return 1 if (has_expr or has_backslash) else 0

    lines = []          # list of (list-of-tokens, has_expr)
    cur = []
    cur_len = 0
    cur_has_expr = False
    cur_has_backslash = False
    for kind, val in tokens:
        tok_len = len(val)
        tok_has_backslash = '\\' in val
        new_has_expr = cur_has_expr or kind == 'expr'
        new_has_backslash = cur_has_backslash or tok_has_backslash
        new_prefix = prefix_len(new_has_expr, new_has_backslash)
        if cur and cur_len + tok_len + new_prefix > avail_base:
            lines.append((cur, cur_has_expr))
            cur, cur_len, cur_has_expr, cur_has_backslash = [], 0, False, False
            new_has_expr = kind == 'expr'
            new_has_backslash = tok_has_backslash
        cur.append((kind, val))
        cur_len += tok_len
        cur_has_expr = new_has_expr
        cur_has_backslash = new_has_backslash
    if cur:
        lines.append((cur, cur_has_expr))

    pad = ' ' * indent
    out = []
    for line_tokens, has_expr in lines:
        if has_expr:
            body = ''.join(
                v if k == 'expr' else escape_for_single_quote(v)
                for k, v in line_tokens
            )
            out.append(pad + "f'" + body + "'")
        else:
            text = ''.join(v for _, v in line_tokens)
            out.append(pad + smart_repr(text))
    return out


def find_descr_nodes(tree):
    """
    Find every ``descr=`` keyword argument in the parsed module whose
    value spans more than one source line.

    Parameters
    ----------
    tree : ast.Module
        The parsed source, from `ast.parse`.

    Returns
    -------
    list of ast.expr
        Each entry is the ``descr`` keyword's value node (a
        :class:`ast.Constant` or :class:`ast.JoinedStr`).
    """
    nodes = []
    for node in ast.walk(tree):
        if isinstance(node, ast.keyword) and node.arg == 'descr':
            v = node.value
            if v.lineno != v.end_lineno:
                nodes.append(v)
    return nodes


def process_file(path, width=99, dry_run=True):
    """
    Reflow every multi-line ``descr=`` value in the Python file at
    `path`.

    Parameters
    ----------
    path : str or pathlib.Path
        The file to process.
    width : int, optional
        Target total column width (indentation included).
    dry_run : bool, optional
        If True (the default), only print what would change; if False,
        write the changes to `path`.

    Returns
    -------
    int
        The number of ``descr`` values changed (or that would change, in
        a dry run).
    """
    with open(path, 'r') as f:
        src = f.read()
    tree = ast.parse(src)

    replacements = []
    for v in find_descr_nodes(tree):
        if isinstance(v, ast.JoinedStr):
            tokens = tokenize_joined(v, src)
        else:
            try:
                text = ast.literal_eval(v)
            except Exception:
                continue
            if not isinstance(text, str):
                continue
            tokens = tokenize_plain(text)

        indent = v.col_offset
        new_lines = pack_lines(tokens, indent, width=width)
        # The first physical line is prepended, at apply time, to the
        # existing source text up through col_offset -- which is itself
        # that same indentation -- so it must not carry its own copy.
        new_lines = [new_lines[0][indent:]] + new_lines[1:]
        new_src = '\n'.join(new_lines)

        old_src = ast.get_source_segment(src, v)
        if old_src is None or new_src == old_src:
            continue

        replacements.append((v.lineno, v.col_offset, v.end_lineno, v.end_col_offset,
                              old_src, new_src))

    if dry_run:
        for (sl, sc, el, ec, old, new) in replacements:
            print(f'--- {path}, lines {sl}-{el} ---')
            print('BEFORE:')
            print(old)
            print('AFTER:')
            print(new)
            print()
        print(f'{len(replacements)} descr value(s) would change in {path}')
        return len(replacements)

    # Apply replacements from bottom to top so earlier offsets stay valid.
    src_lines = src.split('\n')
    by_position = sorted(replacements, key=lambda r: (r[0], r[1]), reverse=True)
    for (sl, sc, el, ec, old, new) in by_position:
        if sl == el:
            line = src_lines[sl - 1]
            src_lines[sl - 1] = line[:sc] + new + line[ec:]
        else:
            first_line = src_lines[sl - 1][:sc] + new
            last_line_tail = src_lines[el - 1][ec:]
            new_block = first_line.split('\n')
            new_block[-1] = new_block[-1] + last_line_tail
            src_lines[sl - 1:el] = new_block
    with open(path, 'w') as f:
        f.write('\n'.join(src_lines))
    print(f'{len(replacements)} descr value(s) changed in {path}')
    return len(replacements)


def iter_python_files(paths):
    """
    Expand a list of file/directory arguments into the ``.py`` files
    they name, recursing into directories.

    Parameters
    ----------
    paths : list of str
        Files and/or directories, as given on the command line.

    Returns
    -------
    list of pathlib.Path
        Every ``.py`` file found, in sorted order, with `__pycache__`
        contents excluded.
    """
    files = []
    for raw in paths:
        p = Path(raw)
        if p.is_dir():
            files.extend(
                f for f in p.rglob('*.py') if '__pycache__' not in f.parts
            )
        else:
            files.append(p)
    return sorted(set(files))


def main():
    """Parse command-line arguments and reflow the requested files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'paths', nargs='+',
        help='Python files and/or directories (searched recursively) to reflow.'
    )
    parser.add_argument(
        '--width', type=int, default=99,
        help='Target column width for reflowed descr= text.'
    )
    parser.add_argument(
        '--apply', action='store_true',
        help='Write changes to disk. Without this, only reports what would change.'
    )
    args = parser.parse_args()

    total = 0
    for path in iter_python_files(args.paths):
        total += process_file(path, width=args.width, dry_run=not args.apply)
    if not args.apply:
        print(f'\n{total} descr value(s) total would change (dry run; pass --apply to write).')
    else:
        print(f'\n{total} descr value(s) total changed.')


if __name__ == '__main__':
    main()
