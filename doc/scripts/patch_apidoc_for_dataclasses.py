#!/usr/bin/env python3
"""
Patch generated sphinx-apidoc module pages so selected dataclasses are
documented with ``:no-undoc-members:`` while preserving
``:undoc-members:`` for the rest of each module.

This script is intended to run *after* ``sphinx-apidoc`` regenerates the
``doc/api/*.rst`` files.

For each entry in ``PATCHES``:

- the module-level ``.. automodule::`` directive is updated to exclude the
  listed dataclass names via ``:exclude-members:``
- a corresponding ``.. autoclass::`` block is appended for each dataclass,
  using ``:members:`` and ``:no-undoc-members:``

The goal is to avoid ugly/duplicated dataclass-field rendering in API docs
while still keeping project-wide ``undoc-members`` enabled.
"""

from __future__ import annotations

from pathlib import Path
import re
import sys

from importlib import resources

# Mapping of @dataclass names that should be documented separately with
#   :no-undoc-members:.
PATCHES: dict[str, list[str]] = {
    "pypeit.coadd2d": [
        "CoAdd2dStack",
    ],
    # Other @dataclasses go here...
}


# Define the template
AUTOMODULE_RE_TEMPLATE = r"""^
(?P<directive>\.\.\s+automodule::\s+{module}\s*\n)
(?P<options>(?:^[ \t]+:[^\n]*\n)*)
"""

# Define the API directory
API_DIR = resources.files("pypeit").parent / "doc" / "api"


def normalize_class_names(class_names: list[str]) -> list[str]:
    """Return class names in stable order with duplicates removed

    _extended_summary_

    Parameters
    ----------
    class_names : list
        List of class names to be normalized

    Returns
    -------
    list
        List in stable order with duplicates removed
    """
    seen: set[str] = set()
    out: list[str] = []
    for name in class_names:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def parse_exclude_members(options_block: str) -> tuple[list[str], list[str]]:
    """
    Split an automodule options block into its non-exclude lines and parsed
    exclude-members entries.
    """
    other_lines: list[str] = []
    exclude_members: list[str] = []

    for line in options_block.splitlines(keepends=True):
        stripped = line.strip()
        if stripped.startswith(":exclude-members:"):
            _, value = stripped.split(":", 2)[1:]
            # After splitting on ':', `value` is " member1, member2"
            members_str = value.strip()
            if members_str:
                exclude_members.extend(
                    [m.strip() for m in members_str.split(",") if m.strip()]
                )
        else:
            other_lines.append(line)

    return other_lines, exclude_members


def build_options_block(
    other_lines: list[str],
    exclude_members: list[str],
) -> str:
    """Rebuild the automodule options block

    Reconstruct the automodule block to exclude the @dataclass objects

    Parameters
    ----------
    other_lines : list
        Lines to keep in the options block, as parsed by :func:`parse_exclude_members`
    exclude_members : list
        List of @dataclass objects to exclude

    Returns
    -------
    str
        New automodule options block
    """
    block = "".join(other_lines)
    if exclude_members:
        exclude_line = "   :exclude-members: " + ", ".join(exclude_members) + "\n"
        block += exclude_line
    return block


def build_autoclass_block(module_name: str, class_name: str) -> str:
    """Construct the autoclass block for one patched dataclass

    This displays the proper level of information in the rendered output,
    excluding the ``__init__()`` function, and removing the extra output
    for each "undoc" member.

    Parameters
    ----------
    module_name : str
        Name of the module the class lives in
    class_name : str
        Name of the class

    Returns
    -------
    str
        The output autoclass block to be inserted into the .rst
    """
    return (
        f".. autoclass:: {module_name}.{class_name}\n"
        "   :members:\n"
        "   :no-undoc-members:\n"
        "   :exclude-members: __init__\n"
        "   :show-inheritance:\n"
    )


def patch_module_file(rst_path: Path, module_name: str, class_names: list[str]) -> bool:
    """Patch one generated module rst file

    _extended_summary_

    Parameters
    ----------
    rst_path : Path
        _description_
    module_name : str
        _description_
    class_names : list[str]
        _description_

    Returns
    -------
    bool
        True if the file was modified, False otherwise.

    Raises
    ------
    FileNotFoundError
        _description_
    RuntimeError
        _description_
    """
    if not rst_path.is_file():
        raise FileNotFoundError(f"Generated API file not found: {rst_path}")

    text = rst_path.read_text(encoding="utf-8")
    pattern = re.compile(
        AUTOMODULE_RE_TEMPLATE.format(module=re.escape(module_name)),
        re.MULTILINE | re.DOTALL | re.VERBOSE,
    )
    match = pattern.search(text)
    if match is None:
        raise RuntimeError(
            f"Could not find automodule directive for {module_name} in {rst_path}"
        )

    directive = match.group("directive")
    options_block = match.group("options")

    other_lines, existing_excludes = parse_exclude_members(options_block)
    merged_excludes = normalize_class_names(existing_excludes + class_names)
    new_options_block = build_options_block(other_lines, merged_excludes)

    new_automodule_block = directive + new_options_block
    old_automodule_block = match.group(0)

    # Remove any existing autoclass blocks for these class names so the script
    # is idempotent.
    new_text = text
    for class_name in class_names:
        autoclass_re = re.compile(
            rf"^"
            rf"\.\.\s+autoclass::\s+{re.escape(module_name)}\.{re.escape(class_name)}\s*\n"
            rf"(?:^[ \t]+:[^\n]*\n)*"
            rf"(?:\n)?",
            re.MULTILINE | re.DOTALL,
        )
        new_text = autoclass_re.sub("", new_text)

    new_text = new_text.replace(old_automodule_block, new_automodule_block, 1)

    autoclass_blocks = (
        "\n".join(
            build_autoclass_block(module_name, class_name) for class_name in class_names
        ).rstrip()
        + "\n"
    )

    if autoclass_blocks not in new_text:
        if not new_text.endswith("\n"):
            new_text += "\n"
        new_text += "\n" + autoclass_blocks

    if new_text != text:
        rst_path.write_text(new_text, encoding="utf-8")
        return True

    return False


# Main Driver Function =======================================================#
def main() -> int:
    """Run all configured patches"""
    modified = 0

    for module_name, class_names in PATCHES.items():
        class_names = normalize_class_names(class_names)
        rst_path = API_DIR / f"{module_name}.rst"
        changed = patch_module_file(rst_path, module_name, class_names)
        if changed:
            modified += 1
            print(f"Patched {rst_path}")
        else:
            print(f"No changes needed for {rst_path}")

    print(f"Done. Modified {modified} file(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
