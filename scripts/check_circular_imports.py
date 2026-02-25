#!/usr/bin/env python3
"""Detect circular imports within the pypeit package.

Usage: python scripts/check_circular_imports.py
"""
import ast
import os
from pathlib import Path
from collections import defaultdict, deque

ROOT = Path(__file__).resolve().parents[1] / "pypeit"


def module_name_from_path(path: Path) -> str:
    rel = path.relative_to(ROOT.parent)
    parts = rel.with_suffix("").parts
    return ".".join(parts)


def parse_imports(path: Path):
    src = path.read_text()
    tree = ast.parse(src)
    imports = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for n in node.names:
                imports.append((n.name, None, None))
        elif isinstance(node, ast.ImportFrom):
            mod = node.module
            level = node.level
            names = [n.name for n in node.names]
            imports.append((mod, level, names))
    return imports


def resolve_import(current_mod: str, imp_mod, level):
    # imp_mod: may be None (e.g., "from . import x"), level: int or None
    if level and level > 0:
        cur_parts = current_mod.split('.')
        # drop 'pypeit' parent context to allow relative resolution
        base = cur_parts[:]
        for _ in range(level):
            if base:
                base = base[:-1]
        if imp_mod:
            return ".".join(base + imp_mod.split('.'))
        else:
            return ".".join(base)
    else:
        return imp_mod


def build_graph(root: Path):
    graph = defaultdict(set)
    modules = {}
    for p in root.rglob("*.py"):
        mod = module_name_from_path(p)
        modules[mod] = p

    for mod, p in modules.items():
        for imp_mod, level, names in parse_imports(p):
            if imp_mod is None:
                resolved_base = resolve_import(mod, None, level)
            else:
                resolved_base = resolve_import(mod, imp_mod, level)
            if not resolved_base:
                continue
            # If the import specifies names (from X import a, b), try to resolve
            # to submodules X.a if such a module file exists in our tree. Otherwise
            # fall back to the base module/package.
            if names:
                for name in names:
                    # skip wildcard imports
                    if name == '*':
                        # treat as dependency on the base package
                        target = resolved_base
                        if target == 'pypeit' or target.startswith('pypeit.'):
                            graph[mod].add(target)
                        continue
                    candidate = f"{resolved_base}.{name}" if resolved_base else name
                    if candidate in modules:
                        target = candidate
                    else:
                        target = resolved_base
                    if target and (target == 'pypeit' or target.startswith('pypeit.')):
                        graph[mod].add(target)
            else:
                if resolved_base == 'pypeit' or resolved_base.startswith('pypeit.'):
                    graph[mod].add(resolved_base)

    return graph, modules


def find_cycles(graph):
    # detect cycles using DFS
    visited = set()
    stack = []
    cycles = []

    def dfs(node, pathset):
        visited.add(node)
        pathset.add(node)
        stack.append(node)
        for nb in graph.get(node, []):
            if nb not in visited:
                dfs(nb, pathset)
            elif nb in pathset:
                # found cycle
                try:
                    i = stack.index(nb)
                    cyc = stack[i:] + [nb]
                except ValueError:
                    cyc = [nb]
                cycles.append(cyc)
        pathset.remove(node)
        stack.pop()

    for n in list(graph.keys()):
        if n not in visited:
            dfs(n, set())

    # canonicalize cycles
    uniq = set()
    out = []
    for c in cycles:
        key = tuple(c)
        if key not in uniq:
            uniq.add(key)
            out.append(c)
    return out


def main():
    graph, modules = build_graph(ROOT)
    cycles = find_cycles(graph)
    if not cycles:
        print("No cycles detected in pypeit package imports.")
        return 0
    print(f"Detected {len(cycles)} cycle(s):\n")
    for cyc in cycles:
        print(" -> ".join(cyc))
        print()
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
