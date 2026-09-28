#!/usr/bin/env python3
"""Copy the shared integration into independent pipeline repositories."""
import argparse
from pathlib import Path
import shutil

SOURCE = Path(__file__).resolve().parents[1] / 'integrations' / 'nextflow'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('pipelines', nargs='+', type=Path)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    for pipeline in args.pipelines:
        for source in sorted(SOURCE.rglob('*')):
            if not source.is_file():
                continue
            target = pipeline / source.relative_to(SOURCE)
            if args.check:
                if not target.is_file() or target.read_bytes() != source.read_bytes():
                    parser.exit(1, f'Integration differs: {target}\n')
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source, target)


if __name__ == '__main__':
    main()
