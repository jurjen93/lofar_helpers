import argparse
import re
import shutil
from pathlib import Path


def main(source_path, out_path):
    pattern = r'^([0-9.]+)_([0-9.]+)_([0-1]+)_(.+)$'

    for img in Path(source_path).glob('*.png'):
        name = img.stem
        match = re.match(pattern, name)

        if not match:
            print(f"Regex didn't match on {name}!")
            continue

        uncertaincy, pred, label, name = match.groups()

        if abs(float(pred) - float(label)) > 0.4:
            print(f"yes! {pred}, {label}, {name}")
            shutil.copy(img, out_path)
        else:
            print(f"no! {pred}, {label}, {name}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # Add the argument
    parser.add_argument('source_path', type=str, help='The path to the classified files')
    args = parser.parse_args()

    out_path = Path(args.source_path)
    out_path = out_path.with_name(out_path.name + '_misclassifications')
    out_path.mkdir()
    print(out_path)

    main(args.source_path, out_path)