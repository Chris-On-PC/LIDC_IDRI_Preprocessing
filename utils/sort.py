import os
import h5py
import shutil
import sys

def sort_data(path, min_annotations = 4):
    full_path = os.path.join(path, 'selected')

    if not os.path.isdir(full_path):
        os.mkdir(full_path)

    for name in os.listdir(path):
        if '.h5' in name:
            filename = os.path.join(path, name)
            f = h5py.File(filename, 'r')
            num_keys = len(f.keys()) 
            f.close()
            print(f"{num_keys=}")
            if num_keys == min_annotations+1:
                shutil.move(filename, os.path.join(full_path, name))

if __name__ == "__main__":
    path = sys.argv[1]
    min_annoations = sys.argv[2]

    sort_data(path, min_annoations)