import os

files = os.listdir('.')
for filename in files:
    portion = os.path.splitext(filename)
    print(portion)
    if portion[1] == ".ply":
        newname = portion[0]+'.3d'
        os.rename(filename, newname)

