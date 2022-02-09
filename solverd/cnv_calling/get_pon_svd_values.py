import os
import sys
import h5py


indirfiles = os.listdir(sys.argv[1])
ponfiles = [ponfile for ponfile in indirfiles if ponfile.endswith(".hdf5")]

for ponfile in ponfiles:
    ponname = ponfile.split(".")[0]
    outfilepath = f"{sys.argv[2]}/svd_{ponname}.txt"
    hdf5pon = h5py.File(f"{sys.argv[1]}/{ponfile}", 'r')
    pon_svd = hdf5pon["/panel/singular_values"].value

    try:
        with open(outfilepath, 'w') as outfile:
            for svdvalue in pon_svd:
                outfile.write(f"{svdvalue}\n")
    except IOError:
        print("AAP")
