import numpy as np
from astropy.io import fits
from argparse import ArgumentParser
import os
from glob import glob

def flipY(image):
    """
    Utility function to mirror an image along its Y-direction (across the X-axis). The argument is modified in-place.

    Args:
        image(2D array): The image data to be flipped.
    """
    for i in range(len(image)//2):
        tmp = np.copy(image[i,:])
        image[i,:] = image[-i-1,:]
        image[-i-1,:] = tmp

#Assumes that the files you are attempting to crop have not been previously cropped.
def cropFitsFile(hdu, x, y, width, height):
    hdu.data[np.isnan(hdu.data)] = -1
    newHDU = fits.ImageHDU(hdu.data[y:y+height,x:x+width], hdu.header)
    newHDU.header["NAXIS1"] = width
    newHDU.header["NAXIS2"] = height
    newHDU.header["TRIMSEC"] = f"[{y}:{y+height},{x}:{x+width}]"
    return newHDU

def trimFiles(x = None, y = None, width = None, height = None, root = None, output = None, recursive = False):
    if root is None:
        root = os.getcwd()
    if output is None:
        output = os.getcwd()
    files = glob(("**/" if recursive else "") + "*.fits", root_dir=root, recursive=recursive)
    if len(files) == 0:
        print("Warning: no files found to be trimmed in directory: " + root)
        return
    if not os.path.isdir(output):
        os.mkdir(output)
    print(f"Found {len(files)} files to be trimmed, trimming and saving files...")

    useDefaults = x is None or y is None or width is None or height is None

    for f in files:
        with fits.open(os.path.join(root,f)) as hdul:
            TrimmedHDUL = fits.HDUList()
            for hdu in hdul:
                if useDefaults:
                    if hdu.header["FILTER2"].strip().lower() == "vph30" and hdu.header["CAMNAME"].strip().lower() == "hk":
                        TrimmedHDUL.append(cropFitsFile(hdu, 2070, 1600, 200, 500))
                    if hdu.header["FILTER2"].strip().lower() == "vph300" and hdu.header["CAMNAME"].strip().lower() == "hk":
                        TrimmedHDUL.append(cropFitsFile(hdu, 1000, 1600, 2000, 500))
                else:
                    TrimmedHDUL.append(cropFitsFile(hdu, int(x), int(y), int(width), int(height)))
                #if ((hdu.header["FILTER2"].strip().lower() == "vph30" and hdu.header["FILTER1"].strip().lower() == "vph30")
                #    or (hdu.header["CAMNAME"].strip() == "HK" and hdu.header["FILTER2"].strip().lower() == "vph30")
                #    or (hdu.header["CAMNAME"].strip() == "YJ" and hdu.header["FILTER1"].strip().lower() == "vph30")):
               
            TrimmedHDUL.writeto(os.path.join(output, os.path.basename(f)), overwrite=True)



def main():
    parser = ArgumentParser()
    parser.add_argument("-r", "--root_dir", default=None, help="Directory containing files to be trimmed. If ommitted, the current working directory is used.")
    parser.add_argument("-o", "--output_dir", default=None, help="Directory in which to place the trimmed files. If omitted, the current working directory is used.")
    parser.add_argument("-W", "--width", default=None, help="Width of cropped region in pixels.")
    parser.add_argument("-H", "--height", default=None, help="Height of cropped region in pixels.")
    parser.add_argument("-x", default=None, help="Starting x-index of the cropped region in pixels.")
    parser.add_argument("-y", default=None, help="Starting y-index of the cropped region in pixels.")
    parser.add_argument("-R", "--recursive", action="store_true", default=False, help="Trim files in all sub-directories of root directory.")
    args = parser.parse_args()
    trimFiles(args.x, args.y, args.width, args.height, args.root_dir, args.output_dir, args.recursive)

if __name__ == '__main__':
    main()