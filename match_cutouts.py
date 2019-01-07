#!/usr/bin/env python3

import os
import sys

import pyfits
import math
from astLib import astWCS
import numpy
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("ref_file", #nargs=1,
                        type=str,
            help="Reference file")

    parser.add_argument("file2match", #nargs=1,
                        type=str,
            help="File to match")

    parser.add_argument("output_fn", #nargs=1,
                        type=str,
            help="Output filename")

    parser.add_argument("--keep", default=False,
                        action="store_true",
                        help="Keep intermediate frame")
    parser.add_argument("--interpol", default="LANCZOS3",
                        help="interpolateion scheme (use BILINEAR for weight maps)")

    args = parser.parse_args()


    ref_file = args.ref_file #sys.argv[1]
    file2match = args.file2match #sys.argv[2]

    output_fn = args.output_fn #sys.argv[3]


    #
    # open ref file and find center and pixelscale
    #
    ref_hdu = pyfits.open(ref_file)
    hdr = ref_hdu[0].header

    out_nx = numpy.max([hdr['CRPIX1'], hdr['NAXIS1']-hdr['CRPIX1']])
    out_ny = numpy.max([hdr['CRPIX2'], hdr['NAXIS2']-hdr['CRPIX2']])
    
    center_ra = hdr['CRVAL1'] + 0.5*hdr['CD1_1']
    center_dec = hdr['CRVAL2'] + 0.5*hdr['CD1_1']

    # center_x = hdr['NAXIS1']//2
    # center_y = hdr['NAXIS2']//2

    # center_ra = (center_x + 1 - hdr['CRPIX1']) * hdr['CD1_1'] + hdr['CRVAL1']
    # center_dec = (center_y + 1 - hdr['CRPIX2']) * hdr['CD2_2'] + hdr['CRVAL2']
    # print("ME: ", center_ra, center_dec)


    wcs = astWCS.WCS(hdr, mode='pyfits')
    # print("WCS", wcs.getCentreWCSCoords())

    pixelscale = math.fabs(hdr['CD1_1']) * 3600
    print("pixelscale: %.4f arcsec" % (pixelscale))

    # center_ra, center_dec = wcs.getCentreWCSCoords()

    memory_setup = """
    -VMEM_MAX 16636
    -MEM_MAX 4096
    -COMBINE_BUFSIZE 2048
    """

    swarp_tmp_fn = output_fn+".tmp"

    swarp_cmd = """
    swarp 
    -IMAGEOUT_NAME %(outfile)s
    -WEIGHTOUT_NAME /dev/null
    -SUBTRACT_BACK N
    -CENTER_TYPE MANUAL
    -CENTER %(center_ra).8f,%(center_dec).8f
    -PIXELSCALE_TYPE MANUAL
    -PIXEL_SCALE %(pixelscale).9f
    -IMAGE_SIZE %(nx)d,%(ny)d
    -RESAMPLING_TYPE %(interpol)s
    %(memory_setup)s
    %(input)s
    """ % dict(
        outfile=swarp_tmp_fn,
        center_ra=center_ra, center_dec=center_dec,
        pixelscale=pixelscale,
        nx=2*out_nx,
        ny=2*out_ny,
        input=file2match,
        memory_setup=memory_setup,
        interpol=args.interpol,
    )
    run_swarp = " ".join(swarp_cmd.split())
    print(run_swarp)
    os.system(run_swarp)
    

    tmp_hdu = pyfits.open(swarp_tmp_fn)
    # print(tmp_hdu[0].header)
    out_wcs = astWCS.WCS(tmp_hdu[0].header, mode='pyfits')
    
    ra_ll, dec_ll = wcs.pix2wcs(0.5, 0.5) 
    # first pixel is 1,1, not 0,0, and center of that pixel is 0.5/0.5
    _llx, _lly = out_wcs.wcs2pix(ra_ll, dec_ll)
    print(_llx, _lly)

    llx, lly = int(numpy.floor(_llx)), int(numpy.floor(_lly))
    # print(llx, lly)

    final_data = tmp_hdu[0].data[lly:lly+hdr['NAXIS2'], llx:llx+hdr['NAXIS1']]
    final_hdu = pyfits.PrimaryHDU(header=tmp_hdu[0].header, data=final_data)
    final_hdu.header['CRPIX1'] -= llx
    final_hdu.header['CRPIX2'] -= lly
    final_hdu.writeto(output_fn, clobber=True)

    print("output image size:", final_data.shape, hdr['NAXIS1'], hdr['NAXIS2'])

    if (not args.keep):
        os.remove(swarp_tmp_fn)
