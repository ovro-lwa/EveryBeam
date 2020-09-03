from astropy.io import fits

with fits.open("template-image.fits") as img:

    N = img[0].data.shape[-1]

    offset = N//4

    img[0].data[:] = 0.0
    img[0].data[0,0,offset,offset] = 1.0
    img[0].data[0,0,offset,N-offset] = 1.0
    img[0].data[0,0,N-offset,offset] = 1.0
    img[0].data[0,0,N-offset,N-offset] = 1.0
    img.writeto('wsclean-I-model.fits', overwrite = True)

    img[0].data[:] = 0.0
    img[0].data[0,0,offset,N-offset] = 1.0
    img.writeto('wsclean-Q-model.fits', overwrite = True)

    img[0].data[:] = 0.0
    img[0].data[0,0,N-offset,offset] = 1.0
    img.writeto('wsclean-U-model.fits', overwrite = True)

    img[0].data[:] = 0.0
    img[0].data[0,0,N-offset,N-offset] = 1.0
    img.writeto('wsclean-V-model.fits', overwrite = True)

