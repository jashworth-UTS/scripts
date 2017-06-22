# Fiji/ImageJ Python Macro script for creating composites of paired fluoresecence images (e.g. autofluorescence and YFP)
# J. Ashworth 2017 University of Technology Sydney

import ij
from ij.process import ImageConverter
from ij.plugin import RGBStackMerge
from loci.plugins import BF
from loci.plugins.in import ImporterOptions
import os,sys

path = '/Users/justin/UTS/Research/microscopy/'
filetab = '%s/ijfiles' %path

for l in open(filetab):
	f = l.strip().split()
	out = f[0]
	ch1 = '%s/%s' %(path,f[1])
	ch2 = '%s/%s' %(path,f[2])
	print('%s: opening %s and %s' %(out,ch1,ch2))
	ch1 = ij.IJ.openImage(ch1)
	conv = ImageConverter(ch1)
	conv.convertToGray16()
#	ch1.show()
	ch2 = ij.IJ.openImage(ch2)
	conv = ImageConverter(ch2)
	conv.convertToGray16()
#	ch2.show()
#	ch1.close()
#	ch2.close()
	mrg = RGBStackMerge.mergeChannels([ch1, ch2], False)
#	mrg.show()
	mrgf = '%s/%s.merged.tiff' %(path,out)
	ij.IJ.saveAs(mrg,"Tiff",mrgf)
	mrg.close()
	opts = ImporterOptions()
	opts.setId(mrgf)
	opts.setColorMode(ImporterOptions.COLOR_MODE_COMPOSITE)
	bio = BF.openImagePlus(opts)[0]
#	bio.show()
	ij.IJ.saveAs(bio.flatten(), "Tiff", '%s/%s.flat.tiff' %(path,out))
#	bio.close()
#	while(True): continue
