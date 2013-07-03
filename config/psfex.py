import astromatic.psfex.psfexStarSelector
import astromatic.psfex.psfexPsfDeterminer

root.calibrate.measurePsf.starSelector.name = "psfex"
root.calibrate.measurePsf.psfDeterminer.name = "psfex"

root.calibrate.measurePsf.starSelector["psfex"].maxFwhmVariability = 0.1
root.calibrate.measurePsf.starSelector["psfex"].maxbadflag = False
