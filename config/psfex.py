import astromatic.psfex.psfexStarSelector

root.calibrate.measurePsf.starSelector.name = "psfex"
#root.calibrate.measurePsf.starSelector["psfex"].fluxName = "initial.flux.sinc"
root.calibrate.measurePsf.starSelector["psfex"].maxFwhmVariability = 0.1
root.calibrate.measurePsf.starSelector["psfex"].maxbadflag = False
