import grass.script as grass

print(grass)

## v.extract sa_quaternaries out=liebenbergsvlei where="TERTIARY = 'C83'"
## v.dissolve liebenbergsvlei output=lieb_dissolve column=TERTIARY
## v.to.rast input=lieb_dissolve output=lieb_mask use=val
