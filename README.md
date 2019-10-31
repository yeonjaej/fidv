# fidv

# to test space charge boundary
```
root -l
.L cutflow.C
sce()
```

# to test fiducial volume 
```
root -l
.L cutflow.C
cuts()
```

cuts() can take input variables,
one can also do
```
root -l
.L cutflow.C
cuts(5., 0., true, false)
```
: 5cm fid. cut, 0cm trk end cut, deltarad file, save cutflow to pdf

# to throw random vertices and count them
```
root -l
.L cutflow.C
random_vtx()
```
