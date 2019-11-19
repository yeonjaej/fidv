# fidv

# Generate random vertices and check each volumes
```
root -l
.L cutflow.C
random_vtx()
```


# to test space charge boundary
literature for scb (ub DocDB: 26423)
```
root -l
.L cutflow.C
sce()
```

# to count passing events in your analysis file (change the input file and tree input as needed)
```
root -l
.L cutflow.C
cuts()
```

cuts() can take input variables,
one can also do
(5cm fid. cut, 0cm trk end cut, deltarad file?, save cutflow to pdf?)
```
root -l
.L cutflow.C
cuts(5., 0., true, false)
```



