gerg2004-f
==========

```
The partial GERG-2004 implementation used to calculate the
compressibility factors of oxygen, helium, nitrogen and their
mixtures in the gas phase was inspired by open source code written
by David de Marneffe in December 2007. David can be contacted at
daviddemarneffe at yahoo dot com.
```

Ported to Fortran by Eskil Olsen (eskil at eskil dot org) for the heck of it.

To run, compile the .f file and run the resulting binary. The program reads
it's arguments from ```input```, which should contain 1 line ala

```
FractionO2,FractionHe,Pressure,Temperatue,Delta,System
```

where

* ```FractionO2``` is the fraction of oxygen in the gas, eg. 0.21 for air.
* ```FractionHelium``` is the fraction of helium in the gas, eg. 0.35 for 35% helium.
* ```Pressure``` the pressure of the gas in gauge reading (not absolute).
* ```Temperature``` the temperature of the gas.
* ```Delta``` the minimum delta to converge on when computing the Z value.
* ```System``` the system of measurement for pressure and temperature, use ```imperial``` or ```metric```.


Todo
* Add range checks on inputs.