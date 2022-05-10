[![Actively Maintained](https://img.shields.io/badge/Maintenance%20Level-Actively%20Maintained-green.svg)](https://gist.github.com/cheerfulstoic/d107229326a01ff0f333a1d3476e068d)


JWcs
======================

JWcs is a pure JAVA library, that implements the "World Coordinate System" (WCS) standard in FITS (Flexible Image Transport System). 

## Description
This library allows to convert pixel <--> world coordinates. In addition, it handles sky system conversion.

## Requirements
JAVA >= 8

## Installation

### Getting the sources

```console
$ git clone https://github.com/esdc-esac-esa-int/JWcs.git
```
	
### Building the sources

Build the sources

```console
$ cd JWcs
$ mvn install
```
	
### Building the documentation

Build the documentation

```console
$ mvn site
```
	
The documentation will be generated in the directory target/site/


### Utility tool

```console
	$ java -jar target/JWcs-1.2.2.jar -h
	Usage: java -jar JWcs.jar -g PROG [OPTIONS]
	    or java -jar JWcs.jar --file HDR_FILE --project X,Y [OPTIONS]
	    or java -jar JWcs.jar --file HDR_FILE --unproject RA,DEC [OPTIONS]
	    or java -jar JWcs.jar --file HDR_FILE --convert RA,DEC --to SYS_TARGET [OPTIONS]
	    or java -jar JWcs.jar --convert RA,DEC --from SYS_ORGIN --to SYS_TARGET [OPTIONS]
	           where:
	               - PROG: either projection or converter
	               - HDR_FILE: Header FITS or FITS file
	               - X: pixel coordinate along X axis on the camera (starts to 1) 
	               - Y: pixel coordinate along Y axis on the camera (starts to 1) 
	               - RA: sky coordinate
	               - DEC: sky coordinate
	               - SYS_ORIGIN: sky system of the sky coordinates
	               - SYS_TARGET: convert sky coordinates to the SYS_TARGET
	
	           SYS_ORIGIN or SYS_TARGET can be:
	               - GALACTIC
	               - SUPER_GALACTIC
	               - EQUATORIAL
	               - EQUATORIAL(ICRS())
	               - EQUATORIAL(J2000())
	               - EQUATORIAL(FK5())
	               - EQUATORIAL(FK5(<equinox>))
	               - EQUATORIAL(FK4())
	               - EQUATORIAL(FK4(<equinox>))
	               - EQUATORIAL(FK4(<equinox>,<epoch>))
	               - EQUATORIAL(FK4_NO_E())
	               - EQUATORIAL(FK4_NO_E(<equinox>))
	               - EQUATORIAL(FK4_NO_E(<equinox>,<epoch>))
	               - ECLIPTIC
	               - ECLIPTIC(ICRS())
	               - ECLIPTIC(J2000())
	               - ECLIPTIC(FK5())
	               - ECLIPTIC(FK5(<equinox>))
	               - ECLIPTIC(FK4())
	               - ECLIPTIC(FK4(<equinox>))
	               - ECLIPTIC(FK4(<equinox>,<epoch>))
	               - ECLIPTIC(FK4_NO_E())
	               - ECLIPTIC(FK4_NO_E(<equinox>))
	               - ECLIPTIC(FK4_NO_E(<equinox>,<epoch>))
	
	Projection and sky conversion library
	
	Mandatory arguments to long options are mandatory for short options too.
	  -p, --project            Project a pixel to the sky
	  -u, --unproject          Unproject a point on the sky to 2D 
	  -f, --file               Header file or Fits file starting by a scheme (ex: file://, http://)
	  -s, --from               Origin sky system
	  -t, --to                 Target sky system
	  -c, --convert            Convert a sky coordinate from a sky system to antoher one
	  -g, --gui                Display projection or converter with a GUI
	  -h, --help               Display this help and exit
	
	OPTIONS are the following:
	  -d, --debug              Sets the DEBUG level : ALL,CONFIG,FINER,FINEST,INFO,OFF,SEVERE,WARNING
	  -e, --extension          HDU number starting at 0 when --file argument is used. If not set, 0 is default
	  -r, --precision          Precision such as %.6f. By default, precision is set to %.15f


