# _Z/p^s-Additive Codes_ Package

_Z/p^s-Additive Codes_ is a package that generalizes the functionality
for linear codes over Z4 in Magma to linear codes over Z/p^s with p prime. 
Specifically, there are functions to construct some families of these 
codes, functions related to generalized Gray maps, information sets, 
the process of encoding, and decoding using permutation decoding, among 
others. 

The _Z/p^s-Additive Codes_ package consists of files written in the 
Magma language.  Please, send your  bug reports to Combinatorics, 
Coding and  Security Group (CCSG)  at [support.ccsg@uab.cat](mailto:support.ccsg@uab.cat)
or if it is a Magma problem to magma-trouble (magma@maths.usyd.edu.au). See section [Bug reports](#bug-reports) below.

The "Z/p^s-Additive Codes" package has been written by several 
assistant developers (supervised by Mercè Villanueva)
as a support for some research projects on Z/p^s-linear codes developed
by the Combinatorics, Coding and  Security Group (CCSG)  within the 
Department of Information and  Communications Engineering (dEIC) at 
the Autonomous University of Barcelona (UAB).

This version of the package has been developed in Magma version V2.27-5.


## Composition of the package

The _Z/p^s-Additive Codes_ package is composed of five directories:

* [/src](src): The files to be attached to Magma: 
      ZpAdditiveCodes_Core.m, ZpAdditiveCodes_Constructions.m,
      ZpAdditiveCodes_Distances.m, ZpAdditiveCodes_Encoding.m and
      ZpAdditiveCodes_Decoding.m
* [/license](license): The license of the package.
* [/doc](doc): The manual of the package, in pdf format.
* [/test](test): Several files to test the package.
       They can be loaded in Magma as soon as
       the package is attached.
* [/examples](examples): Examples from the manual. They can be loaded in
           Magma as soon as the package is attached.			



## Using/Installing "Z/p^s-Additive Codes"

To use  "Z/p^s-Additive Codes"  temporally  (as a Magma Package)
unpack  the  archive  file in a directory.   Enter to the ./src
directory. Call Magma and then write:
```
   AttachSpec(“ZpAdditiveCodes.spec”);
```
in the Magma command line.

To install "Z/p^s-Additive Codes" permanently (as a Magma Package) on Linux:

1. Unpack the archive file in a directory <code>$DIR</code>.

2. If you do not have a directory to store user-defined packages, create one in your preferred location, for instance <code>$HOME</code>:

   ```
      mkdir UserPackages
   ```

   If you already have such a directory, proceed with the instructions changing <code>$HOME/UserPackages</code> by the path to your directory.

3. Create a new directory in the <code>$HOME/UserPackages</code> directory:

   ```
      cd UserPackages
      mkdir ZpAdditiveCodes
   ```

4. Copy all files in <code>$DIR/src/</code> into this new directory <code>ZpAdditiveCodes</code>:

   ```
      cp $DIR/src/* $HOME/UserPackages/ZpAdditiveCodes/
   ```

5. Create a file named <code>spec</code> in the directory <code>$HOME/UserPackages</code>:

   ```
      touch spec
   ```

   Edit the <code>spec</code> file and add the following content:

   ```
      ZpAdditiveCodes
      {
         +ZpAdditiveCodes.spec
      }
   ```

   In case that the <code>spec</code> file already exists, add the lines above at the end of the old <code>spec</code> file.

6. Ensure that all files have the correct permissions:

   ```
      chmod -R a+rX .
   ```

7. Set the environment variable <code>MAGMA_USER_SPEC</code> to the <code>spec</code> file. Change to the directory where Magma is installed and edit the <code>magma</code> script. Locate the line

   ```
      export MAGMA_SYSTEM_SPEC
   ```

   and add the following lines just after that:

   ```
      MAGMA_USER_SPEC="$HOME/UserPackages/spec"
      export MAGMA_USER_SPEC
   ```

To do the installation on Windows OS, follow the above items from 1 to 5 (in item 4, use "copy" instead of "cp"; and in item 5, use "type nul >>" instead of "touch" if a spec file does not exist). Finally, edit the system environment variables by adding MAGMA_USER_SPEC and set it to the spec file.

In order to check that the package has been installed correctly, run Magma in a terminal window and try to run the following lines:

```
   ZpAdditiveCodes_Core_version();
   ZpAdditiveCodes_Constructions_version();
   ZpAdditiveCodes_Distances_version();
   ZpAdditiveCodes_Encoding_version();
   ZpAdditiveCodes_Decoding_version();
```

If the installation has been successful, Magma should return the following lines, one for each function, respectively:

```
   [2,1]
   [2,0]
   [1,0]
   [2,0]
   [2,2]
```

If the numbers appear but are different from the ones shown above, then the respective files have not been installed correctly and they may correspond to a previous version.



## Bug reports

When  sending a  bug  report  to support-ccsg@deic.uab.cat or to
magma@maths.usyd.edu.au,    remember  we will need to be able to
reproduce the problem; so please include:

 * The  version  of  Magma  you  are  using; either look  at the
   header when you start up Magma.
 * The  operating  system you are using e.g. Linux,  SunOS 5.8 =
   Solaris 2.8, IRIX 6.5, Windows, ...
 * A script that demonstrates the bug, along  with a description
   of why it's a bug (e.g.  by  adding  comments to  the  script - 
   recall  comments  in Magma  begin  with  a  //  or  between
   /*  */).


February 06, 2024
