# Wurtzite symmetry

These images and animationa are rendered using the [Rayshade software](https://sourceforge.net/projects/rayshade/) software which is in the public domain and can be downlowded from the sourceforge link.
The default image format is not in the public domain which is why we use here output to an intermediary _\*.mtv_ format from which one can convert, using ImageMagick for instance, to the desire image type.
The _*.ray_ inout files can be compiled on a Linux machine after the installation of Rayshade (and making sure ImageMagick is available)  by typing:
```
rayshade [filePNG].ray > [filePNG].mtv  
convert [filePNG].mtv [file].png
```
where [file] is the name of the file used. _\*.mtv_ files can also be converted to _\*.png_ or _\*.gif_ with ImageMagick . All the output images can be found in the pngs folder and the animation in the gifs folder. To render the _\*.gif_ s one need more patience and
has to do:
```
rayshade [fileGIF].ray > [fileGIF].mtv
convert [fileGIF].mtv [file].gif
```

The input _\*.ray_ files for symmetry are closely based on those developed by Marc De Graef[^1] with the purpose of teaching group symmetry in crystallography[^2]. A link to his database can be found [here](http://som.web.cmu.edu/frames2.html) .

[^1]: M. De Graef. “_A novel way to represent the 32 crystallographic pointgroups_”. In: J. Matter: Educ 20 (1998).
[^2]: M. De Graef. [ _Teaching crystallographic and magnetic point group symmetry using three-dimensional rendered visualizations. International Union of Crystallography._ ](https://www.iucr.org/education/pamphlets/23/full-text) 2008