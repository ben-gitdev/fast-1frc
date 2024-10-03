# fast-1frc
Implementation of fast single image Fourier ring correlation based on the paper
Single image Fourier ring correlation by Bernd Rieger et. al.
https://opg.optica.org/oe/fulltext.cfm?uri=oe-32-12-21767&id=551322
with massive speed improvement.
Takes 25 ms to process one 16 bit 4608x2592 pixels gray scale image
on a Intel Xeon Gold 6226R processor.
Speed independent of image intensity.
