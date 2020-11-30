# RayTrace-Go

This is an implementation of a primitive ray tracer that I built as I was learning Go.

`go run raytrace.go` compiles and runs the program which outputs a bunch of .png images. Each image has one of the light sources in a slightly different location to generate the following .gif.

![](animated.gif)

For reference, the .gif was compiled with [ImageMagick](https://imagemagick.org/script/download.php)

`magick convert -delay 5 -loop 0 *.png animated.gif`