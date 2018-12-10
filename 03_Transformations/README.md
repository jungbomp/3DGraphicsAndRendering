# Transformation

When resterize a 3D model onto screen, it is necessary to convert 3D coordination into 2D coordination of screen. A 3D model usually has a origin coordination (0,0,0) as the center point at the bottom face. This 3D coordinate systems is model space. Our eyes or a camera capture the 3D model as 2D space, which has origin coordination is left-top or left-bottom point. Thus, when the 3D model is rendered onto a screen, first, it converts model space into world space, camera space, image space, perspective space and screen space consequently. The perspective space applies the perspective of model. For example, if a triangle in the 3D model is far away from camera focal point, it should be expressed smaller than its original size. Thus, when it applis transformation such as scaling, roatation, or transition, it is better to apply before perspective space and image space is one of good space to do this. This source code also applies transformation at image space.


## The repository includes:
* Source code
* Sample datasets
* Sample screen shots


## Datasets
### Input

The input file named rects contains rectangles information. Each line is consisted of seven numbers. The first two numbers represent top left corner cordination (X,Y) of rectangle, following two numbers mean bottom right corner cordination (x,y) of rectangle, and following three numbers are the RGB color values. The media contents consist of series of image files. Each image file is .RGB where the resolution is 352x288 containing 352x288 red bytes, followed by 352x288 green bytes, floowed by 352x288 blue bytes.

### Output

The output file is .PPM file format which has an ascii header followed by 8-bit binary pixel color values in raster order (top-left to bottom-right). For example the header is `P6 255 256 255\n` `RGBRGBRGB...` produces a 256x256 image.


## Result Screen shot
![Sample screen shot](demo.png)


## Building Environment
* Microsoft Window 10
* Microsoft Visual Studio Community 2015 Version 14.0.25431.01 Update 3


## Compile and Run
```bash
Set running configuration to Release
Build - Build Solution
> %(Solution Dir)\Release\CS580HW2.exe
Render - RunRender
```


## Status

This is the 3rd assignment of CSCI-580 3D Graphics and Rendering, 2018 fall

Version 1.0
