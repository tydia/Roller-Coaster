
CSCI 420 Computer Graphics
Assignment 2: Roller Coaster
Tong Wang

-------------------------------------------------------------
                          How to run
-------------------------------------------------------------

This program is develeoped for Linux environment.
To compile, in terminal and in project folder, type:
	make
To clean up the executable, type:
	make clean
To run, type:
	./assign2 track.txt

-------------------------------------------------------------
                   Extra Part I Implemented
-------------------------------------------------------------

1. Draw splines using recursive subdivision.
   I used the recursive subdivision algorithm given in the hint slides
   to draw splines. One thing to note is that I'm only using the 
   subdivision method to store u's. Then using these as drawing points
   multiple times for later rendering. The order of the points are 
   from start point to end point because my subdivide() method first 
   recursively call the left half, then the right half. That's the reason
   that I am abling to move the camera directly from the beginning to end
   through the stored vector of drawing points.
   Note that my subdivide does two things:
        1) Calculate location for each drawing point and store them.
        2) Calculate tangent vector of each drawing point and store them.

2. Determine Roller Coaster's Normal Vectors Using Sloans Method
   I implemented Sloan's method mentioned in Level 4 to determine normalized
   normals of the coaster at every drawing point. These normalized normals 
   are set to be the camera's up vector. This indeed gives more realistic
   feel of actually riding on a roller coaster compare to simply set the 
   camera's up vector to be the normalized normal of the ground

3. Render double rail
   The way I achived this is, again, thanks to Sloan's methd. Given tangent
   vectors and their normals, we can calculate their binormal vectors. Then,
   for each point on spline, we get a local coordinate system defined by these
   three unit length vectors.
   Thus, we can then use -N[i] - B[i] as the direction of moving the cross 
   section to left bottom corner for drawing point i, and -N[i] + B[i] as
   the direction of moving the cross section to right bottom corner. Then 
   by rendering these two corss sections, we get the double rail. 
   The camera position never changed in process of rendering rails. It
   always moves along the original generated drawing points.

-------------------------------------------------------------
                          Other Notes
-------------------------------------------------------------

1. Making screenshots is TURNED OFF. If you want to make 
screenshots when testing, you can uncomment codes in the doIdle()
function. 
2. For your convenience, I included the pic folder for you to compile.
