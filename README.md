# blackhole-raytracing 

This project uses raytracing to produce an images of the accretion disk of a rotating black hole. Simulated light rays are projected from a virtual screen toward a virtual black hole using the Kerr Metric and Boyer-Lindquist coordinates. The images produced are used to produce animations demonstrating what it would look like to orbit a black hole and the order in which light rays reach the black hole.

## Usage

This project uses three .cpp files to output data that is then used by the highres-bh-animation.ipynb to create images and animations. 
-raytraceOverAngles.cpp takes snapshots of the blackhole from various elevations and outputs .csv files representing the light intensity seen at each pixel on the virtual screen.
-raytraceForTimes.cpp outputs a .csv file that gives the time taken by each light ray to reach termination starting from the virtual screen.
-raytraceForLightPaths.cpp outputs .csv files containing the cartesian coordinate positions of the light rays over the course of their simulation

## Animations

To see the end result of the project's pipeline, check out the .mp4 files
-BH_Anim_Long.mp4, BH_Anim_Long_HD, and BH_Anim_Short.mp4 demonstrate the order in which light rays reach the black hole with their fade in and then demonstrate what it would look like to orbit the black hole. These animations are all played out from the perspective of the virtual screen
-BH_Extra_Anim.mp4 is an extra animiation that provides intuition for how the black hole affects the light ray paths by showing a side veiw of the raytracing light paths. The light paths are displayed in batches according to how long it takes them to reach the black hole.
