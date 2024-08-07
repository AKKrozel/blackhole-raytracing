# blackhole-raytracing 

This project uses raytracing to produce images of the accretion disk of a rotating black hole. Simulated light rays are projected from a virtual screen toward a virtual black hole according to the Kerr Metric and Boyer-Lindquist coordinates using an implimentation of the rk45dp integration method. The images produced are used to produce animations demonstrating what it would look like to orbit a black hole and the order in which light rays reach the black hole throughout the raytracing procedure

## Usage

Three separate .cpp files are used to various output data that is then used by the .ipynb files to create images and animations.
- raytraceOverAngles.cpp takes snapshots of the blackhole from various elevations and outputs .csv files representing the light intensity seen at each pixel on the virtual screen.
- raytraceForTimes.cpp outputs a .csv file that gives the time taken by each light ray to reach termination starting from the virtual screen.
- raytraceForLightPaths.cpp outputs .csv files containing the cartesian coordinate positions of the light rays over the course of their simulation.
- highres-bh-animation.ipynb takes .csv files containing light intensities and times taken and produces videos from the perspective of the virtual screen
- RaytraceTrajectoriesAnimation.ipynb takes .csv files containing light intensities, times, and light paths to produce a video of the ray trajectory for each pixel on the virtual screen 

 All of these files used should be fairly easy to run from the command line or a Jupyter Notebook. Be sure to allow the .cpp and .h files to access eachother by placing them in the same directory. It will also be necessary to provide appropriate file paths in highres-bh-animation.ipyb.
 
## Animations

To see the end result of the project's pipeline, check out the .mp4 files
- BH_Anim_Long.mp4, BH_Anim_Long_HD, and BH_Anim_Short.mp4 demonstrate the order in which light rays reach the black hole with their fade in and then demonstrate what it would look like to orbit the black hole. These animations are all played out from the perspective of the virtual screen
- BH_Extra_Anim.mp4 is an extra animiation that provides intuition for how the black hole affects the light ray paths by showing a side veiw of the raytracing light paths. The light paths are displayed in batches according to how long it takes them to reach the black hole.
