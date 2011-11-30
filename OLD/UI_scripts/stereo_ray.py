from pymol import cmd
 
def stereo_ray(output='', width='', height=''):
   '''
DESCRIPTION
   "stereo_ray" ray-traces the current scene twice (separated by 
   a six-degree rotation around the y axis)
   and saves a pair of images that can be combined in any image
   manipulation software to form a stereoimage.
   The first argument, the output file name, is mandatory.
   The second and third arguments, the size of the image, are not.
   If the width is given, the height will be calculated.
 
USAGE
   stereo_ray filename [, width [, height]]
 
EXAMPLE
   stereo_ray output, 1000, 600
   stereo_ray secondImage.png
   '''
 
   if output == '':
      print 'no output filename defined\n'
      print 'try: \'stereo_ray filename\''
      return -1
      # abort if no output file name given
 
   if width == '':
      width,height = cmd.get_session()['main'][0:2]
      # get dimensions from viewer if not given
 
   elif height == '':
      oldWidth,oldHeight = cmd.get_session()['main'][0:2]
      height = int(width)*oldHeight/oldWidth
      # calculate height from proportions of viewer if
      # only width is given
 
   cmd.ray(width, height, angle=-3)
   cmd.png(output+"_r")
   cmd.ray(width, height, angle=3)
   cmd.png(output+"_l")
 
cmd.extend('stereo_ray',stereo_ray)
