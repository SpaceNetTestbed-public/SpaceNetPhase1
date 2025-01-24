import glob
from PIL import Image


def convert_gif(frame_folder, output_path_with_name, delay):
    frames = [Image.open(image) for image in sorted(glob.glob(f"{frame_folder}/*.jpg"))]
    frame_one = frames[0]
    frame_one.save(output_path_with_name+".gif", format="GIF", append_images=frames,
               save_all=True, duration=delay, loop=0)
