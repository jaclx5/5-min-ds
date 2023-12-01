import math
import os
from dataclasses import dataclass

from PIL import Image
from tqdm import tqdm

from .alignment import Alignment
from .algorithm import Algorithm

COLOR_EXPANDED_BOX = "#00ff00"
COLOR_SOLUTION_BOX = "#0000ff"

BOX_WIDTH = 80
BOX_HEIGHT = 35
H_MARGIN = 2
V_MARGIN = 5

@dataclass
class MovieFrame:
    img: Image
    root_x: int
    root_y: int
    end: bool
    steps: int

    def save(self, file_name):
        self.img.save(file_name)

count = 0

class Movie:
    def __init__(self, frames=None):
        self._frames = frames if frames else []
        self._centered = False

    def center_frames(self):
        def add_margin(img, x, y, width, height):
            """
            from: https://note.nkmk.me/en/python-pillow-add-margin-expand-canvas/
            """
            new_img = Image.new(img.mode, (width, height), (255, 255, 255))
            new_img.paste(img, (x, y))
            
            return new_img

        if self._centered:
            # no need to recenter
            return
            
        self._centered = True
        
        # resizes and aligns all images so they have the same size and the start node appears
        # allways in the same position.
        min_x = min(map(lambda frame: -frame.root_x, self._frames))
        min_y = min(map(lambda frame: -frame.root_y, self._frames))
        max_x = max(map(lambda frame: frame.img.width - frame.root_x, self._frames))
        max_y = max(map(lambda frame: frame.img.height - frame.root_y, self._frames))

        width = max_x - min_x
        height = max_y - min_y

        for frame in self._frames:
            # shifts the coordinates to align with the start node
            x_shift = -frame.root_x - min_x
            y_shift = -frame.root_y - min_y

            # adds a white margin to center the image and make tham all the same size
            frame.img = add_margin(frame.img, int(x_shift), int(y_shift), int(width), int(height))

    def add_frame(self, frame):
        self._centered = False
        self._frames.append(frame)

    def get_frames(self, i):
        if i < len(self._frames):
            return self._frames[i]
        else:
            assert False, f"Frame index {i} must be lower than {len(self._frames)}."

    def frame_count(self):
        return len(self._frames)

    def save(self, image_path, image_name="step_$STEP$.png"):
        self.center_frames()

        for (i, frame) in enumerate(self._frames):
            # save the resized image
            fname = os.path.join(image_path, image_name.replace("$STEP$", f"{i:02d}"))

            frame.save(fname)


class Simulation:
    """
    This classes simulate algorithms and generate frames and movies of specific states.
    """
    def __init__(self, aln: Alignment, algo: Algorithm):
        self._aln = aln
        self._algo = algo

        self._count_steps = None

    def get_steps(self):
        return self._count_steps

    def frame(self, max_steps):
        # run at most `max_steps` from the algorihtm
        end, steps = self._algo.run(self._aln, max_steps=max_steps)

        # generates the image
        img = self._aln.draw(BOX_WIDTH, BOX_HEIGHT, H_MARGIN, V_MARGIN)

        # get the coordinates of the starting node (must be called after draw)
        x, y = self._aln.get_xy()

        return MovieFrame(img, x, y, end, steps)
        
    def movie(self, max_steps, progress=False):
        movie = Movie()
        
        # generate the individual steps of the algorithm
        for i in tqdm(range(max_steps), disable=not progress):

            # run the algorithm up to the i-th step and keep the frame
            frame = self.frame(max_steps=i)

            movie.add_frame(frame)

            if frame.end:
                self._count_steps = i
                break

        movie.center_frames()
        
        return movie