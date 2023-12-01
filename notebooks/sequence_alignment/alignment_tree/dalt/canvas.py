import io

from PIL import Image

import drawsvg as draw

class Canvas:
    """
    This is the canvas drawing board.

    This class makes is easier to draw boxes and lines in a predifined grid.

    Args:
        min_col (int): Coordinate of the minimum column of the grid.
        min_row (int): Coordinate of the minimum row of the grid.
        max_col (int): Coordinate of the maximum column of the grid.
        max_row (int): Coordinate of the maximum row of the grid.
        box_width (int): Width of the boxes to be drawn (in pixels).
        box_height (int): Height of the boxes to be drawn (in pixels).
        h_margin (int): Horizontal distance between the edge of the boxes and the limit of the grid
                        cell (in pixels).
        v_margin (int): Vertical distance between the edge of the boxes and the limit of the grid
                        cell (in pixels).

    Private Attributes:
        _drawing (Drawing): The ` Drawing`  object from the ` drawsvg`  package that is used to
                            actually draw into the image
    """
    TEXT_INDENT_PERC = 0.02

    def __init__(self, min_col: int, min_row: int, max_col: int, max_row: int, box_width: int, box_height: int, h_margin: int, v_margin: int):
        self._min_row = min_row
        self._min_col = min_col

        self._box_width = box_width
        self._box_height = box_height
        self._h_margin = h_margin
        self._v_margin = v_margin

        canvas_width = (max_col - min_col + 1) * (self._box_width + self._h_margin * 2)
        canvas_height = (max_row - min_row + 1) * (self._box_height + self._v_margin * 2)

        self._drawing = draw.Drawing(canvas_width, canvas_height)
        self._drawing.append(draw.Rectangle(0, 0, '100%', '100%', rx=None, ry=None, fill='rgb(255,255,255)'))
    
    def col2x(self, col) -> int:
        """
        Converts columns coordinates into the X value of the box in pixels.

        Args:
            col (int): Column position to convert.

        Retuns:
            int: The X position of the center of the box in column ` col` , in pixels. 
        """
        return (col - self._min_col + 0.5) * (self._box_width + self._h_margin * 2)

    def row2y(self, row: int) -> int:
        """
        Converts row coordinates into the Y value of the mid box pixel coordinates.

        Args:
            row (int): Row position to convert.

        Retuns:
            int: The Y position of the center of the box in row ` row` , in pixels. 
        """
        return (row - self._min_row + 0.5) * (self._box_height + self._v_margin * 2)


    #def get_size(self) -> tuple[int, int]:
    #    """
    #    Size of the canvas.
    #
    #    Returns:
    #        tuple[int int]: Total width and height of the canvas in pixels.
    #    """
    #    return self._drawing.calc_render_size()
    
    def add_box(self, col: int, row: int, text: str=None, color: str=None, font_size: int=10, font_family: str="Monospace"):
        """
        Add a box to the specified coordinates of the grid.

        Args:
            col (int): Column in which to draw the the box.
            row (int): Row in which to draw the the box.
            text (str): Text to display in the box. Can be a multiline text.
            box_color (str): Color of the box square in the "#rrggbb" format. If ` None`  the box
                             will not be drawn
            font_size (int): Size of the font.
            font_family (str): Name of the font.
        """
        mid_x = self.col2x(col)
        mid_y = self.row2y(row)

        left_x = mid_x - (self._box_width / 2)
        top_y = mid_y - (self._box_height / 2)

        # shift left just a little so text doesn't touch the border
        h_shift = (self._box_width * self.TEXT_INDENT_PERC)

        # shift up multiline text to center it vertically in the box
        v_shift = (len(text.split("\n")) - 1) / 2 * font_size

        self._drawing.append(draw.Text(text, font_size, left_x + h_shift, mid_y - v_shift, font_family=font_family, dominant_baseline='middle'))

        if color:
            self._drawing.append(draw.Rectangle(left_x, top_y, self._box_width, self._box_height, stroke_width=2, stroke=f"{color}", fill='none'))

    def add_link(self, start_col: int, start_row: int, end_col: int, end_row: int, width: int=1, color: str="black", opacity: float=0.2):
        """
        Add a line between two boxes.

        The line will connect right edge of the "start" box with the left edge of the "end" box at
        mid height.

        Args:
            start_col (int): Column coordinate of the start box.
            start_row (int): Row coordinate of the start box.
            end_col (int): Column coordinate of the end box.
            end_row (int): Row coordinate of the end box.
            width (int): Width of the line to be drawn in pixels.
            color (str): Color of the line in the "#rrggbb" format.
            opacity (float): Opacity level of the line from 0 (transparent) to 1 (solid).
        """
        x0 = self.col2x(start_col) + self._box_width / 2
        y0 = self.row2y(start_row)
        
        x1 = self.col2x(end_col) - self._box_width / 2 
        y1 = self.row2y(end_row)
        
        self._drawing.append(draw.Line(x0, y0, x1, y1, stroke_width=width, stroke=color, fill=color, fill_opacity=opacity))
        
    def image(self) -> draw.Drawing:
        """
        Returns a `PIL.Image` object that can be used to visualize and save the image.

        Returs:
            draw.Drawing: Drawing object.
        """
        return Image.open(io.BytesIO(self._drawing.rasterize().png_data))
