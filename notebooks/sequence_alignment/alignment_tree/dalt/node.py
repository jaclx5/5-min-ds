from dataclasses import dataclass

from PIL import Image

from .canvas import Canvas

@dataclass
class TreeBoxCoords:
    min_col: int
    min_row: int
    max_col: int
    max_row: int


class Node:
    """
    This class represents a node of a tree that can be drawn in the canvas.

    It implements a recursive positioning algorithm that assigns a column and a row to each node.

    Args:
        text (str): Text to be shown in the tree.
        color (str): Color of the box.

    Private Attributes:
        _children (list[Node]): List of the children of the current node.
        _col (int): Column to be assigned to the `Node` before drwaing the full tree.
        _row (int): Row to be assigned to the `Node` before drwaing the full tree.
    """
    def __init__(self, text: str="", color: str=None):
        self._text = text
        self._color = color

        # by default a node is a root when it's created
        self._id = "*"


        self.reset()

    def reset(self):
        self._children = []
        self._col, self._row = None, None
        self._xy = None
        self._color = None

    def _get_color(self):
        return self._color

    def _set_color(self, color):
        self._color = color

    color = property(fget=_get_color, fset=_set_color)

    def _get_text(self):
        return self._text

    def _set_text(self, text):
        self._text = text
        
    text = property(fget=_get_text, fset=_set_text)

    def _get_id(self):
        return self._id

    id = property(fget=_get_id)

    def add_child(self, child: type["Node"]):
        """
        Add the `child` node as a child of the current node.

        As any tree this should be a Direct Acyclic Graph. The code however does not check for it.
        It's the responsibility of the user to avoid cycles.

        Args:
            child (Node): Node to be added as a child.
        """
        self._children.append(child)

        # set the id of the child based on it's own id
        child._id = f"{self._id}.{len(self._children)}"
        
        return self
            
    def is_leaf(self) -> bool:
        """
        Returns True if the current node is a leaf node, i.e., does not contain any children.
        """
        return not bool(self._children)
       

    def _rearrange_leaf_children(self, first_leaf: int, count_leaf: int, start_non_leaf_row: int, end_non_leaf_row: int):
        """
        Compute the ideal position of the leaf nodes (LNs) between non leaf nodes (NLNs). The idea
        is that LNs should be distributed uniformly between the surrounding NLNs.

        For example: if we have a start NLN in row 2 followed by 2 LNs ending in a NLN in row 8 the
        LNs should occupy rows 4 and 6 to leave 1 space row between all of them.

        Args:
            first_leaf (int): Index of the first LN in the list of children.
            count_leaf (int): Number of LNs in the internal.
            start_non_leaf_row (int): Row coordinate of the starting NLN.
            end_non_leaf_row (int): Row coordinate of the end NLN.
        """
        if start_non_leaf_row is None:
            start_non_leaf_row = end_non_leaf_row - count_leaf - 1

        elif end_non_leaf_row is None:
            end_non_leaf_row = start_non_leaf_row + count_leaf + 1

        step = int((end_non_leaf_row - start_non_leaf_row - 1) / (count_leaf + 1)) + 1

        for i in range(count_leaf):
            self._children[first_leaf + i]._row = start_non_leaf_row + (step * i) + 1
            
    def _arrange_all(self, level: int=0, _next_free_row: list[int]=None, _min_free_row: int=0):
        """
        Computes the (row, col) coordinates of the grid in which the node will be drawn. This
        coordinates will be later converted to absolute (x, y) canvas coordinates.
        """
        _next_free_row = _next_free_row if _next_free_row else []
        
        if len(_next_free_row) <= level:
            # initialize the next free row for this depth level
            _next_free_row.append(0)

        # for now the row value of the node is the next free on the respective column
        self._col = level
        self._row = max(_next_free_row[level], _min_free_row)

        # the tree box starts with it's own coordinates
        self._box = TreeBoxCoords(self._col, self._row, self._col, self._row)

        if not self.is_leaf():
            # get the index of the middle children to help balance the tree
            mid_child_ndx = int(len(self._children) / 2)
            
            for i, child in enumerate(self._children):
                child_col, child_row, child_box = child._arrange_all(level + 1, _next_free_row, self._row - mid_child_ndx)

                # adjusts the own node row to align with the middle children
                #if mid_child_ndx == i:
                #    self._row = max(self._row, child_row)

                # update the tree box coordinates given the new children tree box.
                self._box = TreeBoxCoords(
                    min(self._box.min_col, child_box.min_col),
                    min(self._box.min_row, child_box.min_row),
                    max(self._box.max_col, child_box.max_col),
                    max(self._box.max_row, child_box.max_row)
                )

            """
            Rebalance the leaf nodes.

            This is a complicated way to guarantee that leaf nodes are not too close to their upper
            siblings.
            """
            first_leaf, count_leaf, start_non_leaf_row = None, 0, None

            for i, child in enumerate(self._children):
                if child.is_leaf():
                    first_leaf = i if first_leaf is None else first_leaf
                    count_leaf += 1
                else:
                    if count_leaf:
                        self._rearrange_leaf_children(first_leaf, count_leaf, start_non_leaf_row, self._children[i]._row)
                        
                    # reinitialize
                    first_leaf, count_leaf, start_non_leaf_row = None, 0, self._children[i]._row

            if count_leaf and start_non_leaf_row:
                self._rearrange_leaf_children(first_leaf, count_leaf, start_non_leaf_row, None)

            # re-adjusts the own node row to align with the middle children (if changed)
            for i, child in enumerate(self._children):
                if mid_child_ndx == i:
                    self._row = max(self._row, child._row)
            # -----------------------------------------

        # update the free row for the next node in this column.
        _next_free_row[level] = self._row + 1

        return self._col, self._row, self._box

    def _draw_all(self, canvas):           
        """
        TODO: Document
        """
        canvas.add_box(self._col, self._row, self.text, self.color)

        for child in self._children:
            canvas.add_link(self._col, self._row, child._col, child._row)
            child._draw_all(canvas)

        self._xy = canvas.col2x(self._col), canvas.row2y(self._row)

    def draw(self, box_width: int, box_height: int, h_margin: int, v_margin: int) -> Image:
        """
        This method returns an image with a tree rooted in the current node. Additionally it can
        save the image into a file.

        When calling this method it is assumed that the current node is the root node. This way we
        can draw any subtree as a isolated tree on is own.

        Args:
            box_width, box_height (int): Dimensions of the box to be drawn for each node, in pixels.
            h_margin, v_margin (int): Space between the boxes, in pixels.
            image_file (str): Name of the file where o save the tree.

        Returns:
            A `PIL.Image` object with the image.

        """

        # positions all nodes in a grid
        _, _, box = self._arrange_all()
        
        # initializes the drawing canvas
        self._canvas = Canvas(box.min_col, box.min_row, box.max_col, box.max_row, box_width, box_height, h_margin, v_margin)

        # adds the node (and children) to the canvas
        self._draw_all(self._canvas)

        # generates the image
        return self._canvas.image()

    def get_xy(self):
        """
        """


        assert self._xy is not None, "Absolute XY coordinates not available yet.\nRun `draw` method before calling this function."
        
        return self._xy

    def _draw_text_all(self, _level: int=0):
        """
        Prints the information for each node and, recursively, traverse all children.

        Args:
            _level (int): Depth of recursion. Used to ident the text.
        """

        txt = self.text.replace('\n', '<nl>')
        color = f"[{self.color}]" if self.color else ""

        print(f"{self._id}: {'  ' * _level}{txt} ({self._col} x {self._row}) {color}")
        
        for child in self._children:
            child._draw_text_all(_level + 1)

    def draw_text(self):
        """
        Prints a textual representation of the tree. Assuming the current node as the tree's root.
        """

        # positions all nodes in a grid
        col, row, box = self._arrange_all()

        print("------------------------")
        print(f"Root node: {col} x {row}")
        print(f"Box limits: [{box.min_col} x {box.min_row}] - [{box.max_col} x {box.max_row}]")
        print("------------------------")

        # call the recursive function that will actually print
        self._draw_text_all()

    def _count_children_all(self, root):
        # not to count the node calling the function
        count = 0 if root else 1

        for child in self._children:
            count += child._count_children_all(False)

        return count

    def count_children(self):
        return self._count_children_all(True)

    def _get_by_level_all(self, level):
        nodes = []

        if level == 0:
            nodes.append(self)
        else:
            nodes = []

            for child in self._children:
                nodes += child._get_by_level_all(level - 1)

        return nodes

    def get_by_level(self, level):
        return self._get_by_level_all(level - 1)
