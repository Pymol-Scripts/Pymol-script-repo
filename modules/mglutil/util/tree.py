#
# Last modified on Wed Oct 10 15:07:59 PDT 2001 by lindy
#
# $Header: /opt/cvs/python/packages/share1.5/mglutil/util/tree.py,v 1.3 2001/10/10 22:27:01 lindy Exp $
#


class TreeNode:
    """Base class of generic tree nodes.

    This class will work all by itself with data being
    whatever you want your nodes to be. It also could
    be used as the superclass of specific subclasses.
    """
    def __init__(self, parent=None, data=None):
        self.parent = parent
        self.children = []
        if parent:
            parent.add_child(self)
        if data:
            self.data = data


    def add_child(self, child):
        """Add given node to children of self.

        parent.add_child( node) adds node to children of parent.
        This is done in the __init__ if parent is given when node
        in instanciated (that's the prefered method).
        """
        self.children.append(child)


    def pre_traverse(self, f, *args, **kw):
        """Apply f to yourself and then your children recursively
        The function f must take the treeNode as it's first argument.
        """
        args = list(args)
        args[0] = self
        a = tuple([f] + args)

        apply(f, args, kw)
        for c in self.children:
            apply(c.pre_traverse, a, kw)


    def post_traverse(self, f, *args, **kw):
        """Traverse children with f and args, then visit parent.
        The function f must take the treeNode as it's first argument.
        """
        args = list(args)
        args[0] = self
        a = tuple([f] + args)
        for c in self.children:
            apply(c.post_traverse, a, kw)
        apply(f, args, kw)





    def get_iterator(self):
        """over-ride me to supply an appropriate subclass of TreeIterator"""
        raise NotImplementedError



class TreeIterator:
    """This iterator class is not finished yet.
    """
    def __init__(self, node):
        self.iterRoot = node
        self.currentNode = None # set by subclass
        self.done = None


    def current(self):
        """return the currently-visited node"""
        return self.currentNode


    def done(self):
        """Returns false (None) until the traversal is finished"""
        return self.done


    def first(self):
        """Reset the currentNode to the initally specified node
        """
        self.currentNode = self.iterRoot


    def next(self):
        """Move currentNode on to the next item in the traversal.

        Over-ride this method to provide a specific type of traversal
        """
        # move on to next item
        raise NotImplementedError




class PostTreeIterator(TreeIterator):
    """This iterator class is not finished yet.
    """
    def __init__(self, node):
        TreeIterator.__init__(self)
        self.nodeList = []
        self.iterRoot.post_traverse(self.nodeList.append)
        self.currentIx = 0


    def next(self):
        self.currentNode = self.nodeList[self.currentIx]
        self.currentIx = self.currentIx + 1
        if self.currentIx == len(self.nodeList):
            self.done = 1



class PreTreeIterator(TreeIterator):
    """This iterator class is not finished yet.
    """
    def __init__(self, node):
        TreeIterator.__init__(self)
        self.nodeList = []
        self.iterRoot.pre_traverse(self.nodeList.append)
        self.currentIx = 0


    def next(self):
        self.currentNode = self.nodeList[self.currentIx]
        self.currentIx = self.currentIx + 1
        if self.currentIx == len(self.nodeList):
            self.done = 1

