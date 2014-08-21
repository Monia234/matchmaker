#!/usr/bin/env python

""" A collection of helper functions that I've come to develop to help with
    functional programming in Python.

    Some conventions used in this file:
        * functions of the form ``F_c'' are produced from ``curry2(F)''.
          In other words, the curried version of a function is obtained by
          suffixing the function name with ``_c''.
        * flipped versions of arbitrary functions are of the form ``F_f''.
          Some flipped versions have a more readable name, generally involving
          the word ``for''. For example, flip(map) is called ``for_each''.

    Protips:
        * If you find yourself mapping functions over the same iterable many
          times, consider making a binding with for_each_c. For example, Let
          xs be some interable:
            for_each_x = for_each_c(xs)
          now we can write:
            ys = for_each_x(x_to_y)
            zs = for_each_x(x_to_z)
            etc.
    """

from itertools import imap

def with_file(f, path, mode='r'):
    """ Run a function on a file handle, using a with-statement. Useful in lambdas. """
    with open(path, mode) as handle:
        return f(handle)

def mkdir_p(dir):
    """ Mimics the functionality of Unix `mkdir -p`, which will create a
        directory tree only if it does not exist.
        """
    if not path.exists(dir):
        os.makedirs(dir)
        return True
    else:
        return False

def is_iterable(object):
    """ A synonym-function for ``lambda x: hasatter(x, "__iter__")'', but more
        human-readable.
        """
    return hasattr(object, "__iter__")

def compose(left, right, star=False):
    """ Compose a function with another. The left function should take only one
        argument, but the right function can be of any number of arguments or
        keyword arguments. If the parameter `star` is True, then the return
        value of the right function will be exploded with a star (*). This
        allows the left function to be of any number of arguments, provided
        that the right function return an iterable of matching length.
        """
    if star:
        left_ = lambda parameter: left(*parameter)
    else:
        left_ = left

    return lambda *args, **kwargs: left_(right_(*args, **kwargs))

def curry2(function):
    """ Convert a function of two parameters into a function of the first
        parameter, yielding another function, but of the second parameter.

        Law:
            If f is a function of two parameters, then
            let g = curry2(f),
            f(x, y) = g(x)(y)
            (Provided that f and g are executed in the same environment, should
            they have side effects.)
        """
    return lambda x: lambda y: function(x, y)

def uncurry2(function):
    """ Convert a function of one parameter that yields another function of one
        parameter into a function of two parameters.

        Law:
            If f is a curried function, then
            let g = uncurry2(f),
            f(x)(y) = g(x, y)
            (Provided that f and g are executed in the same environment, should
            they have side effects.)
        """
    return lambda x, y: function(x)(y)

def curry3(function):
    """ Based on the same principle as curry2, but extended to three arguments. """
    return lambda x: lambda y: lambda z: function(f, y, z)

def uncurry3(function):
    """ Based on the same principle as uncurry2, but extended to three
        arguments.
        """
    return lambda x, y, z: function(x)(y)(z)

def flip(function):
    """ Reverse the arguments for a function of two arguments. """
    return lambda x, y: function(y, x)

def file_as_lines(handle):
    """ Produce a generator from a file handle that yields successive lines
        from the file, stripping the trailing newline.
        """
    return (line[:-1] for line in handle)

def partition(predicate, seq):
    """ Separate a sequence into two, according to the satisfaction of a given
        predicate. The returned tuples has the satisfying elements on the left,
        and the failing elements on the right.
        """
    satisfy = []
    fail = []
    for elem in seq:
        (satisfy if predicate(elem) else fail).append(elem)
    return (satisfy, fail)

def ipairs(seq):
    """ Convert an arbitrary sequence S into a list of pairs Ps, such that each
        element in Ps is of the form (S_i, S_i+1). Ps has length one less than
        S. This will effectively duplicate N-2 elements of S, where N is the
        length of S.
        Note: The result is lazy. To get the strict version, use ``pairs''.
        """
    for (i, x) in enumerate(seq):
        yield (x, seq[i+1])

def succ(n):
    """ Get the successor of a number, i.e. increment it, i.e. add one to it.
        """
    return n + 1

imap_c      = curry2(imap)
ifor_each   = flip(imap)
ifor_each_c = curry2(ifor_each)

map_c      = curry2(map)
for_each   = flip(map)
for_each_c = curry2(for_each)

with_file_c   = curry2(with_file)
with_file_f   = flip(with_file)
with_file_f_c = curry2(with_file_f)

for_file   = with_file_f
for_file_c = curry2(for_file)

