'''
http://pymolwiki.org/index.php/movie_fade

(c) 2011 Jason Vertrees
(c) 2013 Thomas Holder

License: BSD-2-Clause
'''

from pymol import cmd, CmdException


def movie_fade(setting, startFrame, startVal, endFrame, endVal=None, selection=""):
    """
DESCRIPTION

    Fades representations in movies with their transparency settings.

USAGE

    movie_fade setting, startFrame, startVal, endFrame, endVal [, selection ]

EXAMPLE

    fetch 1rx1, async=0
    as cartoon
    show surface
    mset 1x80
    movie.roll
    movie_fade transparency,  1, 0., 40, 1.
    movie_fade transparency, 41, 1., 80, 0.

SEE ALSO

    mdo, mappend, set
    """
    startFrame, endFrame, startVal = int(startFrame), int(endFrame), float(startVal)
    endVal = abs(1.0 - startVal) if endVal is None else float(endVal)

    if startFrame == endFrame:
        raise CmdException("start == end")

    if startFrame > endFrame:
        startFrame, endFrame = endFrame, startFrame
        startVal, endVal = endVal, startVal

    for frame in range(startFrame, endFrame + 1):
        frac = float(frame - startFrame) / (endFrame - startFrame)

        value = (1.0 - frac) * startVal + frac * endVal
        cmd.mappend(frame, "/cmd.set(%s, %f, %s)" % (repr(setting), value, repr(selection)))

cmd.extend("movie_fade", movie_fade)
cmd.auto_arg[0]["movie_fade"] = cmd.auto_arg[0]["set"]
