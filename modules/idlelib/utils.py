def callback(func, *myargs):
    def w(*args):
        return func(*(args + myargs))
    return w
