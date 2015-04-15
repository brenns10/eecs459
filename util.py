"""Utility functions that we may use."""

import sys
import inspect


def _convert_args(l):
    """Take a list of args, find kwargs, and return an arg,kwarg list."""
    args = []
    kwargs = dict()
    for arg in l:
        if arg.count('=') > 0:
            name, arg = arg.split('=', maxsplit=1)
            if name != '':
                kwargs[name] = arg
            else:
                args.append(arg)
        else:
            args.append(arg)
    return args, kwargs


def quick_main(module='__main__'):
    """A quick main function for a module.

    This is a quick-n-dirty main function that allows you to pick the function
    and arguments to run.  Running this will parse the command line arguments.
    The first argument must be the name of a function defined in the given
    module.  The remaining args are the arguments to that function.  You can
    pass keyword arguments as 'name=value'.  This calls the named function with
    the arguments given, passed as strings.  If no matching function is found,
    it displays an error message to stderr.

    """
    # Get all (name, function) bindings in module.
    functions = inspect.getmembers(sys.modules[module], inspect.isfunction)

    # Filter to those that are actually defined in the module (rather than
    # imported).
    functions = [x for x in functions if x[1].__module__ == module]

    # Execute the function they specified.
    for name, function in functions:
        if sys.argv[1] == name and function.__module__ == module:
            args, kwargs = _convert_args(sys.argv[2:])
            function(*args, **kwargs)
            return

    # Error message for invalid function.
    print('Invalid function "%s".' % sys.argv[1], file=sys.stderr)
    print('Choices: [%s].' % ', '.join([n for n, f in functions]),
          file=sys.stderr)
