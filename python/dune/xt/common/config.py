from importlib import import_module
import sys


def _can_import(module):
    def _can_import_single(m):
        try:
            import_module(m)
            return True
        except ImportError:
            pass
        return False
    if not isinstance(module, (list, tuple)):
        module = [module]
    return all((_can_import_single(m) for m in module))


_PACKAGES = {
    'K3D': lambda: _can_import('k3d'),
}


# copied from https://github.com/pymor/pymor/blob/main/src/pymor/core/config.py and adapted
class Config:

    def __init__(self):
        self.PYTHON_VERSION = f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}'

    @property
    def version(self):
        from dune.xt import __version__
        return __version__

    def __getattr__(self, name):
        if name.startswith('HAVE_'):
            package = name[len('HAVE_'):]
        elif name.endswith('_VERSION'):
            package = name[:-len('_VERSION')]
        else:
            raise AttributeError

        if package in _PACKAGES:
            try:
                version = _PACKAGES[package]()
            except ImportError:
                version = False

            if version is not None and version is not False:
                setattr(self, 'HAVE_' + package, True)
                setattr(self, package + '_VERSION', version)
            else:
                setattr(self, 'HAVE_' + package, False)
                setattr(self, package + '_VERSION', None)
        else:
            raise AttributeError

        return getattr(self, name)

    def __dir__(self, old=False):
        keys = set(super().__dir__())
        keys.update('HAVE_' + package for package in _PACKAGES)
        keys.update(package + '_VERSION' for package in _PACKAGES)
        return list(keys)

    def __repr__(self):
        status = {p: (lambda v: 'missing' if not v else 'present' if v is True else v)(getattr(self, p + '_VERSION'))
                  for p in _PACKAGES}
        key_width = max(len(p) for p in _PACKAGES) + 2
        package_info = [f"{p+':':{key_width}} {v}" for p, v in sorted(status.items())]
        separator = '-' * max(map(len, package_info))
        package_info = '\n'.join(package_info)
        info = f'''
dune-xt Version {self.version}

Python: {self.PYTHON_VERSION}

External Packages
{separator}
{package_info}
'''[1:]
        return info


config = Config()
