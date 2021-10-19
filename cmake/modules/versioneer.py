# Version: 0.19

"""The Versioneer - like a rocketeer, but for versions.

The Versioneer
==============

* like a rocketeer, but for versions!
* https://github.com/python-versioneer/python-versioneer
* Brian Warner
* License: Public Domain
* Compatible with: Python 3.6, 3.7, 3.8, 3.9 and pypy3
* [![Latest Version][pypi-image]][pypi-url]
* [![Build Status][travis-image]][travis-url]

This is a tool for managing a recorded version number in distutils-based
python projects. The goal is to remove the tedious and error-prone "update
the embedded version string" step from your release process. Making a new
release should be as easy as recording a new tag in your version-control
system, and maybe making new tarballs.

## License

To make Versioneer easier to embed, all its code is dedicated to the public
domain. The `_version.py` that it creates is also in the public domain.
Specifically, both are released under the Creative Commons "Public Domain
Dedication" license (CC0-1.0), as described in
https://creativecommons.org/publicdomain/zero/1.0/ .

[pypi-image]: https://img.shields.io/pypi/v/versioneer.svg
[pypi-url]: https://pypi.python.org/pypi/versioneer/
[travis-image]:
https://img.shields.io/travis/com/python-versioneer/python-versioneer.svg
[travis-url]: https://travis-ci.com/github/python-versioneer/python-versioneer

"""

import configparser
import errno
import json
import os
import re
import subprocess
import sys


class VersioneerConfig:
    """Container for Versioneer configuration parameters."""

PROJECT_ROOT = None

def get_root():
    """Get the project root directory.

    We require that all commands are run from the project root, i.e. the
    directory that contains setup.py, setup.cfg, and versioneer.py .
    """
    root = os.path.realpath(os.path.abspath(PROJECT_ROOT))

    if not (os.path.exists(root) or os.path.isdir(root)):
        err = (
            "Versioneer was unable to run the project root directory. "
        )
        raise VersioneerBadRootError(err)
    return root


def get_config_from_root(root):
    """Read the project setup.cfg file to determine Versioneer config."""

    CONFIG = """
    # Hardcoded config for all our modules

    [versioneer]
    VCS = git
    style = pep440
    versionfile_source =
    versionfile_build =
    tag_prefix =
    parentdir_prefix =
    """

    # This might raise EnvironmentError (if setup.cfg is missing), or
    # configparser.NoSectionError (if it lacks a [versioneer] section), or
    # configparser.NoOptionError (if it lacks "VCS="). See the docstring at
    # the top of versioneer.py for instructions on writing your setup.cfg .
    parser = configparser.ConfigParser()
    parser.read_string(CONFIG)
    VCS = parser.get("versioneer", "VCS")  # mandatory

    def get(parser, name):
        if parser.has_option("versioneer", name):
            return parser.get("versioneer", name)
        return None

    cfg = VersioneerConfig()
    cfg.VCS = VCS
    cfg.style = get(parser, "style") or ""
    cfg.versionfile_source = get(parser, "versionfile_source")
    cfg.versionfile_build = get(parser, "versionfile_build")
    cfg.tag_prefix = get(parser, "tag_prefix")
    if cfg.tag_prefix in ("''", '""'):
        cfg.tag_prefix = ""
    cfg.parentdir_prefix = get(parser, "parentdir_prefix")
    cfg.verbose = get(parser, "verbose")
    return cfg


class NotThisMethod(Exception):
    """Exception raised if a method is not valid for the current scenario."""


# these dictionaries contain VCS-specific tools
LONG_VERSION_PY = {}
HANDLERS = {}


def register_vcs_handler(vcs, method):  # decorator
    """Create decorator to mark a method as the handler of a VCS."""

    def decorate(f):
        """Store f in HANDLERS[vcs][method]."""
        if vcs not in HANDLERS:
            HANDLERS[vcs] = {}
        HANDLERS[vcs][method] = f
        return f

    return decorate


def run_command(commands, args, cwd=None, verbose=False, hide_stderr=False, env=None):
    """Call the given command(s)."""
    assert isinstance(commands, list)
    p = None
    for c in commands:
        try:
            dispcmd = str([c] + args)
            # remember shell=False, so use git.cmd on windows, not just git
            p = subprocess.Popen(
                [c] + args, cwd=cwd, env=env, stdout=subprocess.PIPE, stderr=(subprocess.PIPE if hide_stderr else None)
            )
            break
        except EnvironmentError:
            e = sys.exc_info()[1]
            if e.errno == errno.ENOENT:
                continue
            if verbose:
                print("unable to run %s" % dispcmd)
                print(e)
            return None, None
    else:
        if verbose:
            print("unable to find command, tried %s" % (commands,))
        return None, None
    stdout = p.communicate()[0].strip().decode()
    if p.returncode != 0:
        if verbose:
            print("unable to run %s (error)" % dispcmd)
            print("stdout was %s" % stdout)
        return None, p.returncode
    return stdout, p.returncode



@register_vcs_handler("git", "get_keywords")
def git_get_keywords(versionfile_abs):
    """Extract version information from the given file."""
    # the code embedded in _version.py can just fetch the value of these
    # keywords. When used from setup.py, we don't want to import _version.py,
    # so we do it with a regexp instead. This function is not used from
    # _version.py.
    keywords = {}
    try:
        f = open(versionfile_abs, "r")
        for line in f.readlines():
            if line.strip().startswith("git_refnames ="):
                mo = re.search(r'=\s*"(.*)"', line)
                if mo:
                    keywords["refnames"] = mo.group(1)
            if line.strip().startswith("git_full ="):
                mo = re.search(r'=\s*"(.*)"', line)
                if mo:
                    keywords["full"] = mo.group(1)
            if line.strip().startswith("git_date ="):
                mo = re.search(r'=\s*"(.*)"', line)
                if mo:
                    keywords["date"] = mo.group(1)
        f.close()
    except EnvironmentError:
        pass
    return keywords


@register_vcs_handler("git", "keywords")
def git_versions_from_keywords(keywords, tag_prefix, verbose):
    """Get version information from git keywords."""
    if not keywords:
        raise NotThisMethod("no keywords at all, weird")
    date = keywords.get("date")
    if date is not None:
        # Use only the last line.  Previous lines may contain GPG signature
        # information.
        date = date.splitlines()[-1]

        # git-2.2.0 added "%cI", which expands to an ISO-8601 -compliant
        # datestamp. However we prefer "%ci" (which expands to an "ISO-8601
        # -like" string, which we must then edit to make compliant), because
        # it's been around since git-1.5.3, and it's too difficult to
        # discover which version we're using, or to work around using an
        # older one.
        date = date.strip().replace(" ", "T", 1).replace(" ", "", 1)
    refnames = keywords["refnames"].strip()
    if refnames.startswith("$Format"):
        if verbose:
            print("keywords are unexpanded, not using")
        raise NotThisMethod("unexpanded keywords, not a git-archive tarball")
    refs = set([r.strip() for r in refnames.strip("()").split(",")])
    # starting in git-1.8.3, tags are listed as "tag: foo-1.0" instead of
    # just "foo-1.0". If we see a "tag: " prefix, prefer those.
    TAG = "tag: "
    tags = set([r[len(TAG):] for r in refs if r.startswith(TAG)])
    if not tags:
        # Either we're using git < 1.8.3, or there really are no tags. We use
        # a heuristic: assume all version tags have a digit. The old git %d
        # expansion behaves like git log --decorate=short and strips out the
        # refs/heads/ and refs/tags/ prefixes that would let us distinguish
        # between branches and tags. By ignoring refnames without digits, we
        # filter out many common branch names like "release" and
        # "stabilization", as well as "HEAD" and "master".
        tags = set([r for r in refs if re.search(r'\d', r)])
        if verbose:
            print("discarding '%s', no digits" % ",".join(refs - tags))
    if verbose:
        print("likely tags: %s" % ",".join(sorted(tags)))
    for ref in sorted(tags):
        # sorting will prefer e.g. "2.0" over "2.0rc1"
        if ref.startswith(tag_prefix):
            r = ref[len(tag_prefix):]
            if verbose:
                print("picking %s" % r)
            return {
                "version": r,
                "full-revisionid": keywords["full"].strip(),
                "dirty": False,
                "error": None,
                "date": date,
            }
    # no suitable tags, so version is "0+unknown", but full hex is still there
    if verbose:
        print("no suitable tags, using unknown + full revision id")
    return {
        "version": "0+unknown",
        "full-revisionid": keywords["full"].strip(),
        "dirty": False,
        "error": "no suitable tags",
        "date": None,
    }


@register_vcs_handler("git", "pieces_from_vcs")
def git_pieces_from_vcs(tag_prefix, root, verbose, run_command=run_command):
    """Get version from 'git describe' in the root of the source tree.

    This only gets called if the git-archive 'subst' keywords were *not*
    expanded, and _version.py hasn't already been rewritten with a short
    version string, meaning we're inside a checked out source tree.
    """
    GITS = ["git"]
    if sys.platform == "win32":
        GITS = ["git.cmd", "git.exe"]

    out, rc = run_command(GITS, ["rev-parse", "--git-dir"], cwd=root, hide_stderr=True)
    if rc != 0:
        if verbose:
            print("Directory %s not under git control" % root)
        raise NotThisMethod("'git rev-parse --git-dir' returned error")

    # if there is a tag matching tag_prefix, this yields TAG-NUM-gHEX[-dirty]
    # if there isn't one, this yields HEX[-dirty] (no NUM)
    describe_out, rc = run_command(
        GITS, ["describe", "--tags", "--dirty", "--always", "--long", "--match", "%s*" % tag_prefix], cwd=root
    )
    # --long was added in git-1.5.5
    if describe_out is None:
        raise NotThisMethod("'git describe' failed")
    describe_out = describe_out.strip()
    full_out, rc = run_command(GITS, ["rev-parse", "HEAD"], cwd=root)
    if full_out is None:
        raise NotThisMethod("'git rev-parse' failed")
    full_out = full_out.strip()

    pieces = {}
    pieces["long"] = full_out
    pieces["short"] = full_out[:7]  # maybe improved later
    pieces["error"] = None

    # parse describe_out. It will be like TAG-NUM-gHEX[-dirty] or HEX[-dirty]
    # TAG might have hyphens.
    git_describe = describe_out

    # look for -dirty suffix
    dirty = git_describe.endswith("-dirty")
    pieces["dirty"] = dirty
    if dirty:
        git_describe = git_describe[: git_describe.rindex("-dirty")]

    # now we have TAG-NUM-gHEX or HEX

    if "-" in git_describe:
        # TAG-NUM-gHEX
        mo = re.search(r'^(.+)-(\d+)-g([0-9a-f]+)$', git_describe)
        if not mo:
            # unparseable. Maybe git-describe is misbehaving?
            pieces["error"] = "unable to parse git-describe output: '%s'" % describe_out
            return pieces

        # tag
        full_tag = mo.group(1)
        if not full_tag.startswith(tag_prefix):
            if verbose:
                fmt = "tag '%s' doesn't start with prefix '%s'"
                print(fmt % (full_tag, tag_prefix))
            pieces["error"] = "tag '%s' doesn't start with prefix '%s'" % (full_tag, tag_prefix)
            return pieces
        pieces["closest-tag"] = full_tag[len(tag_prefix):]

        # distance: number of commits since tag
        pieces["distance"] = int(mo.group(2))

        # commit: short hex revision ID
        pieces["short"] = mo.group(3)

    else:
        # HEX: no tags
        pieces["closest-tag"] = None
        count_out, rc = run_command(GITS, ["rev-list", "HEAD", "--count"], cwd=root)
        pieces["distance"] = int(count_out)  # total number of commits

    # commit date: see ISO-8601 comment in git_versions_from_keywords()
    date = run_command(GITS, ["show", "-s", "--format=%ci", "HEAD"], cwd=root)[0].strip()
    # Use only the last line.  Previous lines may contain GPG signature
    # information.
    date = date.splitlines()[-1]
    pieces["date"] = date.strip().replace(" ", "T", 1).replace(" ", "", 1)
    pieces["run_number"] = int(
        os.environ.get("GITHUB_RUN_NUMBER", os.environ.get("CI_PIPELINE_IID", pieces["distance"]))
    )

    return pieces


def versions_from_parentdir(parentdir_prefix, root, verbose):
    """Try to determine the version from the parent directory name.

    Source tarballs conventionally unpack into a directory that includes both
    the project name and a version string. We will also support searching up
    two directory levels for an appropriately named parent directory
    """
    rootdirs = []

    for i in range(3):
        dirname = os.path.basename(root)
        if dirname.startswith(parentdir_prefix):
            return {
                "version": dirname[len(parentdir_prefix):],
                "full-revisionid": None,
                "dirty": False,
                "error": None,
                "date": None,
            }
        else:
            rootdirs.append(root)
            root = os.path.dirname(root)  # up a level

    if verbose:
        print("Tried directories %s but none started with prefix %s" % (str(rootdirs), parentdir_prefix))
    raise NotThisMethod("rootdir doesn't start with parentdir_prefix")


SHORT_VERSION_PY = """
# This file was generated by 'versioneer.py' (0.19) from
# revision-control system data, or from the parent directory name of an
# unpacked source archive. Distribution tarballs contain a pre-generated copy
# of this file.

import json

version_json = '''
%s
'''  # END VERSION_JSON


def get_versions():
    return json.loads(version_json)
"""


def versions_from_file(filename):
    """Try to determine the version from _version.py if present."""
    try:
        with open(filename) as f:
            contents = f.read()
    except EnvironmentError:
        raise NotThisMethod("unable to read _version.py")
    mo = re.search(r"version_json = '''\n(.*)'''  # END VERSION_JSON", contents, re.M | re.S)
    if not mo:
        mo = re.search(r"version_json = '''\r\n(.*)'''  # END VERSION_JSON", contents, re.M | re.S)
    if not mo:
        raise NotThisMethod("no version_json in _version.py")
    return json.loads(mo.group(1))


def write_to_version_file(filename, versions):
    """Write the given version number to the given _version.py file."""
    os.unlink(filename)
    contents = json.dumps(versions, sort_keys=True, indent=1, separators=(",", ": "))
    with open(filename, "w") as f:
        f.write(SHORT_VERSION_PY % contents)

    print("set %s to '%s'" % (filename, versions["version"]))


def plus_or_dot(pieces):
    """Return a + if we don't already have one, else return a ."""
    if "+" in pieces.get("closest-tag", ""):
        return "."
    return "+"


def render_pep440(pieces):
    """Build up version string, with post-release "local version identifier".

    Our goal: TAG[+DISTANCE.gHEX[.dirty]] . Note that if you
    get a tagged build and then dirty it, you'll get TAG+0.gHEX.dirty

    Exceptions:
    1: no tags. git_describe was just HEX. 0+untagged.DISTANCE.gHEX[.dirty]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            rendered += plus_or_dot(pieces)
            rendered += "%d.g%s" % (pieces["distance"], pieces["short"])
            if pieces["dirty"]:
                rendered += ".dirty"
    else:
        # exception #1
        rendered = "0+untagged.%d.g%s" % (pieces["distance"], pieces["short"])
        if pieces["dirty"]:
            rendered += ".dirty"
    return rendered


def render_pep440_pre(pieces):
    """TAG[.post2.devRUN_NUMBER] -- No -dirty.

    Exceptions:
    1: no tags. 0.post2.devRUN_NUMBER
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        # this needs to stay distance, run_number is always non-zero
        if pieces["distance"]:
            rendered += ".post2.dev%d" % pieces["run_number"]
    else:
        # exception #1
        rendered = "0.post2.dev%d" % pieces["distance"]
    return rendered


def render_pep440_post(pieces):
    """TAG[.postDISTANCE[.dev0]+gHEX] .

    The ".dev0" means dirty. Note that .dev0 sorts backwards
    (a dirty tree will appear "older" than the corresponding clean one),
    but you shouldn't be releasing software with -dirty anyways.

    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            rendered += ".post%d" % pieces["distance"]
            if pieces["dirty"]:
                rendered += ".dev0"
            rendered += plus_or_dot(pieces)
            rendered += "g%s" % pieces["short"]
    else:
        # exception #1
        rendered = "0.post%d" % pieces["distance"]
        if pieces["dirty"]:
            rendered += ".dev0"
        rendered += "+g%s" % pieces["short"]
    return rendered


def render_pep440_old(pieces):
    """TAG[.postDISTANCE[.dev0]] .

    The ".dev0" means dirty.

    Exceptions:
    1: no tags. 0.postDISTANCE[.dev0]
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"] or pieces["dirty"]:
            rendered += ".post%d" % pieces["distance"]
            if pieces["dirty"]:
                rendered += ".dev0"
    else:
        # exception #1
        rendered = "0.post%d" % pieces["distance"]
        if pieces["dirty"]:
            rendered += ".dev0"
    return rendered


def render_git_describe(pieces):
    """TAG[-DISTANCE-gHEX][-dirty].

    Like 'git describe --tags --dirty --always'.

    Exceptions:
    1: no tags. HEX[-dirty]  (note: no 'g' prefix)
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        if pieces["distance"]:
            rendered += "-%d-g%s" % (pieces["distance"], pieces["short"])
    else:
        # exception #1
        rendered = pieces["short"]
    if pieces["dirty"]:
        rendered += "-dirty"
    return rendered


def render_git_describe_long(pieces):
    """TAG-DISTANCE-gHEX[-dirty].

    Like 'git describe --tags --dirty --always -long'.
    The distance/hash is unconditional.

    Exceptions:
    1: no tags. HEX[-dirty]  (note: no 'g' prefix)
    """
    if pieces["closest-tag"]:
        rendered = pieces["closest-tag"]
        rendered += "-%d-g%s" % (pieces["distance"], pieces["short"])
    else:
        # exception #1
        rendered = pieces["short"]
    if pieces["dirty"]:
        rendered += "-dirty"
    return rendered


def render(pieces, style):
    """Render the given version pieces into the requested style."""
    if pieces["error"]:
        return {
            "version": "unknown",
            "full-revisionid": pieces.get("long"),
            "dirty": None,
            "error": pieces["error"],
            "date": None,
        }

    if not style or style == "default":
        style = "pep440"  # the default

    if style == "pep440":
        rendered = render_pep440(pieces)
    elif style == "pep440-pre":
        rendered = render_pep440_pre(pieces)
    elif style == "pep440-post":
        rendered = render_pep440_post(pieces)
    elif style == "pep440-old":
        rendered = render_pep440_old(pieces)
    elif style == "git-describe":
        rendered = render_git_describe(pieces)
    elif style == "git-describe-long":
        rendered = render_git_describe_long(pieces)
    else:
        raise ValueError("unknown style '%s'" % style)

    return {
        "version": rendered,
        "full-revisionid": pieces["long"],
        "dirty": pieces["dirty"],
        "error": None,
        "date": pieces.get("date"),
    }


class VersioneerBadRootError(Exception):
    """The project root directory is unknown or missing key files."""


def get_versions(verbose=False):
    """Get the project version from whatever source is available.

    Returns dict with two keys: 'version' and 'full'.
    """
    if "versioneer" in sys.modules:
        # see the discussion in cmdclass.py:get_cmdclass()
        del sys.modules["versioneer"]

    root = get_root()
    cfg = get_config_from_root(root)

    assert cfg.VCS is not None, "please set [versioneer]VCS= in setup.cfg"
    handlers = HANDLERS.get(cfg.VCS)
    assert handlers, "unrecognized VCS '%s'" % cfg.VCS
    verbose = verbose or cfg.verbose
    assert cfg.versionfile_source is not None, "please set versioneer.versionfile_source"
    assert cfg.tag_prefix is not None, "please set versioneer.tag_prefix"

    versionfile_abs = os.path.join(root, cfg.versionfile_source)

    # extract version from first of: _version.py, VCS command (e.g. 'git
    # describe'), parentdir. This is meant to work for developers using a
    # source checkout, for users of a tarball created by 'setup.py sdist',
    # and for users of a tarball/zipball created by 'git archive' or github's
    # download-from-tag feature or the equivalent in other VCSes.

    get_keywords_f = handlers.get("get_keywords")
    from_keywords_f = handlers.get("keywords")
    if get_keywords_f and from_keywords_f:
        try:
            keywords = get_keywords_f(versionfile_abs)
            ver = from_keywords_f(keywords, cfg.tag_prefix, verbose)
            if verbose:
                print("got version from expanded keyword %s" % ver)
            return ver
        except NotThisMethod:
            pass

    try:
        ver = versions_from_file(versionfile_abs)
        if verbose:
            print("got version from file %s %s" % (versionfile_abs, ver))
        return ver
    except NotThisMethod:
        pass

    from_vcs_f = handlers.get("pieces_from_vcs")
    if from_vcs_f:
        try:
            pieces = from_vcs_f(cfg.tag_prefix, root, verbose)
            ver = render(pieces, cfg.style)
            if verbose:
                print("got version from VCS %s" % ver)
            return ver
        except NotThisMethod:
            pass

    try:
        if cfg.parentdir_prefix:
            ver = versions_from_parentdir(cfg.parentdir_prefix, root, verbose)
            if verbose:
                print("got version from parentdir %s" % ver)
            return ver
    except NotThisMethod:
        pass

    if verbose:
        print("unable to compute version")

    return {
        "version": "0+unknown",
        "full-revisionid": None,
        "dirty": None,
        "error": "unable to compute version",
        "date": None,
    }


def get_version():
    """Get the short version string for this project."""
    return get_versions()["version"]

if __name__ == "__main__":
    PROJECT_ROOT = sys.argv[1]
    print(get_version())
