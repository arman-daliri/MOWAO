import os
from setuptools.command.build import build as _build
from setuptools import Extension, setup

NAME = "mowao"


class build(_build):
    def finalize_options(self):
        super().finalize_options()

        if os.name == "nt":
            self.compiler = "mingw32"  # msvc, cygwin
        elif os.name == "posix":
            self.compiler = "unix"


def files_get(directory):
    files = os.listdir(directory)
    f = []
    for i in files:
        if i.endswith(".c"):
            f.append(directory + i)
    return f


def main():
    mowao_src = "mowao/mowao.c"
    sources = files_get("src/")
    sources.append(mowao_src)
    module = Extension(
        NAME, sources=sources,
        extra_compile_args=["-O3", "-std=c99"],
        extra_link_args=["-lm"],
    )

    setup(
        ext_modules=[module],
        cmdclass={"build": build},
    )


if __name__ == "__main__":
    main()
